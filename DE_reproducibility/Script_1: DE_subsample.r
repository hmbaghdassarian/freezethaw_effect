cr = 15 # number of cores
# n *n_iter is the number of subsets to be generated per subset size
n = 1000 
n_iter = 2

library(reshape)
library(abind)
library(DESeq2)
library(svMisc)
library(openxlsx)
library(RUVSeq)
library(sva)
library(BiocParallel)
library(foreach)
library(doSNOW)

annot = read.csv('data/tmpsaveall_annot.csv', row.names = 1)
expss = read.csv('data/counts.csv', row.names=1)

par=T

use_cmp = TRUE
plot_on=TRUE

ii=0


print('Setup parameters')

res_out_name=c()
genes_sel = rowMeans(expss)>20 
res_out=list()
res_out_all=list()                    

                    
min = .375-.125
max = 0.875 
step = 0.125                  
p=rep(seq(min,max,step),each=n_iter)
n_p = length(unique(p))

sample_filter = ((abs(scale(annot$concentration))<1.96) & (!is.na(annot$DX)) & (!is.na(annot$RINe)) 
                & (!is.na(annot$final_freeze_thaw)))
annot_ = annot[which(sample_filter),] 
 
w_ASD = which(annot_$DX=='ASD')
w_TD = which(annot_$DX=='TD') 
sel = c(w_TD,w_ASD)
total_samples = length(sel)
sel1 = sel[!duplicated(annot_$SampleName[sel])] 
maximal_size_subset = length(sel1)

# sanity check
if (maximal_size_subset == length(annot$DX[!duplicated(annot$SampleName) & !is.na(annot$DX)])){
    print('Check 1 complete')
}

n_ASD = length(intersect(w_ASD, sel1))
n_TD = length(intersect(w_TD, sel1))
n_remove_variance = 0

non_replicate_sel_list = list()
for (i in (1:n)) {
    non_replicate_sel = c()
    while (length(non_replicate_sel) != maximal_size_subset - (2*n_remove_variance)) {
        
        # subsample in a manner to get a distribution about 1.5 average freeze thaws ('prob=' portion does this)
        nr_sel_ASD = sample(w_ASD,n_ASD-n_remove_variance,replace = F)
        nr_sel_TD = sample(w_TD,n_TD-n_remove_variance,replace = F)
        non_replicate_sel = c(nr_sel_ASD,nr_sel_TD)
        
        # no replicates of the same sample requirement
        non_replicate_sel = non_replicate_sel[!duplicated(annot$SampleName[non_replicate_sel])] 
        
    }
    non_replicate_sel_list[length(non_replicate_sel_list) + 1] = list(non_replicate_sel) # append to list
}

non_replicate_combinations = list()
for (i in seq(1:n)){
    val = annot_$Comment[non_replicate_sel_list[[i]]]
    non_replicate_combinations[length(non_replicate_combinations) + 1] = list(val)
    }

print('Drop duplicates and calculate actual proportions relative to valid set size')
dup=duplicated(lapply(non_replicate_combinations,function(x) paste(sort(x),collapse='_')))
non_replicate_combinations = non_replicate_combinations[dup == F] 
                      
names(non_replicate_combinations) = paste0('C',1:length(non_replicate_combinations))
lengths = c()
for (repN in names(non_replicate_combinations)){
    lengths = c(lengths,length(non_replicate_combinations[[repN]]))}

if (length(non_replicate_combinations) < n){ 
    stop('Not enough n')
}

if (min(lengths) == max(lengths)) {
    max_p = max(lengths)/total_samples
}
actual_p = p*max_p

print('Generate random subsets according to p')
label = list()
non_replicate_subsets = list()
for (j in seq(1,length(names(non_replicate_combinations)))){
    progress(j)
    repN = names(non_replicate_combinations)[[j]]
    
    rep=non_replicate_combinations[[repN]]
    if(length(unique(annot$DX[annot$Comment%in%rep]))<=1){list(NA,NA)}
    rep = annot_[which(annot_$Comment%in%rep),]
    TD = as.vector(rep[which(rep$DX == 'TD'),c('Comment')])
    ASD = as.vector(rep[which(rep$DX == 'ASD'),c('Comment')])
    
    res=list()
    for (i in p){
        rep_i = NULL
        
        while ((is.null(rep_i)) | (list(rep_i) %in% non_replicate_subsets)){ #ensures no duplicate subsamples
            TD_i = sample(TD,round(dim(rep)*i)/2, replace = F)
            ASD_i = sample(ASD,round(dim(rep)*i)/2, replace = F)
            rep_i = c(TD_i,ASD_i)}
#          
        non_replicate_subsets[[length(non_replicate_subsets)+1]] = rep_i
        label[[length(label)+1]] = paste0(unique(actual_p)[match(i,unique(p))],'_',repN)
    }
}
                      
# unit test
for (i in 1:length(non_replicate_subsets)){
    t = table(annot_[annot_$Comment%in%non_replicate_subsets[[i]],]$DX)
    if (t[1] != t[2]){
        stop('Subsampling error')
    }
}       
                      
write.csv(cbind(non_replicate_subsets),paste0("data/iteration_samples.n.",n,".n_iter.",n_iter,".csv"))
write.csv(cbind(label),paste0("data/iteration_labels.n.",n,".n_iter.",n_iter,".csv"))
lapply(actual_p, write, paste0("data/actual_proportions.n.",n,".n_iter.",n_iter,".txt"), append=TRUE, ncolumns=1000)
                    
             

# RUV
annot = read.csv('data/tmpsaveall_annot.csv', row.names = 1)
annot_ = annot[which(sample_filter),]
# section 2.4: http://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf                
x<-as.character(annot_[,c('DX')]) 
filtered<-expss[,rownames(annot_)]
set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(x, row.names=colnames(filtered)))
design = model.matrix(~DX,data=annot_)  # only on filtered data
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))] #negative controls
set2 <- RUVg(set, empirical, k=1)  

if (!identical(as.vector(annot_$Comment), rownames(annot_)) ){stop('RUV metadata processing issue')}  
if (!identical(rownames(design),rownames(annot_))){stop('RUV design matrix issue')}
if (!identical(rownames(pData(set2)),rownames(annot_))){stop('RUV data fitting issue 1')}
if (!identical(pData(set2)$x,annot_$DX)){stop('RUV data fitting issue')}
                      
write.table(pData(set2), file = 'data/RUV.covariate.tsv')
                      
cl <- makeCluster(cr)
registerDoSNOW(cl)
pb <- txtProgressBar(max = (n*n_iter*n_p), style = 3)
progress <- function(n_) setTxtProgressBar(pb, n_) 
opts <- list(progress = progress)  
print('Do parallel DEseq:RUV')
                
par_out = foreach(rep_i=non_replicate_subsets,.packages = c('DESeq2', 'RUVSeq'), .options.snow = opts) %dopar% {
    print(rep_i)
    dat_i = counts(set2)[genes_sel,as.vector(annot_[annot_$Comment%in%rep_i,c('Comment')])]
    mD = pData(set2)[annot_$Comment%in%rep_i,]
    if (!identical(colnames(dat_i), rownames(mD))){stop('RUV DEseq subset match error 1')} 
    if(!identical(sort(rep_i),sort(colnames(dat_i)))){stop('RUV DEseq subset match error 3')}
    if(!identical(sort(rep_i),sort(rownames(mD)))){stop('RUV DEseq subset match error 4')}                                           
    dds=NULL
    try(dds <- DESeqDataSetFromMatrix(countData = dat_i, 
                                      colData = mD,
                                      design = ~ W_1 + x ))
    res_i=NULL
#     disp_i=NULL
    try( dds <- DESeq(dds) )
    try( res_i <- results(dds) )
    try( disp_i <- mcols(dds)$dispFit )
    res<-list(res_i,disp_i)
}            
close(pb)
stopCluster(cl)
                    
index = rownames(par_out[[1]][[1]])

lfc_df = data.frame(matrix(ncol = length(par_out), nrow = length(index)))
padj_df = data.frame(matrix(ncol = length(par_out), nrow = length(index)))
base_mean_df = data.frame(matrix(ncol = length(par_out), nrow = length(index)))                      
lfc_SE_df = data.frame(matrix(ncol = length(par_out), nrow = length(index)))
lfc_dispersion_df = data.frame(matrix(ncol = length(par_out), nrow = length(index)))                

rownames(lfc_df) = index
rownames(padj_df) = index
rownames(base_mean_df) = index
rownames(lfc_SE_df) = index
rownames(lfc_dispersion_df) = index

for (i in seq(1,length(par_out))){
    lfc_df[,i] = par_out[[i]][[1]][['log2FoldChange']]
    padj_df[,i] = par_out[[i]][[1]][['padj']]
    base_mean_df[,i] = par_out[[i]][[1]][['baseMean']]
    lfc_SE_df[,i] = par_out[[i]][[1]][['lfcSE']]
    lfc_dispersion_df[,i] = par_out[[i]][[2]]
}

write.csv(lfc_df, paste0("data/RUV.lfc.n.",n,".n_iter.",n_iter,".csv"))
write.csv(padj_df, paste0("data/RUV.padj.n.",n,".n_iter.",n_iter,".csv"))
write.csv(base_mean_df, paste0("data/RUV.base_mean.n.",n,".n_iter.",n_iter,".csv"))
write.csv(lfc_SE_df, paste0("data/RUV.lfc.SE.n.",n,".n_iter.",n_iter,".csv"))         
write.csv(lfc_dispersion_df, paste0("data/RUV.lfc.dispersion.n.",n,".n_iter.",n_iter,".csv"))