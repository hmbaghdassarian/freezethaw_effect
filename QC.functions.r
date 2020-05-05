
##### Counts QC


# induce p percent (of the sum of reads in each sample) random reads into each sample
induced_noise_controls <- function(dat,p=c(.1,1,10),assumption='uniform',max_count=1e2,n=round(nrow(dat)*.5),scale=1){
    dat_orig = dat
    sample_names = colnames(dat_orig)
    colnames(dat) = paste0('0_',colnames(dat),'_0')
    dat_out = list(dat)
    ic = 1
    for(i in p){
        ic = ic+1
        noise_ij = c()
        iter_ij=c()
        dat = dat_orig
        for(col in 1:ncol(dat)){
            if(assumption=='mean_weighted'){
                sel = table(sample(1:nrow(dat),i*nrow(dat),replace = TRUE,prob=prob)) # randomly add i% of the total reads to genes 
                dat[as.numeric(names(sel)),col] = dat[as.numeric(names(sel)),col] + sel
            }else if(assumption=='uniform'){
                sel = table(sample(1:nrow(dat),i*nrow(dat),replace = TRUE)) # randomly add i% of the total reads to genes
                dat[as.numeric(names(sel)),col] = dat[as.numeric(names(sel)),col] + sel
            }else if(assumption=='poisson'){
                orig=dat[,col]
                sum_orig=sum(orig)
                ith=orig
                divrg=0
                ic2=0
                while( (divrg < i) & (ic2<max_count)){
                    ic2=ic2+1
                    choose=sample(1:length(ith),n)
                    ### instead, select poisson mean from gamma distribution parameterized on dispersion (only one iteration)
                    ### negative-binomial <- poisson(lambda=gamma(scale=(set to match data),shape=mean/scale=expression[i,j]/scale))
                    mu = ith[choose] # primary mu estimation
                    #sigma = disp[choose]
                    #theta= sigma/mu 
                    #k=(mu^2)/theta
                    #mu = round(rgamma(n,shape = k , scale=theta)) # gamma updated mu estimation
                    ith[choose] = rpois(n,lambda=mu)
                    divrg = sum(abs(ith-orig))/sum(orig)
                }
                dat[,col]=ith
                noise_ij = c(noise_ij,signif(divrg,5))
                iter_ij = c(iter_ij,ic2)
            }else if(assumption=='negbinomial'){
                orig=dat[,col]
                sum_orig=sum(orig)
                ith=orig
                divrg=0
                ic2=0
                df_stat = data.frame(log_var=log(apply(dat,1,var)),log_mean=log(apply(dat,1,mean)+1)) ###
                los = loess( log_var ~ log_mean , data=na.omit(df_stat[is.finite(df_stat$log_var),]) ) ###
                while( (divrg < i) & (ic2<max_count)){
                    ic2=ic2+1
                    #rint(ith[1:10])
                    disp = exp(predict(object=los,newdata = log(ith+1))) ####
                    jth = rnbinom(length(ith),mu=ith,size=disp)
                    jth[is.na(jth)] = ith[is.na(jth)]
                    ith=jth
                    #print(ith[1:10])
                    divrg = sum(abs(ith-orig))/sum_orig
                }
                dat[,col]=ith
                noise_ij = c(noise_ij,signif(divrg,5))
                iter_ij = c(iter_ij,ic2)
            }else if(assumption=='negbinomial_indirect'){
                orig=dat[,col]
                sum_orig=sum(orig)
                ith=orig
                divrg=0
                ic2=0
                disp=apply(dat,1,var)
                while( (divrg < i) & (ic2<max_count)){
                    ic2=ic2+1
                    choose=sample(1:length(ith),n)
                    ### instead, select poisson mean from gamma distribution parameterized on dispersion (only one iteration)
                    ### negative-binomial <- poisson(lambda=gamma(scale=(set to match data),shape=mean/scale=expression[i,j]/scale))
                    mu = ith[choose] # primary mu estimation
                    sigma = disp[choose]
                    theta= sigma/mu 
                    k=(mu^2)/theta
                    mu = round(rgamma(n,shape = k , scale=theta)) # gamma updated mu estimation
                    ith[choose] = rpois(n,lambda=mu)
                    divrg = sum(abs(ith-orig))/sum(orig)
                }
                dat[,col]=ith
                noise_ij = c(noise_ij,signif(divrg,5))
                iter_ij = c(iter_ij,ic2)
            }else if(assumption=='something_else'){
                orig=dat[,col]
                sum_orig=sum(orig)
                ith=orig
                divrg=0
                ic2=0
                disp=apply(dat,1,var)
                while( (divrg < i) & (ic2<max_count)){
                    ic2=ic2+1
                    choose=sample(1:length(ith),n)
                    ### instead, select poisson mean from gamma distribution parameterized on dispersion (only one iteration)
                    ### negative-binomial <- poisson(lambda=gamma(scale=(set to match data),shape=mean/scale=expression[i,j]/scale))
                    mu = ith[choose] # primary mu estimation
                    sigma = disp[choose]
                    theta= sigma/mu 
                    #k=(mu^2)/theta
                    k=mu
                    mu = round(rgamma(n,shape = k , scale=theta)) # gamma updated mu estimation
                    ith[choose] = rpois(n,lambda=mu)
                    divrg = sum(abs(ith-orig))/sum(orig)
                }
                dat[,col]=ith
                noise_ij = c(noise_ij,signif(divrg,5))
                iter_ij = c(iter_ij,ic2)
            }
        }
        colnames(dat) = paste(as.character(noise_ij),colnames(dat),as.character(iter_ij),sep='_')
        dat_out[[ic]] = dat
    }
    out = do.call(cbind,dat_out)
    out
}
                
#### distribution of differences

# function
reshape_diff<-function(dat_in,dat_i,i){
    # create matric of original and replicate sampels
    tmp = cbind(dat_in[,i],dat_i)
    colnames(tmp)[1] = colnames(dat_in)[i] # name original appropriately
    m = melt(tmp,id=colnames(dat_in)[i]) # melt while retaining original sample as reference (id)
    colnames(m)[1:3] = c('original_sample','sample_name','original_minus_replicate')
    tmp2 = do.call(rbind,strsplit(as.character(m$sample_name),'_')) #add precent noise and parent sample information for replicates
    m = cbind(m,
              percent_noise=as.numeric(as.character(tmp2[,1]))*100,
              parent_sample=tmp2[,2],
              original_sample_id=samples_in[i]
             )
    m$percent_noise = as.factor(round(m$percent_noise))
    m$sample_type=ifelse(m$percent_noise==0,'technical_replicate','noise_replicate') # distinguish between noise and technical replicates
    m
}

# Comparison Functions

### root mean square error
rmse <- function(x,y) (mean((x-y)^2))^(1/2)
    
### remove diagonal
rm_diag <- function(x){
    x[!(upper.tri(x)|lower.tri(x))]=NA
    x
}
    
# poisson based comparisons ask the probability of sample x given the poisson distributions described by sample y
### poisson probability 
pois_p <- function(x,y,ps=1e-10){
    p = dpois(x,lambda=y)
    prod(p[p>0])
}
    
### poisson probability 
log_lik_pois_p <- function(x,y,ps=1e-10){
    p = dpois(x,lambda=y)
    -sum(log(p[p>0]))
}
    
### similarity matrix reshape functions
similarity<-function(dat,method,rep,plot_on=F,conf=c('final_freeze_thaw','concentration_bins','RINe'),unique_sample_name='Comment',rep_sample_name='SampleName'){
    #print(head(dat))
    if(method=='euclidean'){
        out=as.matrix(dist(t(dat)))
    }else if(method=='pearson'){
        out=cor(dat,method=method)
    }else if(method=='spearman'){
        out=cor(dat,method=method)
    }else if(method=='rmse'){
        out=apply(dat,2,function(x){ apply(dat,2,function(y) rmse(x,y)) }) 
    }else if(method=='poissonProb'){
        out=apply(dat,2,function(x){ apply(dat,2,function(y) log_lik_pois_p(x,y)) })
    }else if(method=='jaccard'){
        out=as.matrix(dist(t(dat),method='binary'))
    }else{stop('method improperly specified')}
        #print(1)          
    out = out/max(out)
    #print(head(out))
    out=simMx_to_long(out,metric=method,rep=rep,conf=conf,unique_sample_name=unique_sample_name,rep_sample_name=rep_sample_name)
    #print(head(out))
    out
}
                           
simMx_to_long <- function(x,metric,rep,conf,unique_sample_name='Comment',rep_sample_name='SampleName'){
    # extract 0% noise induced controls
    x = x[,grepl('^0_',colnames(x))]
    if(prod(dim(x))==0){stop('x must inlcude 0% noise induced controls')}
    if(dim(x)[1]==dim(x)[2]){stop('NICs were not properly removed')}
    # melt
    tmp=unique(na.omit(melt(x,direction='long')))
    colnames(tmp) = c('induced_noise_controls_complete','original_sample_complete','Similarity')
    tmp_i = do.call(rbind,strsplit(as.character(tmp$induced_noise_controls_complete),'_'))
    
    # extract NICs with percent noise and the parent sample (from which they were derived)
    tmp$Percent_Noise = as.numeric(as.character(tmp_i[,1]))
    tmp$Parent_Sample = tmp_i[,2]
    
    tmp$Original_Sample = do.call(rbind,strsplit(as.character(tmp$original_sample_complete),'_'))[,2]
   
    if(length(rep)==1 && !is.null(rep)){
        tmp$technical_replicate = rep
    }
    
    tmp2 = tmp
    #colnames(tmp2) = c('controls','samples','Similarity','Percent_Noise','Parent_Sample','Original_Sample','technical_replicate')
    if(!is.null(conf)){
        tmp2 = merge(tmp2,annot[annot[[rep_sample_name]]==rep,c(conf,unique_sample_name)],by.x='Parent_Sample',by.y=unique_sample_name,all.x=T,all.y=F)
    }
    tmp2$Percent_Noise = as.numeric(as.character(tmp2$Percent_Noise))
    tmp2$Parent_Sample = factor(tmp2$Parent_Sample,levels=sort(unique(tmp2$Parent_Sample)))
    tmp2$self_or_replicate = ifelse(as.character(tmp2$Parent_Sample)==as.character(tmp2$Original_Sample),'Original','Replicate')
    tmp2$Similarity_Metric = factor( rep(metric,nrow(tmp2)) , levels=c('euclidean','rmse','pearson','spearman','jaccard','poissonProb'))
    # minimum replicate value (zero noise)
    tmp2$first_replicate = sapply(tmp2$Original_Sample,function(x){
        median( tmp2$Similarity[tmp2$Original_Sample==x & tmp2$self_or_replicate=='Replicate' & tmp2$Percent_Noise==0] )
#        mean( tmp2$Similarity[tmp2$Original_Sample==x & tmp2$self_or_replicate=='Replicate' & tmp2$Percent_Noise==0] )
    })
    # what is the smallest value greater than the minumum replicate value
    tmp2$min_cross_sim = sapply(tmp2$Original_Sample,function(x){
        a=tmp2$Similarity[tmp2$Original_Sample==x & tmp2$self_or_replicate=='Original']
        b=tmp2$first_replicate[tmp2$Original_Sample==x & tmp2$self_or_replicate=='Original']
        a[which.min(abs(a-b))]
    })
    tmp2$min_cross_noise = sapply(tmp2$Original_Sample,function(x){
        a=tmp2$Similarity[tmp2$Original_Sample==x & tmp2$self_or_replicate=='Original']
        b=tmp2$first_replicate[tmp2$Original_Sample==x & tmp2$self_or_replicate=='Original']
        c=tmp2$Percent_Noise[tmp2$Original_Sample==x & tmp2$self_or_replicate=='Original']
        c[which.min(abs(a-b))]
    })
    tmp2
}
                  
RPKM <- function(count_i,gene_length) (count_i/((sum(count_i)/1e6)*(gene_length/1e3)))


test_NIC <- function(nrow=100,dat=NULL,p=c(.1,.5,1,2,2.5)*.01){
    if(is.null(dat)){ 
        rpi <- function(nrow){ rpois(nrow,sample(1:(nrow*100),nrow)) }
        dat=data.frame(a=rpi(nrow),b=rpi(nrow),c=rpi(nrow))
    }
    dat_out = induced_noise_controls(dat,p=p,assumption='poisson',max_count=1e5)

    tmp_i = as.data.frame( do.call(rbind,strsplit(colnames(dat_out),'_')) )
    tmp_i$iterations = as.numeric(as.character(tmp_i[,3]))
    tmp_i$percent_noise = as.numeric(as.character(tmp_i[,1]))
    tmp_i$col =tmp_i[,2]
    tmp_i
}

################### DE QC

aggregate_DE_sim <- function(p,res_out_all,basic_vis=T){
  ### COLLECT SIMILARITY METRICS
  
  sim = list()
  for(i in as.character(unique(p))){
    print(i)
    #    for(j in names(res_out_all)[1:5]){
    row_names = sapply(names(res_out_all),function(x) ifelse(is.null(res_out_all[[x]][[i]]),NA,x) )
    tmp1 = t(do.call(cbind,lapply(names(res_out_all),function(x) res_out_all[[x]][[i]]$log2FoldChange)) )
    tmp2 = t(do.call(cbind,lapply(names(res_out_all),function(x) res_out_all[[x]][[i]]$padj)) )
    
    tmp1[is.na(tmp1)] = 0
    tmp2[is.na(tmp2)] = 1
    tmp2 = ifelse(tmp2<.1,1,0)
    
    tmp1 = tmp1[,colSums(tmp2,na.rm = T)>0]
    tmp2 = tmp2[,colSums(tmp2,na.rm = T)>0]
    
    rownames(tmp1) = na.omit(row_names)
    rownames(tmp2) = na.omit(row_names)
    
    if(prod(dim(tmp1))==0){next}
    
    #        mx1=as.matrix(dist(tmp1))
    mx1=cor(t(tmp1))
    mx1[lower.tri(mx1) ] = NA
    v1 = na.omit(melt(mx1))
    mx2=as.matrix(1-dist(tmp2,method='binary'))
    #mx2=apply(tmp2,1,function(x){apply(tmp2,1,function(y){
    #    if(sum(x)+sum(y)==0){return(0)}
    #    length(intersect(which(x==1),which(y==1)))/length(union(which(x==1),which(y==1)))})}) # get union size
    mx2[lower.tri(mx2)]=NA
    v2 = na.omit(melt(mx2))
    mx3=apply(tmp2,1,function(x){apply(tmp2,1,function(y){length(intersect(which(x==1),which(y==1)))})}) # get union size
    mx3[lower.tri(mx2)]=NA
    v3 = na.omit(melt(mx3))
    
    # HB new
    print('QC.functions.r: Rename columns start')
    colnames(v1)[3]='v1'
    colnames(v2)[3]='v2'
    colnames(v3)[3]='v3'
    print('QC.functions.r: Rename columns end')
                                  
#     # BK: old
#     comb_i =  merge(cbind(merge(v1,v2,by=c('X1','X2')),as.numeric(i)),v3)
# #     HB: new
    print('QC.functions.r: Var1/Var2 change part 1')
    comb_i =  merge(cbind(merge(v1,v2),as.numeric(i)),v3)
    print('QC.functions.r: Var1/Var2 change part 2')
    colnames(comb_i) = c('cbx1','cbx2','lfc_sim','de_sim','sample_portion','intersect_size')
    sim[[i]] = comb_i
    #        sim[[i]] = data.frame(lfc_sim=v1,de_sim=v2,de_count=mean(rowSums(tmp2)),sample_portion=as.numeric(i)) #,sample_portion=i,combination=j)
    #    }
  }
  sim_out = do.call(rbind,sim)
  #        head(sim_out)
  
  g1 = ggplot(sim_out, aes(x=as.character(sample_portion),y=lfc_sim)) + geom_boxplot() 
  #g2 = ggplot(sim_out[sim_out$intersect_size>0,], aes(x=as.character(sample_portion),y=de_sim,
  #            color=log(intersect_size),alpha=log(intersect_size))) + geom_boxplot() + geom_jitter()
  #                                      grid.arrange(g1,g2,ncol=1)
  if(basic_vis){print( g1+ylim(c(0,1)) )}
  return(sim_out)
}

integrate_MCN <- function(sample_sim,non_replicate_combinations,annot,sim_out){
  tmp = sample_sim
  
  cbx_MCN1 = unlist(lapply(non_replicate_combinations,function(x){
    quantile( sapply(x,function(y){  mean(tmp$min_cross_noise[tmp$Original_Sample==y],na.rm=T) }),na.rm=T)[1] }))
  cbx_MCN2 = unlist(lapply(non_replicate_combinations,function(x){
    quantile( sapply(x,function(y){  mean(tmp$min_cross_noise[tmp$Original_Sample==y],na.rm=T) }),na.rm=T)[2] }))
  cbx_MCN3 = unlist(lapply(non_replicate_combinations,function(x){
    quantile( sapply(x,function(y){  mean(tmp$min_cross_noise[tmp$Original_Sample==y],na.rm=T) }),na.rm=T)[3] }))
  cbx_MCN4 = unlist(lapply(non_replicate_combinations,function(x){
    quantile( sapply(x,function(y){  mean(tmp$min_cross_noise[tmp$Original_Sample==y],na.rm=T) }),na.rm=T)[4] }))
  cbx_MCN5 = unlist(lapply(non_replicate_combinations,function(x){
    quantile( sapply(x,function(y){  mean(tmp$min_cross_noise[tmp$Original_Sample==y],na.rm=T) }),na.rm=T)[5] }))
  cbx_MCN_mean = unlist(lapply(non_replicate_combinations,function(x){
    mean( sapply(x,function(y){  mean(tmp$min_cross_noise[tmp$Original_Sample==y],na.rm=T) }),na.rm=T) }))
  cbx_MCN_sum = unlist(lapply(non_replicate_combinations,function(x){
    sum( sapply(x,function(y){  mean(tmp$min_cross_noise[tmp$Original_Sample==y],na.rm=T) }),na.rm=T) }))
  
  # make error annotation data frame
  df1 = data.frame(
    cbx_freeze_thaw_mean = unlist(lapply(non_replicate_combinations,function(x){ mean(annot$final_freeze_thaw[annot$Comment%in%x]) })),
    #    cbx_concentration_mean = unlist(lapply(non_replicate_combinations,function(x){ mean(annot$concentration[annot$Comment%in%x]) })),
    cbx_RIN_mean = unlist(lapply(non_replicate_combinations,function(x){ mean(annot$RINe[annot$Comment%in%x]) })),
    #cbx_conc_mean = unlist(lapply(non_replicate_combinations,function(x){ mean(annot$concentration[annot$Comment%in%x]) })),
    cbx_MCN1=cbx_MCN1, cbx_MCN2=cbx_MCN2, cbx_MCN3=cbx_MCN3, cbx_MCN4=cbx_MCN4, cbx_MCN5=cbx_MCN5,
    #cbx_MCN_mean=cbx_MCN_mean, #cbx_MCN_sum=cbx_MCN_sum,
    #cbx_MCNpred = unlist(lapply(non_replicate_combinations,function(x){ mean(tmp$min_cross_noise_PRED[tmp$Original_Sample%in%x]) })),
    cbx_size = unlist(lapply(non_replicate_combinations,length)),
    Original_Sample=paste0('C',1:length(non_replicate_combinations))
  )
  #df = merge(df1,sim_out,by.x='Original_Sample',by.y='cbx1') # bk
  df = NULL
  return(list(df1,df))
}

vis_reproducability_plus <- function(df,df1,term='cbx_RIN_mean',signif=3,bins=3){
  # LFC SIMILARITY from quality means
  keep = df$cbx2 %in% df1$Original_Sample[scale(df1$cbx_MCN5)<0 & scale(df1$cbx_MCN3)<0] # only compare to low-moderate noise samples

  df[[termCH<-paste0(term,'_char')]] = as.character(cut(signif(df[[term]],bins),signif))

  g1=ggplot(df[keep,], aes_string(y='lfc_sim',x='jitter(sample_portion)',color=termCH)) + geom_point(alpha=.1) + stat_smooth(method='glm')+ ylim(c(0,1))
  g1
}

######################
        #### noise metric: counts of noise / counts of signal
########################
