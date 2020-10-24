#!/usr/bin/env python
# coding: utf-8

# In[1]:


print('Import packages')
import pandas as pd
import numpy as np

from scipy.special import binom
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

from itertools import combinations

import matplotlib.pyplot as plt
import seaborn as sns

from tqdm import trange
from tqdm import tqdm
from multiprocessing import Pool
import time
import random


# In[3]:


cr = 4 # no. of cores
n, n_iter, alpha, lfc_min_thershold, base_mean_min_threshold, min_percentage, max_cutoff, normalization_type = 1000, 2, 0.1, 0.1, 0.1, 0.0, 0.2, 'RUV'


# In[2]:


d = pd.read_csv('../data/rnaseq_pd.isaac.csv', index_col = 0)
manifest = pd.read_excel('../manifest.xlsx', sheet_name = 2)

# filter out industry standards and missing metadata
polyA = d[(d['extraction'] == 'polyA') & (d['DX'] != 'Clontech-2') & (d['RINe'].notna()) & (d['DX'].notna())
          & (d['final_freeze_thaw'].notna())]
print(polyA.shape[0])
polyA = [i.split('_')[0] for i in polyA.index]

rD = d[(d['extraction'] == 'riboDep') & (d['DX'] != 'Clontech-2') & (d['DX'].notna())
  & (d['RINe'].notna()) & (d['final_freeze_thaw'].notna())]
print(rD.shape[0])
rD = [i.split('_')[0] for i in rD.index]


# In[72]:


print('Parse subset labels')

string_ = "../data/iteration_labels.n."+ str(n) + '.n_iter.' + str(n_iter) + '.csv' 
  
subset_labels = open(string_).read().splitlines()[1:]
subset_labels = [i.split(',')[1] for i in subset_labels]
subset_labels = [i[0:4] + i[i.index('_'):] for i in subset_labels]

t = sorted(set(subset_labels))
count = dict(zip(t,[1]*len(t)))

for i in range(len(subset_labels)):
    val = subset_labels[i]
    subset_labels[i] = val + '_' + str(count[val])
    count[val] += 1

print('Parse subset samples')
string_ = "../data/iteration_samples.n."+ str(n) + '.n_iter.' + str(n_iter) + '.csv'

subset_samples = open(string_).read().splitlines()[1:]
subset_samples = [i.split('c(')[1][:-1].split(',') for i in subset_samples]

i = subset_samples[0]
for i in trange(len(subset_samples)):
    list_ = subset_samples[i]
    for j in range(len(list_)):
        val = list_[j]
        index = [k for k in range(len(val)) if val[k] == '"']
        list_[j] = val[index[0]+1:index[1]]
    subset_samples[i] = list_
subset_samples = [' '.join(i) for i in subset_samples]
    
if len(subset_samples) == len(subset_labels):
    colLabel_to_subset = dict(zip(subset_labels,subset_samples))
else:
    raise ValueError('Samples or labels were parsed incorrectly')

d___ = pd.DataFrame(data = {'colLabel': list(colLabel_to_subset.keys()), 'subset': list(colLabel_to_subset.values())})
d___.to_csv('../data/colLabel_to_subset_n_' + str(n) + '_n_iter_' + str(n_iter) + '.csv')


# In[75]:


print('Open and parse files')
# open and parse


lfc_string = '../data/' + normalization_type + '.lfc.n.'+ str(n) + '.n_iter.' + str(n_iter) + '.csv'
padj_string = '../data/' + normalization_type + '.padj.n.'+ str(n) + '.n_iter.' + str(n_iter) + '.csv'
bm_string = '../data/' + normalization_type + '.base_mean.n.'+ str(n) + '.n_iter.' + str(n_iter) + '.csv'
SE_string = "../data/" + normalization_type + '.' + "lfc.SE.n." + str(n) + ".n_iter." + str(n_iter) + ".csv"
disp_string = "../data/" + normalization_type + '.' + "lfc.dispersion.n." + str(n) + ".n_iter." + str(n_iter) + ".csv"


lfc_df = pd.read_csv(lfc_string, index_col = 0)
padj_df = pd.read_csv(padj_string, index_col = 0)
base_mean_df = pd.read_csv(bm_string, index_col = 0)
lfc_SE_df = pd.read_csv(SE_string, index_col = 0)
lfc_dispersion_df = pd.read_csv(disp_string, index_col = 0)


lfc_df.columns, padj_df.columns, base_mean_df.columns = subset_labels, subset_labels,subset_labels
lfc_SE_df.columns, lfc_dispersion_df.columns = subset_labels, subset_labels


cols = lfc_df.columns.tolist()


p = open("../data/actual_proportions.n."+ str(n) + '.n_iter.' + str(n_iter) + '.txt').read().splitlines()
p = sorted(set([float(i) for i in p]))

n_genes = lfc_df.shape[0]
genes = lfc_df.index

lfc_df.to_csv('../data/' + normalization_type + 'parsed_lfc_df_n' + str(n) + '_n_iter_' + str(n_iter) + '.csv')
padj_df.to_csv('../data/' + normalization_type + 'parsed_padj_df_n' + str(n) + '_n_iter_' + str(n_iter) + '.csv')
base_mean_df.to_csv('../data/' + normalization_type + 'parsed_base_mean_df_n' + str(n) + '_n_iter_' + str(n_iter) + '.csv')

lfc_SE_df.to_csv('../data/' + normalization_type + 'parsed_lfc_SE_df_n' + str(n) + '_n_iter_' + str(n_iter) + '.csv')
lfc_dispersion_df.to_csv('../data/' + normalization_type + 'parsed_lfc_dispersion_df_n' + str(n) + '_n_iter_' + str(n_iter) + '.csv')


# Generate cutoff values and apprporiate gene lists

# In[4]:


print('Threshold genes')

if (lfc_min_threshold == 0) or (base_mean_min_threshold == 0): 
    concordance_genes = genes
    correlation_genes = genes
elif min_percentage == 0 or max_cutoff == 0:
    
    # filter by thresholds
    median_df = pd.DataFrame(pd.np.empty((n_genes, 3)) * pd.np.nan) 
    median_df.columns = ['lfc', 'padj', 'base_mean']
    median_df.index = genes
    median_df['lfc'] = lfc_df.median(axis = 1).tolist()
    median_df['padj'] = padj_df.median(axis = 1).tolist()
    median_df['base_mean'] = base_mean_df.median(axis = 1).tolist()


    min_cutoffs_lfc = median_df.quantile(lfc_min_threshold)
    min_cutoffs_base_mean = median_df.quantile(base_mean_min_threshold)
    concordance_genes = median_df[(median_df['base_mean']>min_cutoffs_base_mean['base_mean'])].index.tolist()
    correlation_genes = median_df[(median_df['lfc']> min_cutoffs_lfc['lfc']) & (median_df['base_mean']>min_cutoffs_base_mean['base_mean'])].index.tolist()
else:
    # filter by thresholds
    median_df = pd.DataFrame(pd.np.empty((n_genes, 3)) * pd.np.nan) 
    median_df.columns = ['lfc', 'padj', 'base_mean']
    median_df.index = genes
    median_df['lfc'] = lfc_df.median(axis = 1).tolist()
    median_df['padj'] = padj_df.median(axis = 1).tolist()
    median_df['base_mean'] = base_mean_df.median(axis = 1).tolist()


    min_cutoffs_lfc = median_df.quantile(lfc_min_threshold)
    min_cutoffs_base_mean = median_df.quantile(base_mean_min_threshold)
    concordance_genes = median_df[(median_df['base_mean']>min_cutoffs_base_mean['base_mean'])].index.tolist()
    genes_ = median_df[(median_df['lfc']> min_cutoffs_lfc['lfc']) & (median_df['base_mean']>min_cutoffs_base_mean['base_mean'])].index.tolist()

    # if a gene has padj < max_cutoff for >= min_percentage*total_iterations, include it in the correlation analysis
    significant_hits = padj_df[padj_df < max_cutoff].count(axis = 1)
    correlation_genes = significant_hits[significant_hits > int(padj_df.shape[1]*min_percentage)].index.tolist()
    correlation_genes = sorted(set(correlation_genes).intersection(genes_))

with open('../data/' + normalization_type + 'correlation_genes_n_' + str(n) + '_n_iter_' + str(n_iter) + '.txt','w') as f:
    for i in correlation_genes:
        f.write(i + '\n')
        
with open('../data/' + normalization_type + 'concordance_genes_n_' + str(n) + '_n_iter_' + str(n_iter) + '.txt','w') as f:
    for i in concordance_genes:
        f.write(i + '\n')


# Initialize dfs to store data in

# In[78]:


print('Initialize dfs')
# initialize dfs
nrows = int(binom(n*n_iter,2))
ncols = int(len(p))

correlation_df = pd.DataFrame(pd.np.empty((nrows, ncols)) * pd.np.nan) 
correlation_df.columns = [str(i)[0:4] for i in p]
correlation_sig_df, subsets_df = correlation_df.copy(), correlation_df.copy()

concordance_df = pd.DataFrame(pd.np.empty((len(concordance_genes), 2*ncols)) * pd.np.nan)
concordance_columns = []
for col in correlation_df.columns:
    concordance_columns.append(col + ' SD')
    concordance_columns.append(col + ' mean abs(LFC)')
concordance_df.columns = concordance_columns
concordance_df.index = concordance_genes


# In[ ]:


def get_concordance(p_, lfc_df = lfc_df, concordance_genes = concordance_genes):
    '''Get gene-wise sd and mean(lfc) for a given proportion p_.'''
    proportion_col = str(p_)[0:4]
    p_cols = [i for i in cols if i[0:4] == proportion_col] 
    
    # filter by thresholds
    lfc_concordance = lfc_df.loc[concordance_genes,p_cols]
    
    conc_sd = lfc_concordance.std(axis = 1).tolist()
    conc_mean = lfc_concordance.abs().mean(axis = 1).tolist()
    
    return [conc_sd, conc_mean]


# In[5]:


def get_correlation_sig(p_,lfc_df = lfc_df, alpha = alpha, correlation_genes = correlation_genes):
    '''Get Spearman correlation with significance for a given proportion p_.'''

    proportion_col = str(p_)[0:4]
    p_cols = [i for i in cols if i[0:4] == proportion_col] 

    # generate subset df
    lfc_corr = lfc_df.loc[correlation_genes,p_cols]

    ## CORRELATION
    lfc_corr = lfc_corr.fillna(0)

    subset_labels = list(combinations(p_cols,2))
    corr_vals = []
    p_vals = []
    # for sc in subset_labels:
    for sc in subset_labels:
        c = list(spearmanr(lfc_corr[sc[0]].values, lfc_corr[sc[1]].values))
        corr_vals.append(c[0])
        p_vals.append(c[1])

    # Benjamini Hochberg correction
    p_vals = list(multipletests(p_vals, alpha=alpha, method='fdr_bh', is_sorted=False, returnsorted=False)[0])

    return list(zip(corr_vals, p_vals))
    

def get_subsets(p_):
    """Get actual subset combinations that data was generated on"""
    
    proportion_col = str(p_)[0:4]
    p_cols = [i for i in cols if i[0:4] == proportion_col] 
    
    subset_labels = list(combinations(p_cols,2))
    S = [colLabel_to_subset[i[0]] + ' ; ' + colLabel_to_subset[i[1]] for i in subset_labels]
    return S 


# In[ ]:


print('Get Correlation similarity with significance')
pool = Pool(processes=cr) 
cor = pool.map(get_correlation_sig, p)
pool.close()

print('Get concordance')
pool = Pool(processes=cr)
con = pool.map(get_concordance, p)
pool.close()

print('Get subsets')
pool = Pool(processes=cr)
s = pool.map(get_subsets,p)
pool.close()


print('Write similarity metrics')
for i in trange(len(p)):
    correlation_df.iloc[:,i] = cor[i]
    subsets_df.iloc[:,i] = s[i]
    
    index = i*2
    concordance_df.iloc[:,index] = con[i][0]
    concordance_df.iloc[:,index+1] = con[i][1]

binarize_dict = {True: 1, False: 0}
for i in trange(len(cor)):
    corr_vals = [k[0] for k in cor[i]]
    p_vals = [binarize_dict[k[1]] for k in cor[i]]
    correlation_df.iloc[:,i] = corr_vals
    correlation_sig_df.iloc[:,i] = p_vals
    

correlation_df.to_csv('../data/' + normalization_type + 'correlation_df_n_' + str(n) + '_n_iter_' + str(n_iter) + '.csv')
correlation_sig_df.to_csv('../data/' + normalization_type + 'correlation_sig_df_n_' + str(n) + '_n_iter_' + str(n_iter) + '.csv')
subsets_df.to_csv('../data/' + normalization_type + 'subsets_df_n_' + str(n) + '_n_iter_' + str(n_iter) + '.csv')
concordance_df.to_csv('../data/' + normalization_type + 'concordance_df_n_' + str(n) + '_n_iter_' + str(n_iter) + '.csv')


# # Get average metrics

# Get avg RIN and freeze thaw per subsets combination, serial

# Get avg RIN and freeze thaw per subsets combination, parallel<--don't use, slower than serial

# In[250]:


def avg_val(annot, subsets_df, row, col):
    """Gets average RIN and freeze thaw for subset combination in subsets_df[row,col]"""
    try:
        comb = subsets_df.loc[row,col]
        subset_1 = comb.split(';')[0].split(' ')[:-1]
        subset_2 = comb.split(';')[1].split(' ')[1:]
    except:
        raise ValueError('Error from parsing, freeze thaw, row {}, col {}'.format(row,col)) 
    try:             
        ft_1 = annot['final_freeze_thaw'][annot['Comment'].isin(subset_1)].sum()
        ft_2 = annot['final_freeze_thaw'][annot['Comment'].isin(subset_2)].sum()
        avg_ft = (ft_1 + ft_2)/(len(subset_1) + len(subset_2))
        if len(subset_1) != len(subset_2):
            raise ValueError('Subset parsing went wrong during avg ft and RIN calculation')
    except:
        raise ValueError('Error from freeze thaw calculation, row {}, col {}'.format(row,col))
    
    try:
        RIN_1 = annot['RINe'][annot['Comment'].isin(subset_1)].sum()
        RIN_2 = annot['RINe'][annot['Comment'].isin(subset_2)].sum()
        avg_RIN = (RIN_1 + RIN_2)/(len(subset_1) + len(subset_2))
    except:
        raise ValueError('Error from RIN calculation, row {}, col {}'.format(row,col))
    return [avg_ft, avg_RIN]

print('Load annot df')                
annot = pd.read_excel('../manifest.xlsx',sheet_name = 5)

def get_avg_FT(row, subsets_df = subsets_df, metric = 'Freeze Thaw', annot=annot):
    """Parallelize attaining average value metrics, iterates on rows of subsets_df"""
    
    if metric == 'Freeze Thaw':
        index = 0
    elif metric == 'RIN':
        index = 1
    
    vals = [avg_val(annot,subsets_df,row,col)[index] for col in subsets_df.columns] 
    return vals

def get_avg_RIN(row, subsets_df = subsets_df, metric = 'RIN', annot=annot):
    """Parallelize attaining average value metrics, iterates on rows of subsets_df"""
    
    if metric == 'Freeze Thaw':
        index = 0
    elif metric == 'RIN':
        index = 1
    
    vals = [avg_val(annot,subsets_df,row,col)[index] for col in subsets_df.columns] 
    return vals

print('Iterate through freeze thaw')
pool = Pool(processes=cr) 
d_ft = pool.map(get_avg_FT, subsets_df.index)
pool.close()
print('Put freeze thaw in df format')
freeze_thaw_df = pd.DataFrame(np.array(d_ft), columns = subsets_df.columns, index = subsets_df.index)
print('Save freeze thaw csv')
freeze_thaw_df.to_csv('../data/' + normalization_type + 'serial_average_freeze_thaw_df_n_' + str(n) + '_n_iter_' + str(n_iter) + '.csv')


print('Iterate through RIN')
pool = Pool(processes=cr) 
d_RIN = pool.map(get_avg_RIN, subsets_df.index)
pool.close()
print('Put RIN in df format')
RIN_df = pd.DataFrame(np.array(d_RIN), columns = subsets_df.columns, index = subsets_df.index)
print('Save RIN csv')
RIN_df.to_csv('../data/' + normalization_type + 'serial_average_RIN_df_n_' + str(n) + '_n_iter_' + str(n_iter) + '.csv')


# In[263]:


inv_map = {v: k for k, v in colLabel_to_subset.items()}
for i in range(10):
    coord = (random.randint(0, subsets_df.shape[0]-1), random.randint(0, subsets_df.shape[1]-1))
    test_sets = subsets_df.iloc[coord].split(' ; ')
    set1, set2 = test_sets[0], test_sets[1]
    actual_r = spearmanr(lfc_df.loc[correlation_genes, inv_map[set1]], lfc_df.loc[correlation_genes, inv_map[set2]])[0]
    if correlation_df.iloc[coord] != actual_r:
        raise ValueError('Problem with correlation df indexing')
    
    set1, set2 = set1.split(' '), set2.split(' ')
    actual_ft = annot['final_freeze_thaw'][annot['Comment'].isin(set1)].sum() + annot['final_freeze_thaw'][annot['Comment'].isin(set2)].sum()
    actual_ft = actual_ft/(len(set1) + len(set2))
    if actual_ft != freeze_thaw_df.iloc[coord]:
        raise ValueError('Freeze thaw parsing error')

    actual_RIN = annot['RINe'][annot['Comment'].isin(set1)].sum() + annot['RINe'][annot['Comment'].isin(set2)].sum()
    actual_RIN = actual_RIN/(len(set1) + len(set2))
    if actual_RIN != RIN_df.iloc[coord]:
        raise ValueError('RIN parsing error')


# # Summarize all data

# In[22]:


print('Generate summary metrics')
# generate summary metrics 
n_combs = correlation_df.shape[0]
proportion_col = [[p_]*n_combs for p_ in p]
proportion_col = [item for sublist in proportion_col for item in sublist]

subsets_col = [subsets_df[col].tolist() for col in subsets_df.columns]
subsets_col = [item for sublist in subsets_col for item in sublist]

correlation_col = [correlation_df[col].tolist() for col in correlation_df.columns]
correlation_col = [item for sublist in correlation_col for item in sublist]

correlation_sig_col = [correlation_sig_df[col].tolist() for col in correlation_sig_df.columns]
correlation_sig_col = [item for sublist in correlation_sig_col for item in sublist]

ft_col = [freeze_thaw_df[col].tolist() for col in freeze_thaw_df.columns]
ft_col = [item for sublist in ft_col for item in sublist] 

RIN_col = [RIN_df[col].tolist() for col in RIN_df.columns]
RIN_col = [item for sublist in RIN_col for item in sublist] 


summary_df = pd.DataFrame(data = {'Subset Combination': subsets_col, 'Proportion': proportion_col, 
                                 'Correlation': correlation_col, 
                                  'Avg Freeze Thaw': ft_col, 'Avg Rin': RIN_col, 
                                  'Correlation Significante': correlation_sig_col})
summary_df.to_csv('../data/' + normalization_type + 'summary_df_n_' + str(n) + '_n_iter_' + str(n_iter) + '.csv')


# In[ ]:


# for glms
sdf = summary_df.copy()
sdf['Proportion'] = sdf['Proportion'].map(dict(zip(sorted(set(summary_df['Proportion'])), [4,6,8,10,12,14])))
sdf.to_csv('../data/' + normalization_type + 'sdf_n_' + str(n) + '_n_iter_' + str(n_iter) + '.csv')

