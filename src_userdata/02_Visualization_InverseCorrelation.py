#!/usr/bin/env python
# coding: utf-8

# # 逆相関の結果を可視化

# In[1]:


import pandas as pd
import os, math
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import statistics


# In[]:


parser = argparse.ArgumentParser()

## Perturbation type.
parser.add_argument(
    '--pert_type', 
    type=str, 
    default='trt_sh.cgs', 
    help='Perturbation type of protein signatures, e.g., trt_sh.cgs or trt_oe'
) # knockdown signatures

## Protein cell type.
parser.add_argument(
    '--cell_type', 
    type=str, 
    default='MCF7', 
    help='Cell type of target gene perturbation signatures, e.g., all or MCF7'
)

## Disease tissue type.
parser.add_argument(
    '--tissue_type', 
    type=str, 
    default='Thyroid', 
    help='Disease tissue type, e.g., Thyroid'
)

## Disease identifier.
parser.add_argument(
    '--foc_d', 
    type=str, 
    default="C0027662", 
    help='Disease identifier, e.g., C0027662'
)

## therapeutic target gene name.
parser.add_argument(
    '--foc_p', 
    type=str, 
    default="FHL2,BECN1", 
    help='Therapeutic target candidates, e.g., FHL2,BECN1'
)

## Input file.
parser.add_argument(
    '--input_twas', 
    type=str, 
    default='../data_userdata/example/twas_C0027662_Thyroid.txt', 
    help='Input twas file path'
)

## Output file.
parser.add_argument(
    '--out_dir', 
    type=str, 
    default='../data_userdata/output/trt_sh.cgs/bar_plot', 
    help='Output file path'
)

## Column name of zscore.
parser.add_argument(
    '--col_zscore', 
    type=str, 
    default='zscore', 
    help='Column name of zscore'
)

## Column name of disease id.
parser.add_argument(
    '--col_disease_id', 
    type=str, 
    default='disease_id', 
    help='Column name of disease_id'
)

## Column name of tissue.
parser.add_argument(
    '--col_tissue', 
    type=str, 
    default='tissue', 
    help='Column name of tissue'
)

## Column name of entrezgene.
parser.add_argument(
    '--col_entrezgene', 
    type=str, 
    default='entrezgene', 
    help='Column name of entrezgene'
)

args = parser.parse_args()


# In[27]:


# ------- Basic setting ---------- #

## Perturbation type.
pert_type = args.pert_type
# pert_type = 'trt_sh.cgs'
# pert_type = 'trt_oe'

## Cell type.
cell_type = args.cell_type
# cell_type = "all"
# cell_type = "MCF7"

## Tissue type.
tissue_type = args.tissue_type
# tissue_type = "Thyroid"

## Input TWAS file path.
input_twas = args.input_twas
# input_twas = "../data_userdata/example/twas_C0027662_Thyroid.txt"

## Output directory path.
out_dir = args. out_dir
# out_dir = f"../data_userdata/output/{pert_type}/bar_plot"


# ------- Column setting ---------- #
col_zscore = args.col_zscore
# col_zscore = "zscore"
col_disease_id = args.col_disease_id
# col_disease_id = "disease_id"
col_tissue = args.col_tissue
# col_tissue = "tissue"
col_entrezgene = args.col_entrezgene
# col_entrezgene = "entrezgene"



# # Association type.
# asso_type = "GDA_AeBmCmGvPm_VDA_CmGv"

# # Clculation type.
# calc_type = "JI" # "JI" or "GO"

# # Process type.
# sim_pro = '_rute2' # 類似度の処理（2乗根など）なしなら""

# # Disease selection method.
# dis_select = "Union" # "Product"：共通の疾患のみを使用, "Union"：共通でない疾患は0で置き換える

# # Selected association type.
# sel_asso_type = "VDA_Cm"

# # Train label type.
# train_lab_type = "original"



# # ---- Unsupervised learning ----- #
# twas_type = 'InverseTwasEachCellDimred' # twas


# # ----- Bayesian optimization ----- #
# merge_type = "bayesian_1on1"


# ----- Focus disease ------ #
foc_d = args.foc_d
# if pert_type == 'trt_sh.cgs':
#     foc_d = "C0027662"
# elif pert_type == 'trt_oe':
#     foc_d = "C1800706"


# ------ Focus protein ------ #
foc_p = args.foc_p
# foc_p = "FHL2,BECN1"
foc_p = foc_p.split(',')
foc_p


# In[3]:


set2_list = sns.color_palette("Set2", 8)
color_list = [set2_list[i] for i in [3] ]


# In[4]:


sns.color_palette("Set2")


# In[ ]:





# In[5]:


# ===== LINCS gene data ===== #
i_f = f"../data_userdata/LINCS/{pert_type}_all_converted_gene.txt"
linc_df = pd.read_csv( i_f, sep = '\t' )

linc_gene_list = sorted(set(linc_df['gene'])) # LINCSの遺伝子のうち、解釈可能な遺伝子（CMAP***のように何のタンパクか分からないものは除く）

linc_df.head()


# In[ ]:





# In[20]:


# =======================================================
# TWAS data.
# =======================================================

## Load data.
# i_f = f"../data_userdata/example/twas_C0027662_Thyroid.txt"
i_f = input_twas
tw_df = pd.read_csv(i_f, sep = '\t')

## Rename columns.
tw_df = tw_df.rename({
    col_zscore: "zscore",
    col_disease_id: "disease_id",
    col_tissue: "tissue",
    col_entrezgene: "entrezgene"
})

## -log(p-value)
tw_df['Mlogp'] = [ -1* math.log(s) if s != 0 else 30 for s in tw_df['pvalue']]

## Select the prespecified disease.
tw_df = tw_df[tw_df['disease_id'] == foc_d]

## Select the prespecified tissue.
tw_df = tw_df[tw_df['tissue'] == tissue_type]


## -------- List --------- ##
## TWAS genes
tw_gene = sorted(set(tw_df['entrezgene']))

## Disease id.
dis_list = sorted(set(tw_df['disease_id']))
print("Number of diseases: {}".format(len(dis_list)))

## Tissue list.
tissue_list = sorted(set(tw_df['tissue']))
print("Number of tissues: {}".format(len(tissue_list)))

tw_df


# In[ ]:





# In[21]:


# =======================================================
# Target gene perturbation signatures (TGPs)
# =======================================================

# Load imputed data.
# i_f = f'../data_userdata/LINCS/{pert_type}_original.txt'
i_f = f'../data_userdata/LINCS/{pert_type}.txt'
df = pd.read_csv(i_f, sep = '\t')


df = df[df['cell_mfc_name'] == cell_type] # 細胞を選択
df = df.set_index( 'cmap_name' ) # indexを指定
df = df.drop('cell_mfc_name', axis=1) # 列を削除



df = df.stack().reset_index() # 縦持ち
df.columns = ['target', 'entrezgene', 'target_value'] # 列名変更

df.head()


# In[22]:


def min_max(l):
    l_min = min(l)
    l_max = max(l)
    if (l_max - l_min) != 0:
        return [(i - l_min) / (l_max - l_min) for i in l]
    else:
        return 0

def standardization(l):
    l_mean = statistics.mean(l)
    l_stdev = statistics.stdev(l)
    return [(i - l_mean) / l_stdev for i in l]


# In[28]:


# ===================================================================
# Barplot of TRESOR and target gene perturbation signatures
# ===================================================================

for p in foc_p:

    ## Select the target.
    tmp_df = df[df['target'] == p]

    ## Merge TRESOR and TGP signatures.
    mg_df = pd.merge( tw_df, tmp_df, on = 'entrezgene', how = 'inner').fillna( {'target':p, 'target_value':0} )

    # ---- Normalize zscores ----- #
    ## Whether zscores are positive or not.
    mg_df['positive'] = [1 if s > 0 else 0 for s in mg_df['zscore']]
    mg_df['zscore'] = mg_df.groupby( by = ['positive'])['zscore'].transform(min_max) # Normalize from 0 to 1
    mg_df['zscore'] = [ s-1 if a == 0 else s   for s,a   in zip(mg_df['zscore'], mg_df['positive'])] # Restore to original direction of values

    # ---- Normalize target_values ----- #
    ## Whether zscores are positive or not.
    mg_df['target_positive'] = [1 if s > 0 else 0 for s in mg_df['target_value']]
    mg_df['target_value'] = mg_df.groupby( by = ['target_positive'])['target_value'].transform(min_max) # Normalize from 0 to 1
    mg_df['target_value'] = [ s-1 if a == 0 else s   for s,a   in zip(mg_df['target_value'], mg_df['target_positive'])] # Restore to original direction of values


    ## Sort genes by zscores.
    mg_df = mg_df.sort_values('zscore', ascending=False)


    # ----- Plot ----- #
    plt.figure(figsize=[20,5])
    plt.rcParams['axes.xmargin'] = 0.01 # 横の余白を小さくする
    sns.barplot(data= mg_df, x='gene_name', y='zscore', color=set2_list[3], errcolor='k', linewidth = 0.7, edgecolor = 'k')
    sns.barplot(data= mg_df, x='gene_name', y='target_value', color='gray', linewidth = 0.7, edgecolor = 'k', alpha = 0.7)
    plt.xlabel('Genes')
    plt.ylabel('value')
    plt.xticks(rotation = 90)
    plt.tick_params(labelsize=14)
    plt.axhline(y=0, color = 'k', linestyle = '--')
    
    # ----- Save the data ----- #
    o_dir =f"{out_dir}"
    os.makedirs( o_dir, exist_ok=True )
    
    o_f = f"{o_dir}/{foc_d}_{p}_barplot.png"
    plt.savefig(o_f, bbox_inches = 'tight')

