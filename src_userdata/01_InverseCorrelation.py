#!/usr/bin/env python
# coding: utf-8

# # Calculating inverse correlations.

# In[38]:


import pandas as pd
import numpy as np
import os
from scipy.spatial.distance import cdist
import itertools
import argparse


# In[]:


parser = argparse.ArgumentParser()

## Perturbation type.
parser.add_argument(
    '--pert_type', 
    type=str, 
    default='trt_oe', 
    help='Perturbation type of protein signatures, e.g., trt_sh.cgs or trt_oe'
) # knockdown signatures

## Disease set.
parser.add_argument(
    '--cell_type', 
    type=str, 
    default='MCF7', 
    help='Disease id set, e.g., all or MCF7'
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
    '--out_file', 
    type=str, 
    default='../data_userdata/output/trt_oe/C0027662_Thyroid_MCF7.txt', 
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


# In[52]:


# =======================================================
# Setting.
# =======================================================

# ------- Basic setting ---------- #

## Perturbation type.
pert_type = args.pert_type
# pert_type = 'trt_sh.cgs'
# pert_type = 'trt_oe'

## Cell type.
cell_type = args.cell_type
# cell_type = "all"
# cell_type = "MCF7"

## Input TWAS file path.
input_twas = args.input_twas
# input_twas = "../data_userdata/example/twas_C0027662_Thyroid.txt"

## Output file path.
out_file = args.out_file
# out_file = "../data_userdata/output/C0027662_Thyroid_MCF7.txt"


# ------- Column setting ---------- #
col_zscore = args.col_zscore
# col_zscore = "zscore"
col_disease_id = args.col_disease_id
# col_disease_id = "disease_id"
col_tissue = args.col_tissue
# col_tissue = "tissue"
col_entrezgene = args.col_entrezgene
# col_entrezgene = "entrezgene"

# ## Association type.
# asso_type = "GDA_AeBmCmGvPm_VDA_CmGv"

# ## Disease selection method.
# dis_select = "Union" # "Product"：共通の疾患のみを使用, "Union"：共通でない疾患は0で置き換える

# ## calculation type.
# calc_type = 'H2geneEachCellDimred' # 一つのgeneに複数のeqtlが紐づく場合の処理

# ## Cell selcection method.
# # cell_select = "All" # "All"：平均化データのみ使用, "Each"：各細胞種を使用
# cell_select = "each" # "All"：平均化データのみ使用, "Each"：各細胞種を使用


# ### Genetically perturbed data.

# In[19]:


# =======================================================
# Target gene perturbation signatures.
# =======================================================

## Load imputed data.
# i_f = f'../data_userdata/LINCS/{pert_type}_original.txt'
i_f = f'../data_userdata/LINCS/{pert_type}.txt'
df = pd.read_csv(i_f, sep = '\t')

## Select cell line.
if cell_type != 'all':
    df = df[df['cell_mfc_name'] == cell_type]

cell_lines = list(set(df.cell_mfc_name)) # Cell lines.
pert_genes = sorted(list(set(df.cmap_name))) # Perturbed genes sorted.


df.head()


# ### Disease TWAS data.

# In[53]:


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

## Disease id.
dis_list = sorted(set(tw_df['disease_id']))
print("Number of diseases: {}".format(len(dis_list)))

## Tissue list.
tissue_list = sorted(set(tw_df['tissue']))
print("Number of tissues: {}".format(len(tissue_list)))

tw_df


# ## Landmark genes.

# In[33]:


# =======================================================
# Land mark genes
# =======================================================

# Landmark genes.
land_genes = sorted(set(df.columns))
land_genes = [s for s in land_genes if 'hsa' in s] # hsaを含むものだけ抽出

# Landmark gene dataframe (疾患プロファイル作成に使用).
land_df = pd.DataFrame()
for d in dis_list:
    tmp_df = pd.DataFrame({'entrezgene':land_genes})
    tmp_df['disease_id'] = d # 疾患id
    land_df = pd.concat( [land_df, tmp_df], axis=0)

print("cell line: ", len(cell_lines))
print("perturbed genes", len(pert_genes))

land_df


# ### Calculating inverse correlations.

# In[25]:


# =======================================================
# Combinations of cell types and tissues.
# =======================================================

comb_list = []

for c in cell_lines:
    for t in tissue_list:
        comb_list.append( (c, t) )
            
comb_list[0]


# In[35]:


# =======================================================
# Inverse correlations
# =======================================================
    
def corDimred(tmp):
    
    c = tmp[0] # cell line
    t = tmp[1] # tissue

    # ----- Disease profiles ----- #
    tmp_gw_gene_df = tw_df[tw_df['tissue'] == t].drop('tissue', axis = 1) # tissueを選択
    tmp_gw_gene_df = pd.merge( tmp_gw_gene_df, land_df, on = ['disease_id', 'entrezgene'], how = 'right') # land mark geneを紐づける
    tmp_gw_gene_df = tmp_gw_gene_df.pivot_table(index='entrezgene', columns= 'disease_id', values ='zscore').fillna(0) # 横持ち
    print( "Number of genes used for calculating correlations：{}".format( len(tmp_gw_gene_df)) )


    # ----- Perturbation profiles ----- #
    if c != 'All':
        tmp_p_df = df[df['cell_mfc_name'] == c].drop('cell_mfc_name', axis = 1).set_index('cmap_name').T # cellを選択
    else:
        tmp_p_df = df.groupby(by = 'cmap_name').mean(numeric_only=True).T # 平均化プロファイルを作成
    tmp_p_df = tmp_p_df.loc[ tmp_gw_gene_df.index.values.tolist() ] # disease profileと合わせる


    # ----- ndarray ----- #
    tmp_p_mt = tmp_p_df.T.values
    tmp_eq_mt = tmp_gw_gene_df.T.values

    # ----- Pearson's correlation coefficient ----- #
    # cdistを使用してPearson's rを計算
    # cdistは (1から相関係数を引いた値) を返す
    # つまり、相関係数は (1 - cdistの結果) となる
    tmp_cor_df = (1 - cdist(tmp_p_mt, tmp_eq_mt, metric='correlation'))
    tmp_cor_df = pd.DataFrame(tmp_cor_df, index=tmp_p_df.columns, columns=tmp_gw_gene_df.columns).fillna(0)
    tmp_cor_df = tmp_cor_df.stack().reset_index() # 縦持ち
    tmp_cor_df.columns = ['gene', 'disease_id', 'correlation'] # 列名
    tmp_cor_df['tissue'] = t # tissue
    tmp_cor_df['cell'] = c # cell line
    tmp_cor_df = tmp_cor_df.drop_duplicates() # 重複を除く

    return list(tmp_cor_df['gene']), list(tmp_cor_df['disease_id']), list(tmp_cor_df['correlation']), list(tmp_cor_df['tissue']), list(tmp_cor_df['cell'])


# In[36]:


# ===== 相関係数 ===== #
res_list = list(map(corDimred, comb_list))


# In[42]:


# =======================================================
# 結果のリストを整理する.
# =======================================================

# 結果のリストを整理する.
length = len(res_list)# 結果の長さ

res_gene_list = [ res_list[i][0] for i in range(length)] # gene
res_dis_list = [ res_list[i][1] for i in range(length)] # disease
res_cor_list = [ res_list[i][2] for i in range(length)] # correlation
res_tis_list = [ res_list[i][3] for i in range(length)] # tissue
res_cell_list = [ res_list[i][4] for i in range(length)] # cell

# ----- 平坦化 ----- #
res_gene_list = list(itertools.chain.from_iterable(res_gene_list))
res_dis_list = list(itertools.chain.from_iterable(res_dis_list))
res_cor_list = list(itertools.chain.from_iterable(res_cor_list))
res_tis_list = list(itertools.chain.from_iterable(res_tis_list))
res_cell_list = list(itertools.chain.from_iterable(res_cell_list))


# In[45]:


# =======================================================
# DataFrame
# =======================================================

# ===== DataFrame ===== #
cor_df = pd.DataFrame( {'gene':res_gene_list, 'disease_id':res_dis_list, 'correlation':res_cor_list, 'tissue':res_tis_list, 'cell': res_cell_list})

## Prediction score: -1 * correlation
cor_df['score'] = [ -1 * s for s in cor_df['correlation'] ]

## Prediction rank.
cor_df = cor_df.copy()
cor_df['prediction_rank'] = cor_df.groupby(by=['disease_id', 'tissue', 'cell'])['score'].rank(ascending = False)

## Sort data by prediction rank.
cor_df = cor_df.sort_values(by=['disease_id', 'tissue', 'cell', 'prediction_rank'])

cor_df.head()


# In[49]:


# =======================================================
# Save
# =======================================================

## Make a directory.
o_dir = out_file.replace(out_file.split('/')[-1], "")
os.makedirs(o_dir, exist_ok=True)

## Save the file.
o_f = out_file
cor_df.to_csv(o_f, sep = '\t', index = None)

print("Saved the results as a {} file.".format(out_file))

