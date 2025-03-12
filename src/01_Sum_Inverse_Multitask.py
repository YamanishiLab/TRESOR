#!/usr/bin/env python
# coding: utf-8

# # Merge prediction socres of multitask learning method and inverse signature method.

# In[1]:


import pandas as pd
import numpy as np
import os, glob, copy
import argparse


# In[2]:

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
    '--disease_set', 
    type=str, 
    default='all', 
    help='Disease id set, e.g., all or C0002736,C0002395'
) # knockdown signatures

## Output file.
parser.add_argument(
    '--out_file', 
    type=str, 
    default='../data/Bayesian/performance_evaluation/score_sum/trt_oe_select.txt', 
    help='Output file path'
) # knockdown signatures

args = parser.parse_args()


# In[24]:


# ========================================================
# Setting.
# ========================================================

## Perturbed type.
pert_type = args.pert_type
# pert_type = 'trt_sh.cgs'
# pert_type = 'trt_oe'

## Disease sets.
disease_set = args.disease_set
# disease_set = "all" ## 全疾患に対して予測
# disease_set = "C0027662,C1292769" ## 指定した疾患に対して予測
if disease_set != "all":
    disease_set = disease_set.split(',')


## Output file path.
out_file = args.out_file
# out_file = f"../data/Bayesian/performance_evaluation/score_sum/trt_oe_select.txt"


## Association type.
asso_type = "GDA_AeBmCmGvPm_VDA_CmGv"


# ---- Multitask learinng method ----- #

## Clculation type.
calc_type = "JI" # "JI" or "GO"

## Process type.
sim_pro = '_rute2' # 類似度の処理（2乗根など）なしなら""

## Disease selection method.
dis_select = "Union" # "Product"：共通の疾患のみを使用, "Union"：共通でない疾患は0で置き換える

## Train label type.
train_lab_type = "original"


# ---- Unsupervised learning ----- #
twas_type = 'InverseTwasEachCellDimred' # twas


# ### Label data.

# In[8]:


# ========================================================
# Gold standard data
# ========================================================

## Read data.
i_f = f"../data/LabelData/{pert_type}/{dis_select}_{asso_type}/label.txt"
gold_df = pd.read_csv(i_f, sep = '\t').rename(columns={'disease_id':'disease'})
gold_df = gold_df[['disease', 'gene']] # 列を選択
gold_df['label'] = 1

## protein list
pro_list = sorted(set(gold_df['gene']))

gold_df.head()


# ### Multitask learning method.

# In[9]:


def min_max(l):
    l_min = min(l)
    l_max = max(l)
    return [(i - l_min) / (l_max - l_min) for i in l]


# In[11]:


# ========================================================
# Multitask learning method
# ========================================================

## Read data.
i_f = f"../data/Bayesian/performance_evaluation/multitask/{pert_type}_merge.txt"
mt_df = pd.read_csv(i_f, sep = '\t')

# ---- min-max scaling ---- #
mt_df['mt_score'] = mt_df.groupby( by = ['disease', 'sim', 'l1', 'l2'])['score'].transform(min_max) # 0-1にノーマライズ
mt_df = mt_df.drop( 'score', axis=1) # scoreを削除

## Select diseases.
if disease_set != "all":
    mt_df = mt_df[mt_df['disease'].isin(disease_set)]

mt_df.head()


# #### Inverse signature method

# In[15]:


# ========================================================
# TWAS profilesがある疾患
# ========================================================

# ----- 疾患特異的かどうかのデータ ----- #
i_f = f"../data/Bayesian/performance_evaluation/inverse/disease_specific.txt"
tmp_df= pd.read_csv(i_f, sep = '\t')
tmp_df = tmp_df[tmp_df['non_specific'] == 0] # specificを選択
twas_dis_list = sorted(set(tmp_df['disease_id'])) # disease


# ----- 疾患と関連する組織・細胞の情報 ----- #
i_f = f"../data/Bayesian/performance_evaluation/inverse/Tissue_Manual_Selected.txt"
tw_tissue_df = pd.read_csv(i_f, sep = '\t', encoding='shift-jis')
tw_tissue_df = tw_tissue_df[tw_tissue_df['disease'].isin(twas_dis_list) ] # specificな疾患を選択

if pert_type == 'trt_sh.cgs':
    tw_tissue_df = tw_tissue_df.drop('oe_cell', axis = 1).rename(columns = {'kd_cell':'cell'}) # overexpressedのcellの列を削除
elif pert_type == 'trt_oe':
    tw_tissue_df = tw_tissue_df.drop('kd_cell', axis = 1).rename(columns = {'oe_cell':'cell'}) # knockdownのcellの列を削除

## 列を選択
tw_tissue_df = tw_tissue_df[['disease', 'tissue', 'cell']]

## Select diseases.
if disease_set != 'all':
    tw_tissue_df = tw_tissue_df[tw_tissue_df['disease'].isin(disease_set)]

## Disease list.
twas_dis_list = sorted(set(tw_tissue_df['disease']))

tw_tissue_df.head()


# In[17]:


# ========================================================
# Prediction scores of inverse signatue method with TRESOR
# ========================================================

## Read data.
i_f = f"../data/Bayesian/performance_evaluation/inverse/{pert_type}_prediction_score_label.txt"
tw_df = pd.read_csv(i_f, sep = '\t')

## Select tissues and cell lines.
tw_df = pd.merge( tw_df, tw_tissue_df, on = ['disease', 'tissue', 'cell'], how = 'inner')
tw_df = tw_df[tw_df['gene'].isin(pro_list)] # proteinを選択

# ---- min-max scaling ---- #
tw_df['tw_score'] = tw_df.groupby( by = ['disease', 'tissue', 'cell'])['score'].transform(min_max) # 0-1にノーマライズ
tw_df = tw_df.drop( ['correlation', 'score'], axis=1) # 列を削除

tw_df.head()


# ### Merge prediction scores of multitask and inverse.

# In[20]:


# ========================================================
# inverse signature methodとinverse signature methodの結果をマージ 
# ========================================================

mg_df = pd.merge( mt_df, tw_df, on = ['gene', 'disease', 'label'], how = 'left' ) # twasがない疾患は0で埋める
mg_df = mg_df.fillna( { 'tissue':np.nan, 'cell': np.nan, 'tw_score':0 } )
mg_df['score'] = [ a+b for a,b in zip( mg_df['mt_score'], mg_df['tw_score'] ) ] # mtとinvの和を計算

mg_df.head()


# In[25]:


# ========================================================
# Save
# ========================================================

## Make directory.
o_dir = out_file.replace(out_file.split('/')[-1],"")
# # o_dir = f"../data/Bayesian/performance_evaluation/score_sum/"
os.makedirs(o_dir, exist_ok=True)

## Save the data.
o_f = f"{out_file}"
# o_f = f"{o_dir}/{pert_type}.txt"
mg_df.to_csv(o_f, sep = '\t', index = None)

