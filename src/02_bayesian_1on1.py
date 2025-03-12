#!/usr/bin/env python
# coding: utf-8

# # Estimate optimal weights of inverse signature and multitask learning methods by Bayesian optimization.

# In[1]:


import pandas as pd
import numpy as np
import os, glob, json, math
# import GPy
import GPyOpt
from sklearn import metrics
# import matplotlib.pyplot as plt
# import datetime
from scipy.stats import rankdata
import sys
import argparse
sys.path.append('../BED_AUC/')
from scr.Enrichment import *


# In[]:

parser = argparse.ArgumentParser()

## Target gene perturbation signatures
parser.add_argument(
    '--pert_type', 
    type=str, 
    default='trt_oe', 
    help='Perturbation type of protein signatures, e.g., trt_sh.cgs or trt_oe'
)

## Disease set.
parser.add_argument(
    '--disease_set', 
    type=str, 
    default='all', 
    help='Disease id set, e.g., all or C0002736,C0002395'
)

## Input file path.
parser.add_argument(
    '--input_file', 
    type=str, 
    default='../data/Bayesian/performance_evaluation/score_sum/trt_oe_select.txt', 
    help='Input file path'
)

## Output file directory path.
parser.add_argument(
    '--out_dir', 
    type=str, 
    default='../data/Bayesian/performance_evaluation/parameters_disease', 
    help='Output file directory path'
)

args = parser.parse_args()


# In[53]:


# =================================================================
# Setting.
# =================================================================

# ----------- Basic setting ------------- #
## Perturbed type.
pert_type = args.pert_type
# pert_type = 'trt_sh.cgs'
# pert_type = 'trt_oe'


## Disease set.
disease_set = args.disease_set
# disease_set = "all" # All diseases
# disease_set = "C0002736,C0002395" ## Selected diseases
if disease_set != "all":
    disease_set = disease_set.split(",")


## Input file path.
input_file = args.input_file
# if pert_type == 'trt_sh.cgs' and disease_set != "all":
#     input_file = f"../data/Bayesian/performance_evaluation/score_sum/trt_sh.cgs_select.txt"
# elif pert_type == 'trt_sh.cgs' and disease_set == "all":
#     input_file = f"../data/Bayesian/performance_evaluation/score_sum/trt_sh.cgs.txt"
# elif pert_type == 'trt_oe' and disease_set != "all":    
#     input_file = f"../data/Bayesian/performance_evaluation/score_sum/trt_oe_select.txt"
# elif pert_type == 'trt_oe' and disease_set == "all":
#     input_file = f"../data/Bayesian/performance_evaluation/score_sum/trt_oe.txt"


## Output file directory path.
out_dir = args.out_dir
# out_dir = "../data/Bayesian/performance_evaluation/parameters_disease"


# Association type.
asso_type = "GDA_AeBmCmGvPm_VDA_CmGv"



# ---- Multitask learinng method ----- #

# Clculation type.
calc_type = "JI" # "JI" or "GO"

# Process type.
sim_pro = '_rute2' # 類似度の処理（2乗根など）なしなら""

# Disease selection method.
dis_select = "Union" # "Product"：共通の疾患のみを使用, "Union"：共通でない疾患は0で置き換える

# Train label type.
train_lab_type = "original"

# Selected association type.
sel_asso_type = "GDA_Cm"



# ------------ Unsupervised learning --------------- #
twas_type = 'InverseTwasEachCellDimred' # twas


# ------------ fold list ------------- #
fold_list = [1, 2, 3, 4, 5]


# ### Label data.

# In[54]:


# =================================================================
# Gold standard data
# =================================================================

i_f = f"../data/LabelData/{pert_type}/{dis_select}_{asso_type}/label.txt"
gold_df = pd.read_csv(i_f, sep = '\t').rename(columns={'disease_id':'disease'})
gold_df = gold_df[['disease', 'gene']] # 列を選択
gold_df['label'] = 1

pro_list = sorted(set(gold_df['gene'])) # protein list

## Disease list.
if disease_set == "all":
    dis_list = sorted( set(gold_df['disease']) ) # disease
else:
    dis_list = disease_set

gold_df.head()


# In[55]:


# =================================================================
# Merged prediction scores of inverse signature and multitask learning methods.
# =================================================================

## Read data.
# i_f = f"../data/Bayesian/performance_evaluation/score_sum/{pert_type}.txt"
i_f = input_file
score_df = pd.read_csv(i_f, sep = '\t')
score_df = score_df.drop("score", axis=1) # 列を削除
score_df = score_df[score_df['sim'] == sel_asso_type] # association typeを選択

asso_list = sorted( set(score_df['sim'])) # association type

score_df.head()


# ### Bayesian integrative method.

# In[56]:


# =================================================================
# Bayesian optimization.
# =================================================================

max_iter = 50
acquisition_type = 'LCB'
bounds = [{'name': 'w0', 'type': 'continuous', 'domain': (0,10)}, # 最適化する変数の名前（適当で良い）、連続値か離散値、値の探索範囲
                  {'name': 'w1', 'type': 'continuous', 'domain': (0,10)}]
random_state = 1
# acquisition_weight = 1000 # 獲得関数のパラメータ、大きくすると活用より探索重視、デフォルト2


# In[57]:


def PerformBayesianOptimization( cond ):
    obj_type = cond[0] # auc, aupr, bed auc
    d = cond[1] # disease

    def calculate_bedauc_score(w):
        w0, w1 = w[:,0], w[:,1]
        auc_list = []

        for i in fold_list:
            tmp_df = sel_score_df[ sel_score_df['fold'] == i ] # foldデータを選択

            # Rank list.
            cor_index = [i for i,s in enumerate(tmp_df['label']) if s == 1] # 正例のインデックス
            y_score = [w0*t + w1*m for t,m in zip(tmp_df['tw_score'], tmp_df['mt_score'] )] # 予測スコア
            tmp_rank = rankdata(-np.array(y_score), method='max') # 予測スコアの順位（降順）, 全体での順位に使用するためmax（minにすると過大評価しすぎる）
            y_rank = sorted( [tmp_rank[i] for i in cor_index] )  # 全体における正例の順位、ソート

            # Count list.
            y_count = rankdata(y_rank, method='min') # 正例の中での順位（昇順）, 正例の割合に使用するためmin.

            # Build BEDROS method.
            total = len(tmp_df) # Total pairs.
            enrich = Enrichment()
            enrich.readListData(y_rank, y_count, total)
            auc = enrich.calculateBEDROC(20.0)
            auc_list.append( auc )

        else:
            auc_mean = np.mean(auc_list)

        return auc_mean
    
    
    
    # ---- 疾患を選択 ----- #
    sel_score_df = score_df[score_df['disease'] == d] # disease
    
    
    # ===== ベイズ最適化 ===== #

    # Bayesian optimization.
    myBopt = GPyOpt.methods.BayesianOptimization(f=calculate_bedauc_score, 
                                                 domain=bounds, 
                                                 acquisition_type=acquisition_type,
                                                maximize=True,
                                                random_state=random_state)
    myBopt.run_optimization(max_iter=max_iter)

    # 最適値の抽出
    opt_fx = myBopt.fx_opt * (-1)
    opt_x = myBopt.x_opt



    # ===== 最適化結果をcsvファイルに保存する ===== #
    columns_name_list = ['tw_score', sel_asso_type] # 列名を作成
    columns_name_list_str = ','.join(columns_name_list) # 要素間にカンマを入れて文字列へ
    columns_name_list_str = columns_name_list_str + ',Y' + '\n' # ターゲット名と改行コードを入れる

    o_dir = f"{out_dir}/{pert_type}/{d}/"
    os.makedirs(o_dir, exist_ok=True)
    o_f = f"{o_dir}/opt.csv"
    with open( o_f, 'w') as f:
        f.write(columns_name_list_str)
        buf_list = []
        for X in opt_x.tolist():
            buf_list.append(X)
        buf_list.append(opt_fx)
        opt_list = list(map(str, buf_list))
        opt_str = ','.join(opt_list)
        f.write(opt_str + '\n')


    # ===== 収束状況をプロット ===== #
    # 隣接する測定点の3次元空間内での距離と、評価関数のbest値（最大化の場合は-1を掛けて考える）
    o_f = f"{o_dir}/{obj_type}_convergenece.png"
    myBopt.plot_convergence(filename=o_f)


    # ===== 探索履歴をcsvファイルに保存する ===== #
    history_Y = pd.DataFrame(myBopt.Y*(-1), columns = [obj_type])
    history_X = pd.DataFrame(myBopt.X, columns=columns_name_list)
    hist_df = pd.concat([history_X, history_Y], axis = 1)

    o_f = f"{o_dir}/{obj_type}_history_data.csv"
    hist_df.to_csv(o_f, sep = '\t', index=None)



    # ===== 最終的な予測スコアを計算 ===== #
    result_df = sel_score_df.copy()

    # Calculate total scores.
    result_df['score'] = [myBopt.x_opt[0]*t + myBopt.x_opt[1]*m   for t,m in zip( result_df['tw_score'], result_df['mt_score'] )]

    # Save.
    o_f = f"{o_dir}/Bayesian{obj_type}_label.txt"
    result_df.to_csv(o_f, sep = '\t', index=None )

    print(myBopt.x_opt)
    print(myBopt.fx_opt)

    return myBopt.x_opt


# In[58]:


if __name__ == "__main__":
    
    # ---- 計算条件 ------ #
    comb_list = []
    for v in ['bedauc']:
        for d in dis_list:
            comb_list.append( (v, d) )
            
#     reulst_list = list(map(PerformBayesianOptimization, ['aupr', 'auc', 'bedauc']))
    reulst_list = list(map(PerformBayesianOptimization, comb_list))

