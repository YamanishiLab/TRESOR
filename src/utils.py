#!/usr/bin/env python
# coding: utf-8

# # Merge parameters across all diseases.

# In[1]:


import pandas as pd
# import numpy as np
import os, glob
# import re
# import itertools
# import argparse
# from datatable import dt, f, fread, join


# In[30]:


# =================================================================
# Merge parameters across diseases.
# =================================================================

class MergeParams():

    def __init__(
        self,
        pert_type,
        disease_set,
        input_dir,
        out_dir
    ):

        super(MergeParams, self).__init__()
        self.pert_type = pert_type
        self.disease_set = disease_set
        self.input_dir = input_dir
        self.out_dir = out_dir


        # ------ Multitask learning method setting --------- #
        
        ## Disease selection method
        self.dis_select = "Union"

        ## Association type
        self.asso_type = "GDA_AeBmCmGvPm_VDA_CmGv"

        ## Disease similarity type of VDAs or GDAs
        self.sel_asso_type = "VDA_Cm"

    
    def GetDiseaseList(self):

        ## Load goldstandard data.
        i_f = f"../data/LabelData/{self.pert_type}/{self.dis_select}_{self.asso_type}/label.txt"
        gold_df = pd.read_csv(i_f, sep = '\t').rename(columns={'disease_id':'disease'})
        gold_df = gold_df[['disease', 'gene']] # 列を選択
        gold_df['label'] = 1
        
        pro_list = sorted(set(gold_df['gene'])) # protein list
        
        ## Disease list.
        if self.disease_set == "all":
            dis_list = sorted( set(gold_df['disease']) ) # disease
        else:
            dis_list = self.disease_set.split(",")

        return dis_list


    def Merge(self, dis_list):

        input_dir = self.input_dir
        out_dir = self.out_dir

        opt_df = pd.DataFrame()
        
        for d in dis_list:

            ## Load estimated parameters.
            i_f = f"{input_dir}/{self.pert_type}/{d}/opt.csv"
            # i_f = f"../data/Bayesian/performance_evaluation/parameters_disease/{pert_type}/{d}/opt.csv"
            tmp_df = pd.read_csv(i_f, sep = ',')
            tmp_df['disease'] = d
        
            opt_df = pd.concat( [opt_df, tmp_df], axis=0) # merge
        else:
            opt_df = opt_df.drop('Y', axis=1) # 列を削除
            opt_df.columns = ['tw_param', self.sel_asso_type + "_param", 'disease']
        
            ## Save.
            o_dir = out_dir
            os.makedirs(o_dir, exist_ok=True)
        
            o_f = f"{o_dir}/parameters.txt"
            # o_f = f"../data/Bayesian/performance_evaluation/parameters/{pert_type}.txt"
            opt_df.to_csv( o_f, sep = '\t', index=None )


# In[9]:


# # =================================================================
# # Setting.
# # =================================================================

# # ----------- Basic setting ------------- #
# ## Perturbed type
# # pert_type = 'trt_sh.cgs' # inhibitory target prediction
# pert_type = 'trt_oe' # activatory target prediction

# ## Disease set.
# # disease_set = "all" # All diseases
# disease_set = "C0002736,C0002395" ## Selected diseases
# if disease_set != "all":
#     disease_set = disease_set.split(",")

# ## Output file directory path.
# input_dir = "../data/Bayesian/performance_evaluation/parameters_disease"

# ## Output file path.
# out_dir = f"../data/Bayesian/new_prediction/output/{pert_type}"



# # ---- Multitask learinng method ----- #

# ## Association type
# asso_type = "GDA_AeBmCmGvPm_VDA_CmGv"

# ## Clculation type
# calc_type = "JI" # "JI" or "GO"

# ## Process type
# sim_pro = '_rute2' # 類似度の処理（2乗根など）なしなら""

# ## Disease selection method
# dis_select = "Union" # "Product"：共通の疾患のみを使用, "Union"：共通でない疾患は0で置き換える

# ## Disease similarity type of VDAs or GDAs
# sel_asso_type = "VDA_Cm"

# ## Train label type
# train_lab_type = "original"

# ## Inverse signature method type
# twas_type = 'InverseTwasEachCellDimred' # twas

# ## Bayesian optimization type
# merge_type = "bayesian_1on1"


# ## Label data.

# In[3]:


# # =================================================================
# # Gold standard data
# # =================================================================

# i_f = f"../data/LabelData/{pert_type}/{dis_select}_{asso_type}/label.txt"
# gold_df = pd.read_csv(i_f, sep = '\t').rename(columns={'disease_id':'disease'})
# gold_df = gold_df[['disease', 'gene']] # 列を選択
# gold_df['label'] = 1

# pro_list = sorted(set(gold_df['gene'])) # protein list

# ## Disease list.
# if disease_set == "all":
#     dis_list = sorted( set(gold_df['disease']) ) # disease
# else:
#     dis_list = disease_set

# gold_df.head()


# ### Merge params.

# In[10]:


# # =================================================================
# # ベイズ最適化のパラメータ
# # =================================================================

# opt_df = pd.DataFrame()

# for d in dis_list:

#     i_f = f"{input_dir}/{pert_type}/{d}/opt.csv"
#     # i_f = f"../data/Bayesian/performance_evaluation/parameters_disease/{pert_type}/{d}/opt.csv"
#     tmp_df = pd.read_csv(i_f, sep = ',')
#     tmp_df['disease'] = d

#     opt_df = pd.concat( [opt_df, tmp_df], axis=0) # merge
# else:
#     opt_df = opt_df.drop('Y', axis=1) # 列を削除
#     opt_df.columns = ['tw_param', sel_asso_type + "_param", 'disease']

#     ## Save.
#     o_dir = input_dir.replace(input_dir.split('/')[-1],"")
#     o_dir = f"{o_dir}/parameters/"
#     os.makedirs(o_dir, exist_ok=True)

#     o_f = f"{o_dir}/{out_file.split('/')[-1]}"
#     # o_f = f"../data/Bayesian/performance_evaluation/parameters/{pert_type}.txt"
#     opt_df.to_csv( o_f, sep = '\t', index=None )
    
# opt_df.head()

