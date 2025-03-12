#!/usr/bin/env python
# coding: utf-8

# # Newly Prediction using Bayesian integrative method
# 
# This code is for predicting new therapeutic targets using Bayesian integrative method that is proposed in the following paper:
# "Therapeutic target prediction for orphan diseases, integrating genome-wide and transcriptome-wide association study data".
# 
# ・Authors:  **    
# 
# ・Affiliations: **  

# In[1]:


import pandas as pd
import numpy as np
import os, glob
import argparse

import sys
sys.path.append('./')
from utils import *


# In[2]:


def get_parser():
    parser = argparse.ArgumentParser(
        description='Script for predicting therapeutic targets using Bayesian integrative method'
    )

    ## Perturbation type.
    parser.add_argument(
        "-pert_type", "--pert_type", 
        default = "trt_oe",
        type = str, help = "set therapeutic target type: inhibitory target, 'trt_sh.cgs'; activatory target, 'trt_oe' "
    )

    ## Disease set.
    parser.add_argument(
        '--disease_set', 
        default='all', 
        type=str, 
        help='Disease id set, e.g., all or C0002736,C0002395'
    )

    ## Input file directory path.
    parser.add_argument(
        '--input_dir', 
        type=str, 
        default='../data/Bayesian/performance_evaluation/parameters_disease', 
        help='Output file directory path'
    )

    ## Output file path.
    parser.add_argument(
        "-output", "--out_dir", 
        default='../data/Bayesian/new_prediction/output/trt_oe',
        type = str, 
        help = "set the output directory path"
    )

    ## Selected association type.
    parser.add_argument(
        "--sel_asso_type", 
        default='VDA_Cm',
        type = str, 
        help = "Seleted associaton type in the multitask learning method"
    )

    ## Disease selection method
    parser.add_argument(
        "--dis_select", 
        default='Union',
        type = str, 
        help = "Seleted disease set used for performance evauation"
    )

    ## Association type
    parser.add_argument(
        "--asso_type", 
        default='GDA_AeBmCmGvPm_VDA_CmGv',
        type = str, 
        help = "GDA and VDA types."
    )

    try:
        args = parser.parse_args()
    except:
        args = parser.parse_args(args = [])
    return args


# In[4]:


# ===== Min-max scaling ===== #
def min_max(l):
    l_min = min(l)
    l_max = max(l)
    return [(i - l_min) / (l_max - l_min) for i in l]


# In[5]:


# ===== Bayesian integrative method ====== #

def BImethod( args, dis_list ):

    pert_type = args.pert_type
    out_dir = args.out_dir
    sel_asso_type = args.sel_asso_type



    # ----- Load prediction scores of inverse signature and multitask learning methods ----- #

    # Multitask learning method
    i_f = f"../data/Bayesian/new_prediction/pre_prediction_scores/{pert_type}_multi_task.txt"
    mt_df = pd.read_csv( i_f, sep = '\t' )

    # Inverse signature method
    i_f = f"../data/Bayesian/new_prediction/pre_prediction_scores/{pert_type}_inverse_signature.txt"
    tw_df = pd.read_csv( i_f, sep = '\t' )

    # Parameters of Bayesian integrative method
    i_f = f"{out_dir}/parameters.txt"
    # i_f = f"../data/Bayesian/performance_evaluation/parameters/{pert_type}.txt"
    opt_df = pd.read_csv( i_f, sep = '\t' )

    ## Disease id–name dict.
    i_f = f"../data/LabelData/{pert_type}/{args.dis_select}_{args.asso_type}/disease_id_name.txt"
    dis_dict = pd.read_csv(i_f, sep = '\t')
    dis_dict = dict(zip(dis_dict.disease, dis_dict.diseaseName))



    # ----- Newly prediction using Bayesian integrative method ----- #
    mg_df = pd.merge( mt_df, tw_df, on = ['gene', 'disease'], how = 'left' ) # merge data frames.
    mg_df = mg_df.fillna( { 'tissue':np.nan, 'cell': np.nan, 'tw_score':0 } )# impute parameters as 0 for diseases without TRESOR.
    mg_df = pd.merge( mg_df, opt_df, on = 'disease', how = 'left') # add parameter information.
    mg_df['score'] = [ s*a+ v*b 
                      for a,b, s,v in zip( mg_df['tw_score'], mg_df[sel_asso_type],  mg_df['tw_param'], mg_df[sel_asso_type + '_param'] ) ] # Calculate prediction scores.
    mg_df['normalized_score'] = mg_df.groupby( by = ['disease'])['score'].transform(min_max) # normalize prediction scores.



    # ----- Select and rename columns ----- #
    mg_df = mg_df[['gene', 'disease', 'tw_param', f'{sel_asso_type}_param', 'normalized_score']] # select columns.
    mg_df.columns = ['gene', 'disease', 'inv_param', 'mul_param', 'prediction_score'] # rename columns
    mg_df['diseaseName'] = [dis_dict[s] for s in mg_df['disease']] ## Add diseaseName.


    # ------ Select disease list ------ #
    mg_df = mg_df[mg_df['disease'].isin(dis_list)]


    # ----- Save the data ----- #
    o_f = f"{out_dir}/prediction_results.txt"
    mg_df.to_csv( o_f, sep = '\t', index=None )


# In[ ]:


if __name__ == '__main__':

    ## Get args.
    args = get_parser()

    ## ----- Merge estimated parameters across diseases ------- ##
    mp = MergeParams(
        args.pert_type,
        args.disease_set,
        args.input_dir,
        args.out_dir
    )
    
    ## Get disease list.
    dis_list = mp.GetDiseaseList()

    ## Merge params.
    mp.Merge(dis_list)
    
    ## New prediction.
    BImethod( args, dis_list )

