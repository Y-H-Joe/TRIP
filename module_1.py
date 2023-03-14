#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 10:24:18 2020

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####
## module 1
Given the NAME column and output_dir, check whether name-folder exists (use it if exists),
otherwise, create it. Createprocessed_tables folder in working dir to store tables processed from
repeatmaster output. Create barplots folder to store 4 types of barplots generated
from proccessed tables. Create filtered_tables to store filtered tables from 
processed tables. Create manual.txt to record failed NAMEs for further manual
procession.
#================================== input =====================================

#================================== output ====================================

#================================ parameters ==================================

#================================== example ===================================
python3 module_1.py infile_dir output_dir
#================================== warning ===================================

####=======================================================================####
"""
import sys
import os
import pandas as pd
#os.system('taskset -p %s' %os.getpid())

infile_dir=sys.argv[1]
output_dir=sys.argv[2].rstrip("\\").rstrip("/")
output_dir=str(output_dir+"/"+"TRIP_results")

## read infile
try:
    infile_df=pd.read_csv(infile_dir,header=0,sep="\t")
except Exception as e:
    print("module 1 error: The input tsv file has issues.")
    print(e)
    sys.exit()
    
## get the list of NAMEs
name_list=list(infile_df.loc[:,'NAME'])

## check whether the parent folder exists
try:
    if os.path.exists(output_dir):
        print("module 1: ",output_dir," exists. Use it.")
    else:
        print("module 1: ",output_dir," doesn't exist. Create it.")
        os.makedirs(output_dir)
except Exception as e:
    print("module 1 error: Creating TR_results error.")
    print(e)
    sys.exit()

## check whether the name folder exists
try:
    for name in name_list:
        name_folder_dir=str(output_dir+"/"+name)
        if os.path.exists(name_folder_dir):
            print("module 1: ",name_folder_dir," exists. Use it.")
        else:
            print("module 1: ",name_folder_dir," doesn't exist. Create it.")
            os.makedirs(name_folder_dir)
except Exception as e:
    print("module 1 error: Creating name folder error.")
    print(e)
    sys.exit()
    

## check whether filtered_tables folder exists
filtered_tables_dir=str(output_dir+"/"+"filtered_tables")
if os.path.exists(filtered_tables_dir):
    print("module 1: ",filtered_tables_dir," exists. Use it.")
else:
    print("module 1: ",filtered_tables_dir," doesn't exist. Create it.")
    os.makedirs(filtered_tables_dir)


## check whether barplots folder exists
barplots_dir=str(output_dir+"/"+"barplots")
if os.path.exists(barplots_dir):
    print("module 1: ",barplots_dir," exists. Use it.")
else:
    print("module 1: ",barplots_dir," doesn't exist. Create it.")
    os.makedirs(barplots_dir)
    


