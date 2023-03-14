#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 13:15:45 2020

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####
## module 2
Given the BIOPROJECT, download the recording file from ENA database to the
name-folder. Then use parsing script to extract FASTQ file FTP locations.
Download the FASTQ files to the name-folder.
#================================== input =====================================

#================================== output ====================================

#================================ parameters ==================================

#================================== example ===================================
python3 current_script_dir/module_2.py name prj name_folder_dir ENA2URL_loc wget_loc
#================================== warning ===================================

####=======================================================================####
"""
import sys
import os
import time
import random
#os.system('taskset -p %s' %os.getpid())

name=sys.argv[1]
prj=sys.argv[2]
name_folder_dir=sys.argv[3]
ENA2URL_loc=sys.argv[4]
wget_loc=sys.argv[5]
skip_subreads=sys.argv[6]
downsample=int(sys.argv[7])
continue_=sys.argv[8]
url_filter=sys.argv[9]

## download records from ENA
ENA_record=str(name+"_"+prj+"_tsv.txt")
print("module 2: downloading ",ENA_record)
url=str("\"https://www.ebi.ac.uk/ena/portal/api/filereport?accession="+prj+"&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true\"")
wget_cmd=str("wget --retry-connrefused -t 10  --timeout=300 -q -c "+url+" -O "+name_folder_dir+"/"+ENA_record)
os.system(wget_cmd)

## parse the FTP locations
cp_cmd=str("cp -t "+name_folder_dir+" "+ENA2URL_loc+" "+wget_loc)
os.system(cp_cmd)
#print("module 2: ",cp_cmd)
## will generate url file under name-folder
parse_cmd=str("bash "+name_folder_dir+"/ENA2URL.sh")
os.system(parse_cmd)
#print("module 2: ",parse_cmd)

## filter the url file and generate url_filtered file under name-folder
if skip_subreads=="True":
    skip_subreads=True
elif skip_subreads=="False":
    skip_subreads=False
else:
    print("module 2: skip_subreads=",skip_subreads," is not in [True,False], change it to True")
    skip_subreads=True
url_file_dir=str(name_folder_dir+"/url")
url_filtered_file_dir=str(name_folder_dir+"/url_filtered")
with open(url_file_dir,'r') as f:
    urls=f.readlines()
    subreads_filtered_urls=[]
    downsample_filtered_urls=[]
    filtered_urls=[]
    # subreads filter
    if skip_subreads:
        for line in urls:
            if "subreads" not in line:
                subreads_filtered_urls.append(line)
    # downsample filter
    if len(subreads_filtered_urls)>0:
        downsample_filtered_urls=subreads_filtered_urls
    else:
        downsample_filtered_urls=urls
    if 0<downsample<len(downsample_filtered_urls):
        filtered_urls=random.sample(downsample_filtered_urls,downsample)
    else:
        filtered_urls=downsample_filtered_urls
        
    if continue_=="True":
        if os.path.isfile(url_filtered_file_dir): ## 
            if url_filter=="remain":
                pass
            elif url_filter=="overwrite":
                with open(url_filtered_file_dir,'w') as ff:
                    for line in filtered_urls:
                        ff.write(line)
            elif url_filter=="append":
                with open(url_filtered_file_dir,'a') as ff:
                    for line in filtered_urls:
                        ff.write(line)
            else:
                print("module 2: url_filter parameter error, need to be in [overwrite,append,remain]: ",url_filter)
        else:
            # no original url_filtered_file_dir, write
            with open(url_filtered_file_dir,'w') as ff:
                for line in filtered_urls:
                    ff.write(line)


## Download the FASTQ files to the name-folder
wget_cmd=str("bash "+name_folder_dir+"/wget.sh")
time.sleep(random.choice(range(2,20))) ## in case IP restriction
#print(wget_cmd)
os.system(wget_cmd)


