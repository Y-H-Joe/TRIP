#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 09:45:30 2020

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####
## module 1
Given the NAME column and output_dir, check whether name-folder exists (delete if exists),
otherwise, create it. Create processed_tables folder in working dir to store tables processed from
repeatmaster output. Create barplots folder to store 4 types of barplots generated
from proccessed tables. Create filtered_tables to store filtered tables from
processed tables. Create manual.txt to record failed NAMEs for further manual
procession.

for each thread:
    ## module 2
    Given the BIOPROJECT, download the recording file from ENA database to the
    name-folder. Then use parsing script to extract FASTQ file FTP locations.
    Download the FASTQ files to the name-folder.
    ## module 3
    Given the ASSEMBLY, web-scrap the "Total ungapped length" from NCBI website.
    Write the genome_size to genome_size file in name-folder.
    ## module 4
    Given the dir to repeatmaster (come with TRIP), process the FASTQ files in name-folder
    and generate intermediates into name-folder.
    ## module 5
    Process the repeatmaster output tables and store filtered tables into filtered_tables.
    Store barplots into barplots folder.

## module 6
Given the filtration parameters to filter the processed tables in processed_tables
and generate the TR_candidates folder to store TR candidates.

Add URL column to the log.
Add GENOME_SIZE to the log.
Save TRIP.log.csv

#================================== input =====================================
1.input tsv file:
NAME	BIOPROJECT	ASSEMBLY
CARCR	PRJNA212889	GCA_000690535.1
HARAX	PRJDB6162	GCA_011033045.1

2.path to repeatmaster

2.1 -r	minimal repeat size (default=1)
2.2 -R	maximal repeat size (default=25)
2.3 -n	minimal #repeats*size (default=24, to control false positives)

3.output dir
#================================== output ====================================

#================================ parameters ==================================

#================================== example ===================================

#================================== warning ===================================
system: Linux
The structure of the folder should not change.
"python3" should be abled to call by your system.
needed python pacakges:
    os
    sys
    pandas
    inspect
    numpy
    matplotlib
    requests
    random
    bs4
    time
    traceback
    multiprocessing
    hashlib
####=======================================================================####
"""
import os
import sys
import getopt
import re
import pandas as pd
import inspect
import time
from glob import glob
##from multiprocessing import Process
from multiprocessing import get_context
from multiprocessing import Pool
from multiprocessing import set_start_method
import hashlib
import traceback

#os.system("taskset -p 0xff %d" % os.getpid())

def findnewestgz(file_path):
    filenames = os.listdir(file_path)
    name_ = []
    time_ = []
    #print("newest_gz, filenames: ",filenames)
    for filename in filenames:
        if '.gz' == filename[-3:]: ## only check gz files
            #print(filename)
            c_time = os.path.getctime(file_path+'/'+filename)
            name_.append(file_path+'/'+filename)
            time_.append(c_time)
    #print("newest_gz, name_: ",name_)
    #print(time_)
    if len(time_)>0:
        newest_file = name_[time_.index(max(time_))]
        #print("newest_gz, newest_file: ",newest_file)
        return newest_file
    else:
        #print("newest_gz, 0")
        return 0

def getmd5sum(file_path):
    return hashlib.md5(file_path.encode('utf-8')).hexdigest()

    
def module2to5(name,prj,ass,output_dir,current_script_dir,ENA2URL_loc,wget_loc,skip_subreads,downsample,continue_,url_filter,RepeatDetector_loc,RepeatSummary_loc,RepeatDetector_r,RepeatDetector_R,RepeatDetector_n,CRPG_loc):
    try:
        print("###=============    %s-%s-%s is processing.    ==============###" % (name,prj,ass))
        sys.stdout.flush()
        name_folder_dir=str(output_dir+"/TRIP_results/"+name)  ## created in module 1
        
        ## module 2 command
        module_2_cmd=str("python3 "+current_script_dir+"/module_2.py "+name+" "+\
                         prj+" "+name_folder_dir+" "+ENA2URL_loc+" "+wget_loc+" "+\
                         skip_subreads+" "+downsample+" "+continue_+" "+url_filter)
        print("module 2 cmd: ",module_2_cmd)
        sys.stdout.flush()
        os.system(module_2_cmd)
        ## check md5_log
        ## Check whether new files downloaded. If "no", then no need to re-run the module 3-5
        md5_log_loc=str(name_folder_dir+"/md5_log")
        if os.path.exists(md5_log_loc):
            print(name," md5_log exits.")
            newest_gz=findnewestgz(name_folder_dir)
            try:
                md5sum=getmd5sum(newest_gz)
                #print("md5sum: ",md5sum)
            except:
                print("Unexpected error . Cannot get md5 value of %s. Skip." % (name))
                return
            with open(md5_log_loc,'r') as f:
                #print("reading old md5_log.")
                old_md5sum=str(f.readline())
            if md5sum!=old_md5sum:
                #print("update new md5_log.")
                with open(md5_log_loc,'w') as f:
                    f.write(str(md5sum))
            else:
                print("%s gz files have no change. Skip." % (name))
                return
        else:
            print(name," md5_log doesn't exist.")
            newest_gz=findnewestgz(name_folder_dir)
            #print("newest gz: ",newest_gz)
            if newest_gz==0:
                print("%s has no downloaded gz files, skip." % (name))
                return
            else:
                md5sum=getmd5sum(newest_gz)
                #print("write new md5_log.")
                with open(md5_log_loc,'w') as f:
                    f.write(str(md5sum))

        ## module 3 command
        module_3_cmd=str("python3 "+current_script_dir+"/module_3.py "+ass+" "+name_folder_dir)
        print("module 3 cmd: ",module_3_cmd)
        sys.stdout.flush()
        os.system(module_3_cmd)
        
        ## module 4 command
        RepeatDetector_O=str(name_folder_dir+"/"+name)
        RepeatDetector_I_list=glob(name_folder_dir+"/*gz")
        RepeatDetector_I=" ".join(RepeatDetector_I_list)
        RepeatSummary_O=str(name_folder_dir+"/"+name+"_repeatsummary.tsv")
        RepeatSummary_I=str(name_folder_dir+"/"+name+".repeat")
        module_4_cmd=str("python3 "+current_script_dir+"/module_4.py "+RepeatDetector_loc\
                         +" "+RepeatSummary_loc+" "+RepeatDetector_O+" "+RepeatDetector_r\
                         +" "+RepeatDetector_R+" "+RepeatDetector_n+" "+RepeatDetector_I\
                         +" "+RepeatSummary_O+" "+RepeatSummary_I)
        print("module 4 cmd: ",module_4_cmd)
        sys.stdout.flush()
        os.system(module_4_cmd)
        
        ## module 5 command
        module_5_cmd=str("python3 "+current_script_dir+"/module_5.py "+CRPG_loc+" "+name+" "+name_folder_dir+" "+output_dir)
        print("module 5 cmd: ",module_5_cmd)
        sys.stdout.flush()
        os.system(module_5_cmd)
        
        ## trick to actually acquire cpu cores
        #j=0
        #for i in range(1000000000):
        #    j = i**2
    except Exception as e:
        print("Exception in %s-%s-%s process." % (name,prj,ass))
        traceback.print_exc()
        print()
        raise e
        
if __name__=='__main__':
    set_start_method("spawn")
    
    short_cmd="hi:o:p:r:R:n:"
    long_cmd=["help","input=","output=","process_num=","RepeatDetector_r=","RepeatDetector_R=",\
              "RepeatDetector_n=","rpt_reads_num=","total_reads_num=","unit_len=",\
              "eff_read_len=","genome_size=","avg_genome_cov=","repeats_len=",\
              "repeats_per_read=","reads_per_genome=","repeats_per_genome=",\
              "repeats_per_million_reads=","repeats_len_per_genome=","repeats_len_per_million_reads",\
              "percent_repeats_len_per_read=","percent_repeats_len_per_genome=",\
              "percent_repeat_unit_in_seqs=","best_candidate_enrichment=","max_qualified_num="]
    
    manual=\
        "Program: TRIP (Telomeric Repeat motif Identificatioin Pipeline) \n\n"\
        "Version: 1.0 \n\n"\
        "Contact: Yihang Zhou <yzz0191@auburn.edu> \n\n"\
        "Usage: python3 dir_to_TRIP/TRIP.py -i dir_to_TRIP/example_input.tsv -o output_dir \n\n"\
        "Warning: TRIP use 'python3' as default cmd. Type 'python3' to test whether your system has 'python3'.\n"\
        "         TRIP use 'bash' as the shell. Make sure 'bash' is your default shell.\n"\
        "         Make sure you have the executing authority of all files in TRIP. Try `chmod 777 -R TRIP/` \n\n"\
        "Commands:\n"\
        "-h [--help]                       Print this manual and exit.\n"\
        "-i [--input]                      The input table, which should be tab seprated format with 'NAME','BIOPROJECT','ASSEMBLY' columns. See example_input.tsv.\n"\
        "-o [--output]                     The output dir.\n"\
        "-p [--process_num]                Multi-process mode.The number of processes (default=32).\n"\
        "-r [--RepeatDetector_r]           RepeatDetector:minimal repeat size (default=1). \n"\
        "-R [--RepeatDetector_R]           RepeatDetector:maximal repeat size (default=25). \n"\
        "-n [--RepeatDetector_n]           RepeatDetector:minimal #repeats*size, to control false positives (default=16).\n"\
        "--rpt_reads_num                   Threshold:the number of repeat-containing reads (default=12000).\n"\
        "--total_reads_num                 Threshold:the number of total sequenced reads (default=0).\n"\
        "--repeats_num                     Threshold:the number of repeats (default=0).\n"\
        "--unit_len                        Threshold:the length of a repeat unit (bp), unincluding. (default=5).\n"\
        "--eff_read_len                    Threshold:the effective length of reads (bp), unincluding. (default=0).\n"\
        "--genome_size                     Threshold:genome assembly length (bp), unincluding. (default=0).\n"\
        "--avg_genome_cov                  Threshold:the average coverage of a genome, unincluding. (default=10).\n"\
        "--repeats_len                     Threshold:the length of repeats (bp), unincluding. (default=0).\n"\
        "--repeats_per_read                Threshold:the average number of repeats per read in repeats containing reads, unincluding. (default=0).\n"\
        "--reads_per_genome                Threshold:the number of repeats containing reads per genome, unincluding. (default=0).\n"\
        "--repeats_per_genome              Threshold:number of repeats per haploid genome, unincluding. (default=0).\n"\
        "--repeats_per_million_reads       Threshold:number of repeats per 1 million repeats containing reads, unincluding. (default=0).\n"\
        "--repeats_len_per_genome          Threshold:total repeat length (Kb) in a haploid genome, unincluding. (default=0).\n"\
        "--repeats_len_per_million_reads   Threshold:the average repeat length (Kb) per million reads, unincluding. (default=4).\n"\
        "--percent_repeats_len_per_read    Threshold:the average percentage of  repeat length in repeat-containing reads, unincluding. (default=0.5).\n"\
        "--percent_repeats_len_per_genome  Threshold:the percentage of total repeats length in a haploid genome, unincluding. (default=0).\n"\
        "--percent_repeat_unit_in_seqs     Threshold:the average percentage of repeat units length in all sequencing reads, unincluding. (default=0). \n"\
        "--best_candidate_enrichment       Call criterion:the length ratio of the most abundant repeat  over the next abundant repeat, including. (default=3).\n"\
        "--max_qualified_num               Call criterion:the maximum number of identified TR candidates from upstream, including. (default=999).\n"\
        "--skip_subreads                   Skip subreads. The subreads are sequenced from Pacbio. (default=True). [True, False].\n"\
        "--downsample                      Downsample number of sequencing files. (default=4). [integer, 0]. 0 means sample all.\n"\
        "--continue                        Continue last interrupted downloads. (defaults=True). [True, False].\n"\
        "--url_filter                      Overwrite, append or remain the url_filtered file, which decides how/whether to download new files. (default=remain). [overwrite,append,remain]."
    
    
    try:
        opts,args=getopt.getopt(sys.argv[1:],short_cmd,long_cmd)
    except getopt.GetoptError as err:
        print(manual)
        print(str(err))
        sys.exit(2)
    
    
    infile_dir="/TRIP/example_input.tsv"
    output_dir="."
    
    current_script_dir=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    ## str to int
    process_num=str(32)
    ENA2URL_loc=str(current_script_dir+"/ENA2URL.sh")
    wget_loc=str(current_script_dir+"/wget.sh")
    RepeatDetector_loc=str(current_script_dir+"/RepeatDetector_v2")
    RepeatSummary_loc=str(current_script_dir+"/RepeatSummary_v2")
    RepeatDetector_r=str(1)
    RepeatDetector_R=str(25)
    RepeatDetector_n=str(16)
    CRPG_loc=str(current_script_dir+"/cal_repeats_per_genome_and_percent_of_len_linux.py")
    rpt_reads_num=str(12000)
    total_reads_num=str(0)
    repeats_num=str(0)
    unit_len=str(4)
    eff_read_len=str(0)
    genome_size=str(0)
    avg_genome_cov=str(10)
    repeats_len=str(0)
    repeats_per_read=str(0)
    reads_per_genome=str(0)
    repeats_per_genome=str(0)
    repeats_per_million_reads=str(0)
    repeats_len_per_genome=str(0)
    repeats_len_per_million_reads=str(4) ##Kb
    percent_repeats_len_per_read=str(0.5)
    percent_repeats_len_per_genome=str(0)
    percent_repeat_unit_in_seqs=str(0)
    best_candidate_enrichment=str(3)
    max_qualified_num=str(999)
    skip_subreads=str("True")
    downsample=str(4)
    continue_=str("True")
    url_filter=str("remain")
    FOT2X_loc=str(current_script_dir+"/filter_output_tables_to_xlsx.py")
    
    for o,a in opts:
        if o in ("-h","--help"):
            print(manual)
            sys.exit()
        elif o in ("-i","--input"):
            infile_dir=a
        elif o in ("-o","--output"):
            output_dir=a
        elif o in ("-p","--process_num"):
            process_num=a
        elif o in ("-r","--RepeatDetector_r"):
            RepeatDetector_r=a
        elif o in ("-R","--RepeatDetector_R"):
            RepeatDetector_R=a
        elif o in ("-n","--RepeatDetector_n"):
            RepeatDetector_n=a
        elif o == "--rpt_reads_num":
            rpt_reads_num=a
        elif o == "--total_reads_num":
            total_reads_num=a
        elif o == "--unit_len":
            unit_len=a
        elif o == "--eff_read_len":
            eff_read_len=a
        elif o == "--genome_size":
            genome_size=a
        elif o == "--avg_genome_cov":
            avg_genome_cov=a
        elif o == "--repeats_len":
            repeats_len=a
        elif o == "--repeats_per_read":
            repeats_per_read=a
        elif o == "--reads_per_genome":
            reads_per_genome=a
        elif o == "--repeats_per_genome":
            repeats_per_genome=a
        elif o == "--repeats_per_million_reads":
            repeats_per_million_reads=a
        elif o == "--repeats_len_per_genome":
            repeats_len_per_genome=a
        elif o == "--repeats_len_per_million_reads":
            repeats_len_per_million_reads=a
        elif o == "--percent_repeats_len_per_read":
            percent_repeats_len_per_read=a
        elif o == "--percent_repeats_len_per_genome":
            percent_repeats_len_per_genome=a
        elif o == "--percent_repeat_unit_in_seqs":
            percent_repeat_unit_in_seqs=a
        elif o == "--best_candidate_enrichment":
            best_candidate_enrichment=a
        elif o == "--max_qualified_num":
            max_qualified_num=a
        elif o == "--skip_subreads":
            skip_subreads=a
        elif o == "--downsample":
            downsample=a
        elif o == "--continue":
            continue_=a
        elif o == "--url_filter":
            url_filter=a
        else:
            print(manual)
            sys.exit()
    
    if len(opts)==0:
        print(manual)
        sys.exit()

    ## read infile
    try:
        infile_df=pd.read_csv(infile_dir,header=0,sep="\t")
    except Exception as e:
        print("error: The input tsv file has issues.")
        print(e)
        sys.exit()

    ## get the list of BIOPROJECTs
    name_list=list(infile_df.loc[:,'NAME'])
    if len(name_list)!=len(set(name_list)):
        print("error: repeat NAME. exit.")
        sys.exit()
    prj_list_temp=list(infile_df.loc[:,'BIOPROJECT'])
    prj_list=[]
    for prj in prj_list_temp:
        prj_split=re.split("[,|:|;| |\\|\|]",prj)
        if len(prj_split)>1:
            print("warning: multiple BIOPROJECT recordings in one species, use the first one. prj: ",prj)
        prj_list.append(prj_split[0])

    ass_list=list(infile_df.loc[:,'ASSEMBLY'])
    name_prj_ass_list=list(zip(name_list,prj_list,ass_list))
    name_num=len(name_list)

    ## module 1 command
    module_1_cmd=str("python3 "+current_script_dir+"/module_1.py "+infile_dir+" "+output_dir)
    print("module 1 cmd: ",module_1_cmd)
    os.system(module_1_cmd)

    ## multi-process
    process_num = len(name_prj_ass_list) if len(name_prj_ass_list) < int(process_num) else int(process_num)
    print("TRIP is using {} process.".format(int(process_num)))

    pool=Pool(int(process_num))
    count=0
    for name_prj_ass in name_prj_ass_list:
        name=name_prj_ass[0]
        prj=name_prj_ass[1]
        ass=name_prj_ass[2]
        sys.stdout.flush()
        time.sleep(5)  ## in case NCBI API restriction
        
        count+=1
        #x1 = pool.apply_async(long_time_task, args=(count,current_script_dir,))
        #x1.get()
        x2 = pool.apply_async(module2to5, args=(name,prj,ass,output_dir,current_script_dir,ENA2URL_loc,wget_loc,skip_subreads,downsample,continue_,url_filter,RepeatDetector_loc,RepeatSummary_loc,RepeatDetector_r,RepeatDetector_R,RepeatDetector_n,CRPG_loc,))
        #x2.get()
        #print("###=============    %s-%s-%s is module 2 to 5 finished.    ==============###" % (name,prj,ass))

    pool.close()
    time.sleep(5)
    pool.join()

    ## module 6 command
    filtered_tables_dir=str(output_dir+"/TRIP_results/"+"filtered_tables")
    TR_candidates_dir=str(output_dir+"/TRIP_results/"+"TR_candidates")
    module_6_cmd=str("python3 "+current_script_dir+"/module_6.py "+filtered_tables_dir+" "+\
                    rpt_reads_num+" "+\
                    total_reads_num+" "+\
                    repeats_num+" "+\
                    unit_len+" "+\
                    eff_read_len+" "+\
                    genome_size+" "+\
                    avg_genome_cov+" "+\
                    repeats_len+" "+\
                    repeats_per_read+" "+\
                    reads_per_genome+" "+\
                    repeats_per_genome+" "+\
                    repeats_per_million_reads+" "+\
                    repeats_len_per_genome+" "+\
                    repeats_len_per_million_reads+" "+\
                    percent_repeats_len_per_read+" "+\
                    percent_repeats_len_per_genome+" "+\
                    percent_repeat_unit_in_seqs+" "+\
                    best_candidate_enrichment+" "+\
                    max_qualified_num+" "+\
                    TR_candidates_dir+" "+\
                    FOT2X_loc)
    print("module 6 cmd: ",module_6_cmd)
    sys.stdout.flush()
    os.system(module_6_cmd)
    ## Add URL column to the log.
    ## Add GENOME_SIZE to the log.
    URL_list=[]
    GENOME_SIZE_list=[]
    for name in name_list:
        url_loc=str(output_dir+"/TRIP_results/"+name+"/url")
        GENOME_SIZE_loc=str(output_dir+"/TRIP_results/"+name+"/genome_size")
        with open(url_loc,'r') as f:
            lines=[x.strip() for x in f.readlines()]
            URL=";".join(lines)
            URL_list.append(URL)
        with open(GENOME_SIZE_loc,'r') as f:
            GENOME_SIZE=int(f.readline())
            GENOME_SIZE_list.append(GENOME_SIZE)

    infile_df['URL']=URL_list
    infile_df['GENOME_SIZE']=GENOME_SIZE_list

    infile_df.to_csv(str(output_dir+"/TRIP_results/"+"TRIP.log.csv"),index=None)

