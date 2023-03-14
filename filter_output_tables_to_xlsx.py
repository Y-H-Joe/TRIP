#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 22:56:01 2020

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ discription ==============================####
#
#================================== input =====================================
#
#================================== output ====================================
#
#================================ parameters ==================================
#
#================================== example ===================================
#
#================================== warning ===================================
#
####=======================================================================####
"""

import pandas as pd
import os
import traceback
import sys
#os.system('taskset -p %s' %os.getpid())

filtered_tables_dir=sys.argv[1]
rpt_reads_num_filter=float(sys.argv[2])
total_reads_num_filter=float(sys.argv[3])
repeats_num_filter=float(sys.argv[4])
unit_len_filter=float(sys.argv[5])
eff_read_len_filter=float(sys.argv[6])
genome_size_filter=float(sys.argv[7])
avg_genome_cov_filter=float(sys.argv[8])
repeats_len_filter=float(sys.argv[9])
repeats_per_read_filter=float(sys.argv[10])
reads_per_genome_filter=float(sys.argv[11])
repeats_per_genome_filter=float(sys.argv[12])
repeats_per_million_reads_filter=float(sys.argv[13])
repeats_len_per_genome_filter=float(sys.argv[14])
repeats_len_per_million_reads_filter=float(sys.argv[15])
percent_repeats_len_per_read_filter=float(sys.argv[16])
percent_repeats_len_per_genome_filter=float(sys.argv[17])
percent_repeat_unit_in_seqs_filter=float(sys.argv[18])
best_candidate_enrichment=float(sys.argv[19])
max_qualified_num_filter=float(sys.argv[20])
TR_candidates_dir=sys.argv[21]


processed_tables=filtered_tables_dir.rstrip("/").rstrip("\\")
table_names=os.listdir(filtered_tables_dir)


## folder for .xlsx format tables with availability to filter
if not os.path.exists(TR_candidates_dir):
    os.makedirs(TR_candidates_dir)

## set up rules.
## rules are based on column names defined in previous python script.
def filter_out_low_quality_candidates(df,rpt_reads_num_filter,\
total_reads_num_filter,repeats_num_filter,unit_len_filter,eff_read_len_filter,genome_size_filter,\
avg_genome_cov_filter,repeats_len_filter,repeats_per_read_filter,reads_per_genome_filter,repeats_per_genome_filter,\
repeats_per_million_reads_filter,repeats_len_per_genome_filter,repeats_len_per_million_reads_filter,\
percent_repeats_len_per_read_filter,percent_repeats_len_per_genome_filter,\
percent_repeat_unit_in_seqs_filter):
    
    ## Recommend filters
    ## data quality guarantee 
    avg_genome_cov=df["avg_genome_cov"].to_numpy()
    repeat_containing_reads_num=df["reads"].to_numpy()
    ## telomere objective law
    ## 1. make sure we sequenced telomere region
    repeats_len_per_million_reads=df["repeats_len_per_million_reads"].to_numpy()
    ## 2. make sure the reads mapped to the real continuous telomere region
    ## exclude the unknown nosiy sparse/discrete repeats region
    percent_repeats_len_per_read=df["percent_repeats_len_per_read"].to_numpy()
    ## 3. TR unit length less than 5bp are excluded
    unit_len=df["unit_len"].to_numpy()
    ## Recommend filters
    avg_genome_cov_scores=[1 if x>avg_genome_cov_filter else 0 for x in avg_genome_cov] ## 10
    repeat_containing_reads_scores=[1 if x>rpt_reads_num_filter else 0 for x in repeat_containing_reads_num ] ## 12000
    ## repeats_len is kb
    repeats_len_per_million_reads_scores=[1 if x > repeats_len_per_million_reads_filter else 0 for x in repeats_len_per_million_reads] # 4
    percent_repeats_len_per_read_scores=[1 if x > percent_repeats_len_per_read_filter else 0 for x in percent_repeats_len_per_read] # 0.5
    unit_len_scores=[1 if x>unit_len_filter else 0 for x in unit_len ] # 4
    
    
    ## Optional filters
    total_reads_num=df["num_reads"].to_numpy()
    repeats_num=df["repeats"].to_numpy()
    eff_read_len=df["eff_read_len"].to_numpy()
    genome_size=df["genome_size"].to_numpy()
    ## repeats_per_read can exclude 2bp repeats, such as AG, repeats region noise
    repeats_per_read=df["repeats_per_read"].to_numpy()
    reads_per_genome=df["reads_per_genome"].to_numpy()
    repeats_per_genome=df["repeats_per_genome"].to_numpy()
    repeats_per_million_reads=df["repeats_per_million_reads"].to_numpy()
    repeats_len_per_genome=df["repeats_len_per_genome"].to_numpy()
    percent_repeats_len_per_genome=df["percent_repeats_len_per_genome"].to_numpy()
    percent_repeat_unit_in_seqs=df["percent_repeat_unit_in_sequences"].to_numpy()
    ## The "dominant" principle.
    ## Only accpet situations which no more than 2 TR candidates appeared at 
    ## the same time. if there are too many qualified TR candidates, we can not guarantee 
    ## it is not due to retrotransposon, sequencing errors or other unexpected issues.
    ## In order to identify the real TR motif, we objectively use
    ## repeats_length. The dominant one must exceed ${best_candidate_enrichment}
    ## times than the second one as for repeats_length. Otherwise, both of them will
    ## be excluded.
    repeats_len=df['repeats_len'].to_numpy()
    ## Optional filters
    total_reads_num_scores=[1 if x>total_reads_num_filter else 0 for x in total_reads_num]
    repeats_num_scores=[1 if x>repeats_num_filter else 0 for x in repeats_num]
    eff_read_len_scores=[1 if x>eff_read_len_filter else 0 for x in eff_read_len]
    genome_size_scores=[1 if x>genome_size_filter else 0 for x in genome_size]
    repeats_per_read_scores=[1 if x>repeats_per_read_filter else 0 for x in repeats_per_read]
    reads_per_genome_scores=[1 if x>reads_per_genome_filter else 0 for x in reads_per_genome]
    repeats_per_genome_scores=[1 if x>repeats_per_genome_filter else 0 for x in repeats_per_genome]
    repeats_per_million_reads_scores=[1 if x>repeats_per_million_reads_filter else 0 for x in repeats_per_million_reads]
    repeats_len_per_genome_scores=[1 if x>repeats_len_per_genome_filter else 0 for x in repeats_len_per_genome]
    percent_repeats_len_per_genome_scores=[1 if x>percent_repeats_len_per_genome_filter else 0 for x in percent_repeats_len_per_genome]
    percent_repeat_unit_in_seqs_scores=[1 if x>percent_repeat_unit_in_seqs_filter else 0 for x in percent_repeat_unit_in_seqs]
    repeats_len_scores=[1 if x>repeats_len_filter else 0 for x in repeats_len]
    
    
    
    sum_scores=[sum(i) for i in zip(
        avg_genome_cov_scores,\
        repeat_containing_reads_scores,\
        repeats_len_per_million_reads_scores,\
        percent_repeats_len_per_read_scores,\
        unit_len_scores,\
        total_reads_num_scores,\
        repeats_num_scores,\
        eff_read_len_scores,\
        genome_size_scores,\
        repeats_per_read_scores,\
        reads_per_genome_scores,\
        repeats_per_genome_scores,\
        repeats_per_million_reads_scores,\
        repeats_len_per_genome_scores,\
        percent_repeats_len_per_genome_scores,\
        percent_repeat_unit_in_seqs_scores,\
        repeats_len_scores
        )]
    ##print(sum_scores)
    qualified=[1 if x==17 else 0 for x in sum_scores]

    return avg_genome_cov_scores,\
        repeat_containing_reads_scores,\
        repeats_len_per_million_reads_scores,\
        percent_repeats_len_per_read_scores,\
        unit_len_scores,\
        total_reads_num_scores,\
        repeats_num_scores,\
        eff_read_len_scores,\
        genome_size_scores,\
        repeats_per_read_scores,\
        reads_per_genome_scores,\
        repeats_per_genome_scores,\
        repeats_per_million_reads_scores,\
        repeats_len_per_genome_scores,\
        percent_repeats_len_per_genome_scores,\
        percent_repeat_unit_in_seqs_scores,\
        repeats_len_scores,\
        qualified




## all possible condidates
## if there're more than 1 candidates, which is ambiguous
TR_candidates_qualified_df_sum=pd.DataFrame()
TR_candidates_qualified_dominant=pd.DataFrame()

for name in table_names:
    print(name)
    if "csv" in name or "tsv" in name: ## filter unexpected files
        processed_table="/".join([filtered_tables_dir,name])
        processed_table_df=pd.read_csv(processed_table)
        try:
            avg_genome_cov_scores,\
            repeat_containing_reads_scores,\
            repeats_len_per_million_reads_scores,\
            percent_repeats_len_per_read_scores,\
            unit_len_scores,\
            total_reads_num_scores,\
            repeats_num_scores,\
            eff_read_len_scores,\
            genome_size_scores,\
            repeats_per_read_scores,\
            reads_per_genome_scores,\
            repeats_per_genome_scores,\
            repeats_per_million_reads_scores,\
            repeats_len_per_genome_scores,\
            percent_repeats_len_per_genome_scores,\
            percent_repeat_unit_in_seqs_scores,\
            repeats_len_scores,\
            qualified=filter_out_low_quality_candidates(processed_table_df,rpt_reads_num_filter,\
                    total_reads_num_filter,repeats_num_filter,unit_len_filter,eff_read_len_filter,genome_size_filter,\
                    avg_genome_cov_filter,repeats_len_filter,repeats_per_read_filter,reads_per_genome_filter,repeats_per_genome_filter,\
                    repeats_per_million_reads_filter,repeats_len_per_genome_filter,repeats_len_per_million_reads_filter,\
                    percent_repeats_len_per_read_filter,percent_repeats_len_per_genome_filter,\
                    percent_repeat_unit_in_seqs_filter)
            #repeats_per_read_filter,\
        except :
            print(name," has problem:")
            traceback.print_exc()
            continue

        filter_results=pd.DataFrame({
            'avg_genome_cov_scores':avg_genome_cov_scores,\
            'repeat_containing_reads_scores':repeat_containing_reads_scores,\
            'repeats_len_per_million_reads_scores':repeats_len_per_million_reads_scores,\
            'percent_repeats_len_per_read_scores':percent_repeats_len_per_read_scores,\
            'unit_len_scores':unit_len_scores,\
            'total_reads_num_scores':total_reads_num_scores,\
            'repeats_num_scores':repeats_num_scores,\
            'eff_read_len_scores':eff_read_len_scores,\
            'genome_size_scores':genome_size_scores,\
            'repeats_per_read_scores':repeats_per_read_scores,\
            'reads_per_genome_scores':reads_per_genome_scores,\
            'repeats_per_genome_scores':repeats_per_genome_scores,\
            'repeats_per_million_reads_scores':repeats_per_million_reads_scores,\
            'repeats_len_per_genome_scores':repeats_len_per_genome_scores,\
            'percent_repeats_len_per_genome_scores':percent_repeats_len_per_genome_scores,\
            'percent_repeat_unit_in_seqs_scores':percent_repeat_unit_in_seqs_scores,\
            'repeats_len_scores':repeats_len_scores,\
            'qualified':qualified
                })

        TR_candidates_df=pd.concat([processed_table_df,filter_results],axis=1)

        TR_candidates_qualified_df=TR_candidates_df[TR_candidates_df['qualified']==1].copy()
        TR_candidates_qualified_num=TR_candidates_qualified_df.shape[0]
        TR_candidates_qualified_df.index=range(TR_candidates_qualified_num)
        
        ## TR_candidates_qualified_df need to be formatted to output TR_candidates_qualified_df_sum
        TR_candidates_qualified_df.loc[:,'species']=[name]*TR_candidates_qualified_num
        TR_candidates_qualified_df.loc[:,'qualified_num']=[TR_candidates_qualified_num]*TR_candidates_qualified_num
        TR_candidates_qualified_df_sum=pd.concat([TR_candidates_qualified_df_sum,TR_candidates_qualified_df],axis=0)
        
        ## TR_candidates_df need to be formatted to output TR_candidates_qualified_dominant
        TR_candidates_df.to_excel(str(TR_candidates_dir+"/"+name+"."+str(TR_candidates_qualified_num)+".filtered.xlsx"),index=None)
        
        ## there cannot be to many qualified TR candidates, otherwise it will be too noisy
        max_qualified_num=max_qualified_num_filter

        ## the TR candidate should satisfy two conditions
        ## one, it should be qualified (${qualified}==1)
        ## two, it should be  dominant (repeats_len is more than ${difference_threold} times than the follower)
        difference_threhold=best_candidate_enrichment
        if max_qualified_num>=TR_candidates_qualified_num>0:
            ## the TR candidate must extremly more doinant than the followers
            ## so set the repeats_len_per_million_reads difference to 5 times
            sorted_TR_qualified_candidates_df=\
            TR_candidates_qualified_df.sort_values("repeats_len",ascending=False,inplace=False).copy()
            sorted_TR_qualified_candidates_df.index=range(sorted_TR_qualified_candidates_df.shape[0])
            count=0
            
            if TR_candidates_qualified_num==1:
                dominant_TR_candidate=pd.DataFrame(sorted_TR_qualified_candidates_df.loc[0]).T.copy()
                dominant_TR_candidate.loc[:,'difference']=1
                TR_candidates_qualified_dominant=pd.concat([TR_candidates_qualified_dominant,dominant_TR_candidate],axis=0)
                
                print(name," : only dominant")                
                
            elif TR_candidates_qualified_num>1:
                dominant=1
                first_repeats_len=sorted_TR_qualified_candidates_df.loc[0,'repeats_len']
                second_repeats_len=sorted_TR_qualified_candidates_df.loc[1,'repeats_len']
                difference=first_repeats_len/second_repeats_len
                
                print(name," : ",difference)
                      
                if difference>difference_threhold:
                    dominant=1
                    dominant_TR_candidate=pd.DataFrame(sorted_TR_qualified_candidates_df.loc[0]).T
                    dominant_TR_candidate.loc[:,'difference']=difference
                    TR_candidates_qualified_dominant=pd.concat([TR_candidates_qualified_dominant,dominant_TR_candidate],axis=0)

        else:
            print(name," has no qualified TR candidate, or there're too many of them.")
            continue


    else:
        continue

TR_candidates_qualified_df_sum.to_excel(str(TR_candidates_dir+"/"+"TR_candidates.qualified.xlsx"),index=None)
TR_candidates_qualified_dominant.to_excel(str(TR_candidates_dir+"/"+"TR_candidates.qualified.dominant.xlsx"),index=None)

