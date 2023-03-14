#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 17:01:06 2021

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####
## module 6
Given the filtration parameters to filter the processed tables in processed_tables
and generate the TR_candidates table and TR_candidates.dominant table.
#================================== input =====================================

#================================== output ====================================

#================================ parameters ==================================

#================================== example ===================================

#================================== warning ===================================

####=======================================================================####
"""
import sys
import os
#os.system('taskset -p %s' %os.getpid())

filtered_tables_dir=sys.argv[1]
rpt_reads_num_filter=sys.argv[2]
total_reads_num_filter=sys.argv[3]
repeats_num_filter=sys.argv[4]
unit_len_filter=sys.argv[5]
eff_read_len_filter=sys.argv[6]
genome_size_filter=sys.argv[7]
avg_genome_cov_filter=sys.argv[8]
repeats_len_filter=sys.argv[9]
repeats_per_read_filter=sys.argv[10]
reads_per_genome_filter=sys.argv[11]
repeats_per_genome_filter=sys.argv[12]
repeats_per_million_reads_filter=sys.argv[13]
repeats_len_per_genome_filter=sys.argv[14]
repeats_len_per_million_reads_filter=sys.argv[15]
percent_repeats_len_per_read_filter=sys.argv[16]
percent_repeats_len_per_genome_filter=sys.argv[17]
percent_repeat_unit_in_seqs_filter=sys.argv[18]
best_candidate_enrichment=sys.argv[19]
max_qualified_num_filter=sys.argv[20]
TR_candidates_dir=sys.argv[21]
FOT2X_loc=sys.argv[22]


## Given the filtration parameters to filter the processed tables in processed_tables
## and generate the TR_candidates table and TR_candidates.dominant table.
FOT2X_cmd=str("python3 "+FOT2X_loc+" "+filtered_tables_dir+" "+\
                str(rpt_reads_num_filter)+" "+\
                str(total_reads_num_filter)+" "+\
                str(repeats_num_filter)+" "+\
                str(unit_len_filter)+" "+\
                str(eff_read_len_filter)+" "+\
                str(genome_size_filter)+" "+\
                str(avg_genome_cov_filter)+" "+\
                str(repeats_len_filter)+" "+\
                str(repeats_per_read_filter)+" "+\
                str(reads_per_genome_filter)+" "+\
                str(repeats_per_genome_filter)+" "+\
                str(repeats_per_million_reads_filter)+" "+\
                str(repeats_len_per_genome_filter)+" "+\
                str(repeats_len_per_million_reads_filter)+" "+\
                str(percent_repeats_len_per_read_filter)+" "+\
                str(percent_repeats_len_per_genome_filter)+" "+\
                str(percent_repeat_unit_in_seqs_filter)+" "+\
                str(best_candidate_enrichment)+" "+\
                str(max_qualified_num_filter)+" "+\
                TR_candidates_dir)
print("module 6: ",FOT2X_cmd)
os.system(FOT2X_cmd)



