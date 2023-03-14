#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 14:45:48 2021

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####
## module 4
Given the dir to repeatmaster (come with TRIP), process the FASTQ files in name-folder
and generate intermediates into name-folder.
#================================== input =====================================

#================================== output ====================================

#================================ parameters ==================================

#================================== example ===================================
RepeatDetector  name_out  name/*gz  -r 1 -R 25 -n 24
RepeatSummary name_25nt.tsv name_out.repeat
#================================== warning ===================================

####=======================================================================####
"""
import sys
import os
#os.system('taskset -p %s' %os.getpid())

RepeatDetector_loc=sys.argv[1]
RepeatSummary_loc=sys.argv[2]
RepeatDetector_O=sys.argv[3]
RepeatDetector_r=sys.argv[4]
RepeatDetector_R=sys.argv[5]
RepeatDetector_n=sys.argv[6]
RepeatDetector_I=" ".join(sys.argv[7:-2])
RepeatSummary_O=sys.argv[-2]
RepeatSummary_I=sys.argv[-1]


RepeatDetector_cmd=str(RepeatDetector_loc+" "+RepeatDetector_O+" "+RepeatDetector_I\
                       +" -r "+RepeatDetector_r+" -R "+RepeatDetector_R+" -n "+\
                       RepeatDetector_n)
print("module 4: RepeatDetector_cmd: ",RepeatDetector_cmd)
os.system(RepeatDetector_cmd)

#print("RepeatSummary_loc: ",RepeatSummary_loc)
#print("RepeatSummary_O: ",RepeatSummary_O)
#print("RepeatSummary_I: ",RepeatSummary_I)
RepeatSummary_cmd=str(RepeatSummary_loc+" "+RepeatSummary_O+" "+RepeatSummary_I)
print("module 4: RepeatSummary_cmd: ",RepeatSummary_cmd)
os.system(RepeatSummary_cmd)
