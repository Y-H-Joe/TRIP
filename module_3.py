#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 16:04:46 2020

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####
## module 3 
Given the ASSEMBLY, web-scrap the "Total ungapped length" from NCBI website. 
Write the genome_size to genome_size file in name-folder.
#================================== input =====================================

#================================== output ====================================

#================================ parameters ==================================

#================================== example ===================================

#================================== warning ===================================

####=======================================================================####
"""
import requests
import random
import time
from bs4 import BeautifulSoup
import sys
import os
#os.system('taskset -p %s' %os.getpid())

assembly=sys.argv[1]
name_folder_dir=sys.argv[2]

time.sleep(random.choice(range(1,99))) ## in case IP restriction

## scrap the genome_size=Total ungapped length
url=str("https://www.ncbi.nlm.nih.gov/assembly/"+assembly+"/")
print("module 3: requesting ",url)
html=requests.get(url)
bf=BeautifulSoup(html.text,features="lxml")
len_str=bf.find("td",string="Total ungapped length").next_sibling.text
len_int=int(len_str.replace(",",""))

## write the genome size
print("module 3: writing genome_size.")
with open(str(name_folder_dir+"/genome_size"),'w') as f:
    f.write(str(len_int))
