# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 14:33:49 2018

@author: Yihang Zhou
"""
"""
CRPG_input.csv:
code,read_len,num_reads,avg_genome_cov,genome_size
SUSSC,151,693665606,42.36,2472461935
"""
def format_output(data,eff_read_len,num_reads,avg_genome_cov,genome_size,output_dir):
    """
    format data
    check data integrity
    """
    if len(data.columns)<4:
        if 'unit_len' not in data.columns:
            list_unit_len=[len(x) for x in data.index]
            df_unit_len=pd.DataFrame(list_unit_len)
            df_unit_len.index=data.index
            df_unit_len.columns=['unit_len']
            data=pd.concat([data,df_unit_len],axis=1)
        if 'repeats_per_read' not in data.columns:
            list_rpr=data['repeats']/data['reads']
            df_rpr=pd.DataFrame(list_rpr)
            df_rpr.columns=['repeats_per_read']
            data=pd.concat([data,df_rpr],axis=1)

    ## remove single bases A,C,G,T
    ## remove repeats_per_read<4
    for index in data.index:
        if data.loc[index,'repeats_per_read']<4:
            data=data.drop(index)
            continue
        if len(index)==1:
            data=data.drop(index)
            continue


    ## known variables
    #reads
    #repeats
    #unit_len
    #avg_genome_cov
    #genome_size
    #eff_read_len
    #num_reads
    reads=np.array(data.loc[:,'reads'])
    repeats=np.array(data.loc[:,'repeats'])
    unit_len=np.array(data.loc[:,'unit_len'])


    ## calculated variables
    #repeats_len------------------> the length of repeats
    #repeats_per_read-------------------> average number of repeats per read in repeats containing reads
    #reads_per_genome-------------------> repeat_containing_reads_num per genome
    #repeats_per_genome-----------------> number of repeats per haploid genome
    #repeats_per_million_reads----------> number of repeats per 1 million sequencing reads
    #repeats_len_per_genome-------------> average total repeat length in haploid genome
    #repeats_len_per_million_reads------> average total repeat length per 1 million sequencing reads
    #percent_repeats_len_per_read-------> average percentage of perfect repetitive unit region in repeat containing reads
    #percent_repeats_len_per_genome-----> average percentage of total repeats length in haploid genome
    #percent_repeat_unit_in_sequences --> average percentage of repeat units in all sequencing reads
    repeats_len=repeats*unit_len.tolist()
    ## repeats_per_read=list(map(int,repeats/reads))
    repeats_per_read=list(repeats/reads)
    reads_per_genome=list(reads/avg_genome_cov)
    ## repeats_per_genome=list(map(int,repeats/avg_genome_cov))
    repeats_per_genome=list(repeats/avg_genome_cov)
    ## repeats_per_million_reads=list(map(int,repeats/num_reads*1000000))
    repeats_per_million_reads=list(repeats/num_reads*1000000)
    ## divide by 1000 to have kb
    ## repeats_len_per_genome=list(map(int,repeats*unit_len/avg_genome_cov/1000))
    repeats_len_per_genome=list(repeats*unit_len/avg_genome_cov/1000)
    ## divide by 1000 to have kb
    ## repeats_len_per_million_reads=list(map(int,repeats*unit_len/num_reads*1000000/1000))
    repeats_len_per_million_reads=list(repeats*unit_len/num_reads*1000000/1000)
    percent_repeats_len_per_read=list(repeats*unit_len/reads/eff_read_len)
    percent_repeats_len_per_genome=list(np.array([x*1000 for x in repeats_len_per_genome])/genome_size)
    percent_repeat_unit_in_sequences=list(repeats*unit_len/num_reads/eff_read_len)

    ## output
    length=len(data.index)
    output=pd.DataFrame({'reads':reads,\
                         'repeats':repeats,\
                         'unit_len':unit_len,\
                         'avg_genome_cov':[avg_genome_cov for x in range(length)],\
                         'genome_size':[genome_size for x in range(length)],\
                         'eff_read_len':[eff_read_len for x in range(length)],\
                         'num_reads':[num_reads for x in range(length)],\
                         'repeats_len':repeats_len,\
                         'repeats_per_read':repeats_per_read,\
                         'reads_per_genome':reads_per_genome,\
                         'repeats_per_genome':repeats_per_genome,\
                         'repeats_per_million_reads':repeats_per_million_reads,\
                         'repeats_len_per_genome':repeats_len_per_genome,\
                         'repeats_len_per_million_reads':repeats_len_per_million_reads,\
                         'percent_repeats_len_per_read':percent_repeats_len_per_read,\
                         'percent_repeats_len_per_genome':percent_repeats_len_per_genome,\
                         'percent_repeat_unit_in_sequences':percent_repeat_unit_in_sequences
                         },list(data.index))
    ## output.to_csv(str("/Users/yihangzhou/CurrentProjects/Telomere/data/repeat_unit_table_after_filtered/"+species+"_output.csv"))
    output.to_csv(str(output_dir+"/"+species+"_output.csv"))
    return output



## matplotlib
## ylabel parameter here is useless
def matplotlib_fig(df,column,index=None,ylabel='',yh_max=1,yh_min=0,log=False,former=30,format_='pdf',dp=""):
    """
    if index!=None, then order the repeats using df[index] value, extract top
    $former repeats, to get the desired top $former df.index, and replace the
    repeats with df[column], and finally sort the dictionary based on index
    length again.
    """
    if index!=None:
        ## dict_ is a dictionary with the key of df.index and value of df[column]
        dict_=dict(zip(list(df.index),list(df.loc[:,column])))
        ## dict_number is a dictionary with the key of df.index and value of df[index]
        dict_number=dict(zip(list(df.index),list(df.loc[:,index])))
        sorted_dict_number=sorted(dict_number.items(),key=lambda kv:kv[1],reverse=True)[:former]
        ## length of keys=former
        keys=[item[0] for item in sorted_dict_number]
        #filtered_dict_=dict_.fromkeys(keys)
        ## filtered_dict_ is the key:value pair of df.index:df[column],whose key is
        ## from the index of top $former df[$index]
        filtered_dict_={k:v for k,v in dict_.items() if k in keys}
        sorted_dict_=sorted(filtered_dict_.items(),key=lambda kv:kv[1],reverse=True)
    else:
        dict_=dict(zip(list(df.index),list(df.loc[:,column])))
        sorted_dict_=sorted(dict_.items(),key=lambda kv:kv[1],reverse=True)[:former]
    ## sorted_dict_2 is the sorted-based-on-the-length-of-index('AG','ATG'...)
    ## version of sorted_dict_
    sorted_dict_2=[sorted_dict_[0]]
    for i in range(former)[1:]:
        length=len(sorted_dict_[i][0])
        if length>=len(sorted_dict_2[-1][0]):
            sorted_dict_2.append(sorted_dict_[i])
            continue
        if length<len(sorted_dict_2[0][0]):
            sorted_dict_2.insert(0,sorted_dict_[i])
            continue
        for j in range(len(sorted_dict_2)):
            length_1=len(sorted_dict_2[j][0])
            length_2=len(sorted_dict_2[j+1][0])
            if length_1<=length<length_2:
                sorted_dict_2.insert(j+1,sorted_dict_[i])
                break
    ## summarize the data and prepare to plot
    repeats=[]
    number=[]
    number_right=[]
    for item in sorted_dict_2:
        repeats.append(item[0])
        number.append(item[1])
    for repeat in repeats:
        number_right.append(df.loc[repeat,'repeats_per_read'])
    data_=pd.DataFrame({'number':number,\
                       'number_right':number_right},repeats)


    ## bar-plots
    width=0.8
    ylabel_right='Average # of repeats per read'
    x_pos=list(range(len(repeats)))
    fig=plt.figure(figsize=(12,4))
    ax1=fig.add_subplot(111)

    ax1.bar(x_pos,height=data_.number,width=width,color='r',log=log,alpha=0.5,tick_label=repeats)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(90)
    ax1.set_ylim([yh_min,yh_max])
    ax1.set_ylabel(column)

    ax2=ax1.twinx()
    ax2.scatter(x=x_pos,y=data_.number_right,color='b',alpha=0.5)
    ax2.set_ylabel(ylabel_right)
    ax2.set_yticks([2*x for x in range(12)])
    plt.title(column)
    #plt.savefig(str('/Users/yihangzhou/CurrentProjects/telomere/figures/repeat_unit_barplots/'+species+'_'+column+"."+format_),format=format_,dpi=600,bbox_inches='tight')
    plt.savefig(str(dp+"/"+species+'_'+column+"."+format_),format=format_,dpi=600,bbox_inches='tight')
    #plt.show()
    plt.close("all")

## output
#matplotlib_fig(output,"repeats_per_genome",log=True,ylabel='Normalized number log10')
#matplotlib_fig(output,"percent_repeats_per_genome")
#matplotlib_fig(output,"repeats_len_per_genome",log=True,ylabel='Normalized length(kb) log10',index='repeats_per_genome')
#matplotlib_fig(output,"percent_repeats_len_per_genome",index='repeats_per_genome')
#matplotlib_fig(output,"repeats_per_million_reads",log=True,ylabel='Normalized number log10')
#matplotlib_fig(output,"percent_repeats_per_million_reads")
#matplotlib_fig(output,"repeats_len_per_million_reads",log=True,ylabel='Normalized length(kb) log10',index='repeats_per_million_reads')
#matplotlib_fig(output,"percent_repeats_len_per_million_reads",index='repeats_per_million_reads')

if __name__=='__main__':
    import pandas as pd
    import numpy as np
    import matplotlib
    ## Agg is the non-interactive mode, so the script won't display fig in command-line linux system
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import sys
    import os
    #os.system('taskset -p %s' %os.getpid())

    ## import data
    RepeatMaster_output=sys.argv[1]
    ## RepeatMaster_output="/analysis1/yihang_analysis/telomere_0602_2020/data/repeat_unit_table_and_status"
    ## RepeatMaster_output="/Users/yihangzhou/CurrentProjects/telomere/data/repeat_unit_table_and_status/"
    RepeatMaster_output=RepeatMaster_output.rstrip("/")
    input_table=sys.argv[2]
    output_tables=sys.argv[3]
    output_figs=sys.argv[4]
    ## input_table="CRPG_input.tsv"
    input_table_df=pd.read_csv(input_table,sep="\t")

    ## folder for output tables
    try:
        if not os.path.exists(output_tables):
            os.makedirs(output_tables)
        ## folder for figs
        if not os.path.exists(output_figs):
            os.makedirs(output_figs)
    except:
        pass

    for index in input_table_df.index:
        species=input_table_df.loc[index,"code"]
        eff_read_len=input_table_df.loc[index,"read_len"]
        num_reads=input_table_df.loc[index,"num_reads"]
        avg_genome_cov=input_table_df.loc[index,"avg_genome_cov"]
        genome_size=input_table_df.loc[index,"genome_size"]

        file_type="tsv"
        if file_type=="csv":
            sep=","
        if file_type=="tsv":
            sep="\t"
        try:
            data=pd.read_csv(str(RepeatMaster_output+"/"+species+"_repeatsummary."+file_type),index_col=0,sep=sep) ### need to modify
        except Exception as e:
            print("Can't find ",str(RepeatMaster_output+"/"+species+"_repeatsummary"+file_type)," to read.")
            print(e)
            continue
        print(species," is being processed by cal_repeats_per_genome_and_percent_of_len_linux.py" )
        output=format_output(data,eff_read_len,num_reads,avg_genome_cov,genome_size,output_tables)
        try:
            matplotlib_fig(output,"percent_repeats_len_per_read",log=False,ylabel='Percentage',index='percent_repeat_unit_in_sequences',dp=output_figs)
            matplotlib_fig(output,"percent_repeat_unit_in_sequences",yh_max=0.04,log=True,ylabel='Percentage',index='percent_repeat_unit_in_sequences',dp=output_figs)
            matplotlib_fig(output,"repeats_len_per_genome",log=True,yh_max=12000,ylabel='Normalized length (kb)',index='repeats_per_genome',dp=output_figs)
            matplotlib_fig(output,"repeats_len_per_million_reads",log=True,yh_max=5000,ylabel='Normalized length (kb)',index='repeats_per_million_reads',dp=output_figs)
        except Exception as e:
            print(species," has plot problem. Maybe it's for lacking of available TR canidates.")
            print(e)
         

