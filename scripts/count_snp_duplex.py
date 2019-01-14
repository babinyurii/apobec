# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 13:20:05 2019

@author: babin
"""

""" apobec_count_duplex.py
this script count exact number of SNPs and their percentage in "fasta" file
and writes results into the excel spreadsheet.
the script compares each sequence with the reference
and counts SNPs in the context of duplexes
it takes 'fasta' file imported from Geneious as an input
'fasta' file must be imported from Geneious with the following parameters:
    - wrap sequence lines every 80 chars: No
    - replace spaces in sequence names with undescores: Yes
    - include sequence description: Yes
    - export sequence in: upper case
    - export missing ends of alignments as: Ns
IMPORTANT NOTES:
    - first record in fasta file must be your reference sequence
    - reference sequence must start and end with nucleotides, not gaps
    f.e.: 'A-A---A---A' is the right string, 
    but: '----A---AAAA---' is not valid
    if reference sequence starts or ends with gap, the scrips will throw: 
    'ValueError: Length mismatch: Expected axis has N elements, 
    new values have n elements'. This is the pandas exception concerning indices.

to use the script:
    1. create folder 'input_data' in the current directory and paste 'fasta' files into it
    2. import script into the current directory using Jupyter or whatever and run it:
        1) from apobec_scripts import apobec_count_duplex
        execute apobec_count_duplex()
        2) %run apobec_count_duplex.py
    or put the script into the current directory and run via shell:
    python run_count_indels.py

NOTE: to speed up the exection instead of 'SeqIO.parse()
low level fasta parser 'SimpleFastaParser' is used 
"""


import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser  # low level fast fasta parser
from ipywidgets import FloatProgress
from IPython.display import display
sns.set()


def get_ref(in_file):
    """deriving the first record 
    from the input fasta
    """
    
    with open("./input_data/" + in_file) as in_handle:
        for record in SimpleFastaParser(in_handle):
            # sequence is under the 2nd index in the record's tuple
            ref_seq = record[1] 
            break
    
    return ref_seq
    


def get_duplex_posits(ref_seq, ref_nuc):
    """takes duplex positions out of the reference string
    """
    nucleotides = ['A', 'T', 'G', 'C']
    start = 0
    
    d = {'start_pos': [], 'stop_pos': [], 'ref_duplex': [] }
    
    while start != len(ref_seq) - 1:
        
        if ref_seq[start] == ref_nuc:
            first_duplex_nuc = ref_seq[start]
            first_duplex_pos = start
            
            # + 1 to get slice from the next nuc
            second_duplex_pos = start + 1
            for nuc in ref_seq[start + 1 : ]:
                if nuc in nucleotides:
                    second_duplex_nuc = nuc
                    break
                else:
                    second_duplex_pos += 1
                        
            d['start_pos'].append(first_duplex_pos)
            d['stop_pos'].append(second_duplex_pos)
            d['ref_duplex'].append(first_duplex_nuc + second_duplex_nuc)
            # incrementing index to move forward the string
            start += 1
    
        else:
            # incrementing anyway
            start += 1
    
    return d

def create_df(duplex_posits, snp_type):

    df = pd.DataFrame(duplex_posits)
    for snp in snp_type:
        df[snp] = np.nan
    
    df.fillna(inplace=True, value=0)
    
    return df
    

def count_snp(in_file, df, nuc_to_watch):
    coverage = 0
    
    with open("./input_data/" + in_file) as in_handle:
        
        reads = SimpleFastaParser(in_handle)
        for record in reads:
            coverage += 1
            
            read = record[1]
            for row in df.index:
                if read[df.loc[row, 'start_pos']] != nuc_to_watch and read[df.loc[row, 'start_pos']] != '-' \
                and read[df.loc[row, 'stop_pos']] == df.loc[row, 'ref_duplex'][1]:
                    
                    snp = read[df.loc[row, 'start_pos']]
                    
                    df.loc[row, snp] += 1
        
    return df, coverage


def create_pivot_df(df_snp):
    
    context_data = df_snp.groupby('ref_duplex').sum()
    context_data.drop(['start_pos', 'stop_pos'], inplace=True, axis=1)
    
    return context_data


def save_df(f, df_container, df_perc_container):
    """saves dataframes into an excel spreadsheet
    """
    if not os.path.exists("./output_apobec"):
        os.mkdir("output_apobec")

    file_name = f.rsplit(".", 1)[0]
    writer = pd.ExcelWriter("./output_apobec/" + file_name +'.xlsx')
    
    for df in df_container:
        tab_name = df.index[0][0]
        df.to_excel(writer, "nuc_" + tab_name)
        
    for df in df_perc_container:
        tab_name = df.index[0][0]
        df.to_excel(writer, "nuc_" + tab_name + "_percentage")

    writer.save()

    
    
def get_input_files_names():
    path_to_input = "./input_data"
    fasta_extensions = ["fasta", "fa", "fas"]
    input_files = os.listdir("./input_data")
    input_files = [f for f in input_files if f.rsplit(".", 1)[-1] in fasta_extensions]
     
    return input_files



def create_bar_chart(file_name, df_perc_container):
    """
    creates bar chart using list of dataframes with 
    snp percentage. Percentage is calculated relatively 
    total number of snp found
    """
    for df in df_perc_container:

        pos = [0, 1, 2, 3]
        width = 0.25 
        
        fig, ax = plt.subplots(figsize=(10,5))
    
        columns = df.columns.tolist()

        plt.bar(pos, 
                df[columns[0]], 
                width, 
                alpha=0.5, 
                color='darkcyan', 
                edgecolor='black',
                lw=2,
                label=columns[0])

        plt.bar([p + width for p in pos], 
                df[columns[1]],
                width, 
                alpha=0.5, 
                color='silver', 
                edgecolor='black',
                lw=2,
                label=columns[1])

        plt.bar([p + width * 2 for p in pos], 
                df[columns[2]], 
                width, 
                alpha=0.5, 
                color='sandybrown', 
                edgecolor='black',
                lw=2,
                label=columns[2])

        ax.set_ylabel('raw count')
        ax.set_xlabel('duplex', fontsize=15)
        ax.set_title("number of SNPs in the context", fontsize=20)
        ax.set_xticks([p + 1.0 * width for p in pos])

        ax.set_xticklabels(df.index)
        plt.legend(loc='upper left')
        plt.ylim(0, 30)
    
        nuc = df.index[0][0]
        fig.savefig("./output_apobec/" + file_name.rsplit(".", 1)[0] + "_bars_" + nuc +  ".png") 


def main():
    
    nucleotides = ["A", "T", "G", "C"]
    
    if os.path.exists("./input_data"):
        input_files = get_input_files_names()
        
        num_files = len(input_files)
        
        print("job started at : ")
        progress_bar = FloatProgress(min=0, max=num_files)
        display(progress_bar)
    
        for f in input_files:
            #print("processing file '{}'".format(f))
            
            df_container = []
            for nuc in nucleotides:
        
                snp_type = nucleotides[:]
                snp_type.remove(nuc)
    
                ref_seq = get_ref(f)
    
                duplex_posits = get_duplex_posits(ref_seq, nuc)
    
                df = create_df(duplex_posits, snp_type)
    
                df_snp, coverage = count_snp(f, df, nuc)
                df_duplex_context = create_pivot_df(df_snp)
                df_container.append(df_duplex_context)

            # counting total number of variants detected
            # and calculating percentage pivot tables
            total_snpes = 0
            for df in df_container:
                total_snpes += df.sum().sum()
            
            df_perc_container = []
            for df in df_container:
                df_perc_container.append(df / total_snpes * 100)
            
            save_df(f,  df_container, df_perc_container)
            
            create_bar_chart(f, df_perc_container)

            progress_bar.value += 1
            
    else:
        os.mkdir("input_data")        
        
        
if __name__ == "__main__":
    main()





