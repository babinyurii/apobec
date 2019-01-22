# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 13:20:05 2019

@author: babin
"""

""" count_snp_duplex.py
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
    1. create folder 'input_data' in the current directory and put input 'fasta' files into it
    2. put folder 'apobec' containing scripts into the current directory
    3. run in jupyter : %run ./scripts/count_snp_duplex.py

NOTE: to speed up the exection instead of 'SeqIO.parse()
low level fasta parser 'SimpleFastaParser' is used 
"""


import os
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from time import time
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
        
    return df, coverage - 1


def create_pivot_df(df_snp):
    
    context_data = df_snp.groupby('ref_duplex').sum()
    context_data.drop(['start_pos', 'stop_pos'], inplace=True, axis=1)
    
    return context_data


def save_df(f, df_raw_snp_container, 
            df_duplex_container, 
            df_perc_container,
            df_cov):
    """saves dataframes into an excel spreadsheet
    """
    if not os.path.exists("./output_apobec"):
        os.mkdir("./output_apobec")

    file_name = f.rsplit(".", 1)[0]
    writer = pd.ExcelWriter("./output_apobec/" + file_name +'.xlsx')
    
    counter = 0
    for df in df_raw_snp_container:
        df.to_excel(writer, "raw_count" + str(counter))    
        counter += 1
        
    for df in df_duplex_container:
        tab_name = df.index[0][0]
        df.to_excel(writer, "nuc_" + tab_name)
        
    for df in df_perc_container:
        tab_name = df.index[0][0]
        df.to_excel(writer, "nuc_" + tab_name + "_percentage")
    
    df_cov.to_excel(writer, "coverage")

    writer.save()

    
    
def get_input_files_names(path_to_input):
    
    fasta_extensions = ["fasta", "fa", "fas"]
    input_files = os.listdir(path_to_input)
    input_files = [f for f in input_files if f.rsplit(".", 1)[-1] in fasta_extensions]
     
    return input_files



def create_bar_chart(file_name, df_perc_container, largest_percent):
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

        ax.set_ylabel('percent')
        ax.set_xlabel('reference duplex', fontsize=15)
        ax.set_title("SNP percent at the first position in the duplex context", fontsize=20)
        ax.set_xticks([p + 1.0 * width for p in pos])

        ax.set_xticklabels(df.index)
        plt.legend(loc='upper left', title="SNP type")
        plt.ylim(0, largest_percent + 5)
    
        nuc = df.index[0][0]
        fig.savefig("./output_apobec/" + file_name.rsplit(".", 1)[0] + "_bars_" + nuc +  ".png") 
        # plt.show() # comment to save figure, otherwise it'll save blank file
        plt.close()  # comment the line to show the figure in the jupyter or wherever


def _get_current_time():
    """just returns time stamp
    """
    time_stamp = datetime.datetime.fromtimestamp(
        time()).strftime('%Y-%m-%d %H:%M:%S')
    
    return time_stamp



def show_report(total_time, total_files):
    """prints out  brief report 
    """

    hours = total_time // 3600
    minutes = (total_time % 3600) // 60
    seconds = total_time % 60

    print("""
    
    file processed: {0}
    time taken: {1} hours {2} minutes {3} seconds
    
    the results are in the folder 'output_apobec'
    
    ---------------
    ... job finished at {4}
    ---------------
    
    """.format(total_files, hours, minutes, int(seconds), _get_current_time()))
 


def main():
    
    start_time = time()
    file_counter = 0
    nucleotides = ["A", "T", "G", "C"]
    
    if os.path.exists("./input_data"):
        input_files = get_input_files_names("./input_data")
        
        num_files = len(input_files)
        print("""
               ---------------
               job started at {0} ...
               ---------------
               """.format(_get_current_time()))
        
        progress_bar = FloatProgress(min=0, max=num_files)
        display(progress_bar)
    
        for f in input_files:
            df_raw_snp_container = []
            df_duplex_container = []
            
            for nuc in nucleotides:
                snp_type = nucleotides[:]
                snp_type.remove(nuc)
                
                try:
                    ref_seq = get_ref(f)
                except PermissionError as permerr:
                    print("""
                      warning: It seems that the file'{0}' is open in some other 
                      application, like 'Geneious' or whatever,
                      which doesn't allow it to be processed
                      Now this file is skipped
                       """.format(f))
                
                duplex_posits = get_duplex_posits(ref_seq, nuc)
                df = create_df(duplex_posits, snp_type)
                
                try:
                    df_snp, coverage = count_snp(f, df, nuc)
                except PermissionError as permerr:
                    print("""
                      warning: It seems that the file'{0}' is open in some other 
                      application, like 'Geneious' or whatever,
                      which doesn't allow it to be processed
                      Now this file is skipped
                       """.format(f))
                
                df_duplex_in_context = create_pivot_df(df_snp)
                df_raw_snp_container.append(df_snp)
                df_duplex_container.append(df_duplex_in_context)
                
            # TODO
            # three loops under should be wrapped as function
            # like : get_summary
            
            total_snpes = 0
            for df in df_duplex_container:
                total_snpes += df.sum().sum()
            
            df_perc_container = []
            for df in df_duplex_container:
                df_perc_container.append(df / total_snpes * 100)
                
            largest_percent_in_df = []
            for df in df_perc_container:
                largest_percent_in_df.append(df.max().max())
            largest_percent = max(largest_percent_in_df)
            
            
            df_cov = pd.DataFrame.from_dict({"coverage": coverage}, orient='index')

            try:
                save_df(f, 
                        df_raw_snp_container, 
                        df_duplex_container, 
                        df_perc_container,
                        df_cov)
            except PermissionError as permerr:
                print("""
                      warning: It seems that excel spreadsheet '{0}' is open in excel.
                      To process the corresponding 'fasta' file, rerun the script
                      now this file is skipped. error: {1}
                      """.format(f, permerr))
            
            try:   
                create_bar_chart(f, df_perc_container, largest_percent)
            except ValueError as valerr:
                print(""""
                      plot can't be constructed. 
                      sequences in the input file may be too short. {0}
                      """.format(valerr))
            
            progress_bar.value += 1
            file_counter += 1
        
        finish_time = time()
        total_time = finish_time - start_time
        show_report(total_time, file_counter)
        
    else:
        os.mkdir("./input_data")
        print(
             """
        Houston, we have a problem...
        --------
        folder 'input_data' doesn't exist in the current directory 
        or may be you've created it but misspelled its name.
        Anyway, it has just been created by this script.
        Paste your 'fasta' files into the folder 'input_data' 
        and run this script again.
        --------
        """
        )      
        
if __name__ == "__main__":
    main()


