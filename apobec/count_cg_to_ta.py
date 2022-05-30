# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 13:20:05 2019

@author: babin
"""
# TODO redo the comment according to c > t and g > a change
# TODO redo all the comments
""" count_snp_duplex.py
this script count exact number of SNPs and their percentage according
to the dinucleotide context in "fasta" files
It writes results into  excel spreadsheets.
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
NOTE: to speed up the exection instead of 'SeqIO.parse()
low level fasta parser 'SimpleFastaParser' is used 
"""

import datetime
import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from Bio.SeqIO.FastaIO import SimpleFastaParser  # low level fast fasta parser
from ipywidgets import IntProgress
from IPython.display import display
from time import time
from loguru import logger
sns.set()


def get_ref_name(in_file):
    """deriving the first record 
    from the input fasta
    """

    with open("./input_data/" + in_file) as in_handle:
        for record in SimpleFastaParser(in_handle):
            # sequence is under the 2nd index in the record's tuple
            ref_name = record[0]
            break

    return ref_name

def get_input_files_names(path_to_input):
    """returns a list of input fasta files
    """
    fasta_extensions = ["fasta", "fa", "fas"]
    input_files = os.listdir(path_to_input)
    input_files = [f for f in input_files if f.rsplit(
        ".", 1)[-1] in fasta_extensions]

    return input_files


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


def get_cg_posits_in_ref(ref_seq):
    """
    """

    d = {'C_pos': [], 'G_pos': []}

    for nuc_index in range(0, len(ref_seq)):
        if ref_seq[nuc_index] == "C":
            d["C_pos"].append(nuc_index)

        elif ref_seq[nuc_index] == "G":
            d["G_pos"].append(nuc_index)

    return d


def get_num_of_reads_in_file(in_file):
    """
    """
    reads_counter = 0
    with open("./input_data/" + in_file) as in_handle:
        for record in SimpleFastaParser(in_handle):
            reads_counter += 1

    return reads_counter - 1

def count_num_of_c_and_g_in_ref(cg_posits):
    """
    """
    num_c_in_ref = len(cg_posits["C_pos"])
    num_g_in_ref = len(cg_posits["G_pos"])

    return num_c_in_ref, num_g_in_ref


def count_total_c_and_g_in_reads(ref_name, cg_posits, in_file):

    c_counter = 0
    g_counter = 0
    nucs = ["A", "T", "G", "C"]

    with open("./input_data/" + in_file) as in_handle:
        for record in SimpleFastaParser(in_handle):
           
            if record[0] == ref_name:
                continue
            for c_index in cg_posits["C_pos"]:
                if record[1][c_index] in nucs:
                    c_counter += 1

            for g_index in cg_posits["G_pos"]:
                if record[1][g_index] in nucs:
                    g_counter += 1

    return c_counter, g_counter


#def count_total_c_and_g_in_reads(num_c_in_ref,
#                                num_g_in_ref,
 #                               num_reads_in_file):
 #   """
#    """
#    total_num_c_in_reads = num_c_in_ref * num_reads_in_file
 #   total_num_g_in_reads = num_g_in_ref * num_reads_in_file
#
#    return total_num_c_in_reads, total_num_g_in_reads


def count_c_t_change(f, cg_posits):
    """
    """
    c_t_change_counter = 0
    with open("./input_data/" + f) as in_handle:
        for record in SimpleFastaParser(in_handle):
            for c_index in cg_posits["C_pos"]:

                if record[1][c_index] == "T":
                    c_t_change_counter += 1

    return c_t_change_counter 



def count_g_a_change(f, cg_posits):
    """
    """
    g_a_change_counter = 0
    with open("./input_data/" + f) as in_handle:
        for record in SimpleFastaParser(in_handle):
            for g_index in cg_posits["G_pos"]:

                if record[1][g_index] == "A":
                    g_a_change_counter += 1

    return g_a_change_counter 


def calculate_percent_c_t_change(num_c_t_change_in_reads, 
                                    total_num_c_in_reads, 
                                    num_c_in_ref):
    """
    """
    if num_c_in_ref: # check if there are Cs nucs in ref
        percent_c_t_change = (num_c_t_change_in_reads/total_num_c_in_reads) * 100
    else:
        percent_c_t_change = "no C in reference"

    return percent_c_t_change


def calculate_percent_g_a_change(num_g_a_change_in_reads, 
                                    total_num_g_in_reads,
                                    num_g_in_ref):
    """
    """
    if num_g_in_ref: # check if there are Gs nucs in ref
        percent_g_a_change = (num_g_a_change_in_reads/total_num_g_in_reads) * 100
    else:
        percent_g_a_change = "no G in reference"

    return percent_g_a_change


def calculate_c_t_and_g_a_change_percent(total_num_c_in_reads,
                                            total_num_g_in_reads,
                                            num_c_t_change_in_reads,
                                            num_g_a_change_in_reads):
    """
    """                                        
    total_gc = total_num_c_in_reads + total_num_g_in_reads
    total_cg_ta_change_in_reads = num_c_t_change_in_reads + num_g_a_change_in_reads

    c_t_and_g_a_change_percent = (total_cg_ta_change_in_reads/total_gc) * 100

    return c_t_and_g_a_change_percent


def create_df(num_reads_in_file, 
                total_num_c_in_reads, 
                total_num_g_in_reads, 
                num_c_t_change_in_reads,
                num_g_a_change_in_reads, 
                percent_c_t_change, 
                percent_g_a_change,
                c_t_and_g_a_change_percent):
    """
    """

    data = {"number of reads" : [num_reads_in_file],
            "C in reads count" : [total_num_c_in_reads],
            "G in reads count" : [total_num_g_in_reads],
            "C>T change count" : [num_c_t_change_in_reads],
            "G>A change count" : [num_g_a_change_in_reads],
            "C>T change percent" : [percent_c_t_change],
            "G>A change percent" : [percent_g_a_change],
            "C>T + G>A change percent" : [c_t_and_g_a_change_percent]}

    df = pd.DataFrame(data)

    return df


def save_df(f, df):
    """saves dataframes into an excel spreadsheet
    """
    if not os.path.exists("./output_count_cg_to_ta_change"):
        os.mkdir("./output_count_cg_to_ta_change")

    file_name = f.rsplit(".", 1)[0]
    writer = pd.ExcelWriter("./output_count_cg_to_ta_change/" + file_name + '.xlsx')

    df.to_excel(writer, "cg_to_ta_change")

    writer.save()


def get_current_time():
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
    
    """.format(total_files, hours, minutes, int(seconds), get_current_time()))




#TODO put messages into separate file, into a dict, for example
# and import it to just call the items to show them
# so as to make main function not so large

# TODO delete all the try inside the main
# or just test them prior to it




def main():

    start_time = time()
    file_counter = 0
    #############################
    # logging
    # TODO make function for logging
    if not os.path.exists("./log"):
        os.mkdir("./log")
    logger.remove() # don't put messages into notebook output
    logger.add("./log/apobec_count_cg_to_ta_log_{time}" + "_" + ".txt", backtrace=False)
    logger.add(sys.stderr, level="CRITICAL")
    ##############################
    

    # TODO: make function of if/else
    if os.path.exists("./input_data"): #TODO make if not and then return
    # to get rid of else block below
        
        input_files = get_input_files_names("./input_data")
        num_files = len(input_files)
        
        #############################
        # a progress bar 
        print("""
               ---------------
               job started at {0} ...
               ---------------
               """.format(get_current_time()))

        progress_bar = IntProgress(min=0, max=num_files, bar_style='success')
        display(progress_bar)
        ##################################
        

        for f in input_files:
            try:
                ref_seq = get_ref(f) 
                ref_name = get_ref_name(f)
                cg_posits = get_cg_posits_in_ref(ref_seq)
                num_reads_in_file = get_num_of_reads_in_file(f)
                num_c_in_ref, num_g_in_ref = count_num_of_c_and_g_in_ref(cg_posits)
                total_num_c_in_reads, total_num_g_in_reads = count_total_c_and_g_in_reads(ref_name, cg_posits, f)

                # C>T and G>A change                                                                                        
                num_c_t_change_in_reads = count_c_t_change(f, cg_posits)
                num_g_a_change_in_reads = count_g_a_change(f, cg_posits)
                
                percent_c_t_change = calculate_percent_c_t_change(num_c_t_change_in_reads, 
                                                                        total_num_c_in_reads, 
                                                                        num_c_in_ref)
               
                percent_g_a_change = calculate_percent_g_a_change(num_g_a_change_in_reads, 
                                                                        total_num_g_in_reads,
                                                                        num_g_in_ref)
    
                c_t_and_g_a_change_percent = calculate_c_t_and_g_a_change_percent(total_num_c_in_reads,
                                                                                    total_num_g_in_reads,
                                                                                    num_c_t_change_in_reads,
                                                                                    num_g_a_change_in_reads)

                df = create_df(num_reads_in_file, 
                                total_num_c_in_reads, 
                                total_num_g_in_reads, 
                                num_c_t_change_in_reads,
                                num_g_a_change_in_reads, 
                                percent_c_t_change, 
                                percent_g_a_change,
                                c_t_and_g_a_change_percent)

                save_df(f, df)

                
                
               
                #df_cov = pd.DataFrame.from_dict(
                #    {"coverage": coverage}, orient='index')
                
        
                
                ##############################
                # again progress bar
                progress_bar.value += 1
                file_counter += 1
            
            except Exception:
                logger.info(
                    """\n\t==================================================
                    file: {0} 
                    sequence id: {1}
                    ====================================================""".format(f, record_id))
                logger.exception("")
                print("exception detected. see log for more details")
                progress_bar.value += 1                

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
