# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 13:20:05 2019

@author: babin
"""

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


def get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, ref_nuc):
    """takes duplex positions out of the reference string, 
    uses only one reference nuc in one call and finds
    all its occurences in the reference string
    """
    nucleotides = ['A', 'T', 'G', 'C']
    start = 0

    d = {'start_pos': [], 'stop_pos': [], 'ref_duplex': []}

    while start != len(ref_seq) - 1:

        if ref_seq[start] == ref_nuc:
            first_duplex_nuc = ref_seq[start] # nuc from the ref at pos 'start'
            first_duplex_pos = start # position of the 

            # + 1 to get slice from the next nuc
            second_duplex_pos = start + 1
            # taking slice to the end of the ref_seq to find the second nuc
            # in the duplex (as there may be gaps)
            for nuc in ref_seq[start + 1 : ]:
                if nuc in nucleotides:
                    second_duplex_nuc = nuc
                    break
                # if there's no a nuc (a gap ('-') for example),
                # incrementing second_duplex_nuc position
                # and going to the next nuc 
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

def get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, ref_nuc):
    """ 
    var names are all vice verca here compared to 
    get_duplex_posits_ref_nuc_at_first_pos, as we count here
    second nuc in duplex
    -------
    takes duplex positions out of the reference string, 
    uses only one reference nuc in one call and finds
    all its occurences in the reference string
    """
    nucleotides = ['A', 'T', 'G', 'C']
    start = 1 # CHANGED : start from the second char of the string as
    # we count second mismatch to reference in duplex context
    # so there may be a mismatch of the first nuc at [0] index

    d = {'start_pos': [], 'stop_pos': [], 'ref_duplex': []}

    while start != len(ref_seq): # CHANGED : to grasp the penultimate nuc mismatch

        if ref_seq[start] == ref_nuc:
            second_duplex_nuc = ref_seq[start] # nuc from the ref at pos 'start'
            second_duplex_pos = start # position of the 

            # + 1 to get slice from the next nuc
            first_duplex_pos = start - 1
            # taking slice to the START of the ref_seq to find the second nuc
            # in the duplex (as there may be gaps)
            for nuc in ref_seq[start - 1 : : -1 ]:
                if nuc in nucleotides:
                    first_duplex_nuc = nuc
                    break
                # if there's no a nuc (a gap ('-') for example),
                # DECREMENTING second_duplex_nuc position
                # and going to the next nuc TOWARDS THE START
                else:
                    first_duplex_pos -= 1 # CHANGED : GO TO START OF THE REF SEQ

            
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
    """constructs DataFrame 
    columns in the df will be: 
        - start_pos, stop_pos, ref_duplex, and three snp cols
        according to snp_type parameter
    rows will be pandas indices starting from 0
    an example:
       | start_pos | stop_pos | ref_duplex | T | G | C |
    ------------------------------------------------------
     0 |   1       |    2     |      AT    | 2 | 37 | 8 |
    """
    
    df = pd.DataFrame(duplex_posits)
    for snp in snp_type:
        df[snp] = np.nan

    df.fillna(inplace=True, value=0)

    return df


def count_snp_with_ref_nuc_at_first_pos(in_file, df, nuc_to_watch):
    """counts SNPs which occur in the first
    positions of reference duplexes. SNP is counted 
    only when the second positions in a duplex in a read
    remains as in the reference
    """
    coverage = 0

    with open("./input_data/" + in_file) as in_handle:

        reads = SimpleFastaParser(in_handle)
        for record in reads:
            coverage += 1
            read_id = record[0]

            read = record[1]
            for row in df.index:
                # three conditions must be to judge that there's a snp (mismatch between ref and read at the same pos):
                # we search here certain position in the current read
                # 1. nuc at the first duplex pos in read (start_pos in df) is not the nuc we currently use as a reference. 
                # this is the snp. because if they are the same the duplexes are the same (we currently search for
                # mismatch at the first position)
                # 2. nuc at the first duplex pos in read (start_pos in df and ref) is not a gap (-)
                # 3. nuc at the 2nd pos of the duplex in read (stop_pos in df and ref) is the same as the 
                # nuc at the reference duplex at the 2nd position: here we check that the context isn't changed 
                if read[df.loc[row, 'start_pos']] != nuc_to_watch and read[df.loc[row, 'start_pos']] != '-' \
                        and read[df.loc[row, 'stop_pos']] == df.loc[row, 'ref_duplex'][1]:
                    # if there's a mismatch between start_pos in ref and start_pos in nuc (1st pos of duplex)
                    # we take nuc from the read to find the corresponding column in df and increment the number of snps
                    snp = read[df.loc[row, 'start_pos']]

                    df.loc[row, snp] += 1

    return df, coverage - 1, read_id


def count_snp_with_ref_nuc_at_second_pos(in_file, df, nuc_to_watch):
    """counts SNPs which occur in the first
    positions of reference duplexes. SNP is counted 
    only when the second positions in a duplex in a read
    remains as in the reference
    """
    coverage = 0

    with open("./input_data/" + in_file) as in_handle:

        reads = SimpleFastaParser(in_handle)
        for record in reads:
            coverage += 1
            read_id = record[0]

            read = record[1]
            for row in df.index:
                # three conditions must be to judge that there's a snp (mismatch between ref and read at the same pos):
                # we search here certain position in the current read
                # 1. nuc at the SECOND duplex pos in read (stop in df) is not the nuc we currently use as a reference. 
                # this is the snp. because if they are the same the duplexes are the same (we currently search for
                # mismatch at the SECOND position)
                # 2. nuc at the SECOND duplex pos in read (stop_pos in df and ref) is not a gap (-)
                # 3. here we check that the context isn't changed : nuc at the 1st pos of the duplex in a read 
                # (start_pos in df and ref) is the same as the nuc at the reference duplex at the 1st position: 
                if read[df.loc[row, 'stop_pos']] != nuc_to_watch and read[df.loc[row, 'stop_pos']] != '-' \
                        and read[df.loc[row, 'start_pos']] == df.loc[row, 'ref_duplex'][0]:
                    # if there's a mismatch between start_pos in ref and start_pos in nuc (1st pos of duplex)
                    # we take nuc from the read to find the corresponding column in df and increment the number of snps
                    snp = read[df.loc[row, 'stop_pos']]

                    df.loc[row, snp] += 1

    return df, coverage - 1, read_id


def create_pivot_df(df_snp):
    """constructs pivot DataFrame
     grouping data by the duplex variant
     """
    context_data = df_snp.groupby('ref_duplex').sum()
    context_data.drop(['start_pos', 'stop_pos'], inplace=True, axis=1)
    #print(context_data)
    return context_data


def save_df_first_in_duplex(f,
                            df_raw_snp_container,
                            df_duplex_container,
                            df_perc_container,
                            df_cov,
                            nuc_pos_in_duplex_to_count):
    """saves dataframes into an excel spreadsheet
    """
    if not os.path.exists("./output_apobec"):
        os.mkdir("./output_apobec")

    file_name = f.rsplit(".", 1)[0]
    writer = pd.ExcelWriter("./output_apobec/" + file_name + "_" + nuc_pos_in_duplex_to_count + '.xlsx')

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

def save_df_second_in_duplex(f, 
                            df_raw_snp_container,
                            df_duplex_container,
                            df_perc_container,
                            df_cov,
                            nuc_pos_in_duplex_to_count):
    """saves dataframes into an excel spreadsheet
    """
    if not os.path.exists("./output_apobec"):
        os.mkdir("./output_apobec")

    file_name = f.rsplit(".", 1)[0]
    writer = pd.ExcelWriter("./output_apobec/" + file_name + "_" + nuc_pos_in_duplex_to_count + '.xlsx')

    counter = 0
    for df in df_raw_snp_container:
        df.to_excel(writer, "raw_count" + str(counter))
        counter += 1

    for df in df_duplex_container:
        tab_name = df.index[0][1] # changed to get second nuc
        df.to_excel(writer, "nuc_" + tab_name)

    for df in df_perc_container:
        tab_name = df.index[0][1] # changed to get second nuc
        df.to_excel(writer, "nuc_" + tab_name + "_percentage")

    df_cov.to_excel(writer, "coverage")

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


def convert_pivot_df_into_percent(df_duplex_container):
    total_snpes = 0
    for df in df_duplex_container:
        total_snpes += df.sum().sum()
    df_perc_container = []
    for df in df_duplex_container:
        df_perc_container.append(df / total_snpes * 100)
    
    return df_perc_container


#TODO put messages into separate file, into a dict, for example
# and import it to just call the items to show them
# so as to make main function not so large

# TODO delete all the try inside the main
# or just test them prior to it

def main():
    #print(sys.argv[1])

    start_time = time()
    file_counter = 0
    nucleotides = ["A", "T", "G", "C"]

    if sys.argv[1].strip() == "count_first_in_duplex" or \
        sys.argv[1].strip() == "count_second_in_duplex":
        nuc_pos_in_duplex_to_count = sys.argv[1].strip()
        #print(nuc_pos_in_duplex_to_count)
    #############################
    # logging
    if not os.path.exists("./log"):
        os.mkdir("./log")
    logger.remove() # don't put messages into notebook output
    logger.add("./log/apobec_count_snp_duplex_log_{time}" + "_" + nuc_pos_in_duplex_to_count + ".txt", backtrace=False)
    logger.add(sys.stderr, level="CRITICAL")
    #############################################
    
    if os.path.exists("./input_data"): #TODO make if not and then return
    # to get rid of else block below
        input_files = get_input_files_names("./input_data")

        num_files = len(input_files)
        print("""
               ---------------
               job started at {0} ...
               ---------------
               """.format(get_current_time()))

        progress_bar = IntProgress(min=0, max=num_files, bar_style='success')
        display(progress_bar)
        
        for f in input_files:
            try:
                df_raw_snp_container = []
                df_duplex_container = []
                ref_seq = get_ref(f)
                ####################################################
                # nuc we use in loop is very important
                # 1. it's the nuc we use to find all duplexes in reference 
                # which start from it
                # 2. then we search in count_snp func to which nucs it changes 
                # in the same position in the read and in the same context as in the read
                # (i.e. the second nuc in duplex is the same as in reference but the fitst nuc (
                # which is "nuc" var in the for loop below) changes)
                for nuc in nucleotides:
                    snp_type = nucleotides[:]
                    # remove nuc from snp type
                    snp_type.remove(nuc)
                     
                    if nuc_pos_in_duplex_to_count == "count_first_in_duplex":
                        duplex_posits = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, nuc)
                        df = create_df(duplex_posits, snp_type)        
                        df_snp, coverage, record_id = count_snp_with_ref_nuc_at_first_pos(f, df, nuc)
                        #print(df_snp)
                    elif nuc_pos_in_duplex_to_count == "count_second_in_duplex":
                        duplex_posits = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, nuc)
                        df = create_df(duplex_posits, snp_type)        
                        df_snp, coverage, record_id = count_snp_with_ref_nuc_at_second_pos(f, df, nuc)
                        #print(df_snp)
                    
                    df_duplex_in_context = create_pivot_df(df_snp)
                    #print(df_duplex_in_context)
                    #print("index [0][0]", df_duplex_in_context.index[0][0])
                    #print("index [0][1]", df_duplex_in_context.index[0][1])
                    df_raw_snp_container.append(df_snp)
                    df_duplex_container.append(df_duplex_in_context)
                    
                df_perc_container = convert_pivot_df_into_percent(df_duplex_container)
                df_cov = pd.DataFrame.from_dict(
                    {"coverage": coverage}, orient='index')
                
                # two function with condition, so as not to have conditions inside a single function
                if nuc_pos_in_duplex_to_count == "count_first_in_duplex":
                    save_df_first_in_duplex(f,
                                            df_raw_snp_container,
                                            df_duplex_container,
                                            df_perc_container,
                                            df_cov,
                                            nuc_pos_in_duplex_to_count)
                elif nuc_pos_in_duplex_to_count == "count_second_in_duplex":
                    save_df_second_in_duplex(f,
                                            df_raw_snp_container,
                                            df_duplex_container,
                                            df_perc_container,
                                            df_cov,
                                            nuc_pos_in_duplex_to_count)
            
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
