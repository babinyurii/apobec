# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 10:32:44 2019
@author: babin
"""

""" snp_rate.py
this script count exact number of SNPs in each read from the fasta alignment
It writes results into the excel spreadsheet and constructs distribution plot.
the script compares each sequence with the reference
and counts number of SNPs.
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
    3. run in jupyter : %run ./scripts/snp_rate.py

NOTE: to speed up the exection instead of 'SeqIO.parse()
low level fasta parser 'SimpleFastaParser' is used 
"""

import datetime
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio.SeqIO.FastaIO import SimpleFastaParser  # low level fast fasta parser
from time import time
from IPython.display import display
from ipywidgets import IntProgress
sns.set()


def get_input_files_names(path_to_input):
    """returns list of fasta files
    """
    fasta_extensions = ["fasta", "fa", "fas"]
    input_files = os.listdir(path_to_input)
    input_files = [f for f in input_files if f.rsplit(
        ".", 1)[-1] in fasta_extensions]

    return input_files


def get_ref(in_file):
    """returns reference sequence
    which is the first record in a fasta file
    """
    with open("./input_data/" + in_file) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            ref_seq = seq
            break

    return ref_seq


def count_snp_per_read(ref_seq, in_file):
    """returns a list of integers:
    - len of the list equals coverage, i.e. num of reads in file
    - each int value represents a num of snp in a read
    """
    snpes_per_read = []

    with open("./input_data/" + in_file) as in_handle:
        reads = SimpleFastaParser(in_handle)
        # coverage_counter = 0

        for record in reads:
            read = record[1]
            snp_counter = 0
            for ind in range(0, len(ref_seq)):
                if read[ind] != ref_seq[ind] and read[ind] != "-":
                    if ref_seq[ind] != "-":  # skip gaps, insertions are not counted
                        snp_counter += 1

            #coverage_counter += 1
            snpes_per_read.append(snp_counter)

        series_name = in_file.rsplit(".", 1)[0]
        snpes_per_read = pd.Series(snpes_per_read, name=series_name)

    return snpes_per_read


def get_mode_median(df_raw_count):
    """calcualtes mode and median
    using raw snp count
    """
    mode = df_raw_count.mode()
    # put in list, otherwise an error "a collection of some kind"
    mode.index = ["mode"]
    median = pd.DataFrame(df_raw_count.median()).T
    median.index = ["median"]

    return mode, median


def calculate_non_mutated(pivot_df, df_raw_count):
    """counts zero values in columns of 
    df_raw_count. Zero values represent reads 
    without SNPS
    """
    cols = []
    non_mutated = []
    for col in df_raw_count:
        col_in_list = df_raw_count[col].tolist()
        cols.append(col)
        non_mutated.append(col_in_list.count(0))

    # list is in [] as nested, to get values in a row in df
    df_non_mutated = pd.DataFrame([non_mutated], columns=cols)
    df_non_mutated.rename(
        index={0: "number of non-mutated reads"}, inplace=True)
    pivot_df = pd.concat([pivot_df, df_non_mutated])

    # calculating percent of non-mutated reads
    pivot_df.loc['percent non-mutated reads'] = pivot_df.loc["number of non-mutated reads"] / \
        pivot_df.loc["reads total"] * 100

    return pivot_df


def create_strip_plot(df, time_stamp):
    """constructs distribution plot
    using seaborn 'stripplot'
    """
    fig = plt.figure(figsize=(15, 10))
    sns.stripplot(data=df.loc[:, :], jitter=True, size=4)
    plt.xticks(rotation=90)

    time_stamp = time_stamp.split(":")
    time_stamp = "-".join(time_stamp)
    plt.xticks(fontsize=10)
    plt.ylabel("number of SNPs in a read", fontsize=15)
    plt.xlabel("input file name", fontsize=15)
    plt.title("SNP distribution", fontsize=20)
    plt.tight_layout()
    fig.savefig("./output_apobec/mutation_rate_distribution_" +
                time_stamp + "_.png")
    # plt.show() # comment to save figure, otherwise it'll save blank file
    plt.close()  # comment the line to show the figure in the jupyter or wherever


def save_df(pivot_df, df_raw_count, time_stamp):
    """saves dataframes into an excel spreadsheet
    """
    if not os.path.exists("./output_apobec"):
        os.mkdir("./output_apobec")

    time_stamp = time_stamp.split(":")
    time_stamp = "-".join(time_stamp)

    writer = pd.ExcelWriter(
        "./output_apobec/mutation_rate_" + time_stamp + '_.xlsx')
    pivot_df.to_excel(writer, "pivot spreadsheet")
    df_raw_count.to_excel(writer, "snp raw count")

    writer.save()


def get_current_time():
    """ returns time stamp
    """
    time_stamp = datetime.datetime.fromtimestamp(
        time()).strftime('%Y-%m-%d %H:%M:%S')

    return time_stamp


def show_report(total_time, total_files):
    """prints out brief report 
    """

    hours = total_time // 3600
    minutes = (total_time % 3600) // 60
    seconds = total_time % 60

    print("""
    
    file processed: {0}
    time taken: {1} hours {2} minutes {3} seconds
    
    the results are in the folder 'output_apobec'
    
    note: difference in the output files names 
    is denoted by the time when script finished its work.
    
    ---------------
    ... job finished at {4}
    ---------------
    
    """.format(total_files, hours, minutes, int(seconds), get_current_time()))


def main():

    start_time = time()
    file_counter = 0

    if os.path.exists("./input_data"):
        input_files = get_input_files_names("./input_data")

        num_files = len(input_files)

        print("""
               ---------------
               job started at {0} ...
               ---------------
               """.format(get_current_time()))

        print("collecting snp...")
        progress_bar = IntProgress(min=0, max=num_files, bar_style='success')
        display(progress_bar)

        series_container = []
        for f in input_files:
            try:
                ref_seq = get_ref(f)
                snpes_per_read = count_snp_per_read(ref_seq, f)
            except PermissionError as permerr:
                print("""
                      warning: It seems that the file'{0}' is open in some other 
                      application, like 'Geneious' or whatever,
                      which doesn't allow it to be processed
                      Now this file is skipped
                       """.format(f))

            series_container.append(snpes_per_read)
            progress_bar.value += 1
            file_counter += 1

        print("snp collected...")
        df_raw_count = pd.DataFrame(series_container).T
        # drop first raw which is a reference to itself comparison
        df_raw_count.drop([0], inplace=True)
        pivot = df_raw_count.describe()
        mode, median = get_mode_median(df_raw_count)
        pivot_df = pd.concat([pivot, mode, median])

        # renaming index values in pivot df
        pivot_df.rename(index={'count': 'reads total',
                               'mean': 'mean snp per read',
                               'min': 'min number of snp',
                               'max': 'max number of snp',
                               'mode': 'mode snp',
                               'median': 'median snp'
                               }, inplace=True)

        pivot_df = calculate_non_mutated(pivot_df, df_raw_count).T

        time_stamp = get_current_time()
        
        try:
            create_strip_plot(df_raw_count, time_stamp)
            print("plot consctructed...")
        except:
            print("something went wrong, plot construction is skipped...")

        df_raw_count = df_raw_count.T  # transpose after plotting
        save_df(pivot_df, df_raw_count, time_stamp)
        print("spreadsheet saved...")

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
