# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:02:52 2019

@author: babin
"""

import datetime
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from ipywidgets import IntProgress
from IPython.display import display
from time import time
sns.set()


def get_input_files_names(path_to_input):
    """returns a list of input  files
    """
    input_files = os.listdir(path_to_input)
    input_files = [f for f in input_files if f.rsplit(
        ".", 1)[-1] == "xlsx"]
    
    return input_files

def get_current_time():
    """just returns time stamp
    """
    time_stamp = datetime.datetime.fromtimestamp(
        time()).strftime('%Y-%m-%d %H:%M:%S')

    return time_stamp



def get_largest_percent(df_perc_container):
    largest_percent_in_df = []
    for df in df_perc_container:
        largest_percent_in_df.append(df.max().max())
    largest_percent = max(largest_percent_in_df)
    
    return largest_percent

def create_bar_chart(file_name, df_perc_container, largest_percent):
    """
    creates bar chart using list of dataframes with 
    snp percentage. Percentage is calculated relatively 
    total number of snp found
    """
    for df in df_perc_container:
        n_samples = len(df.index)
        pos = [x for x in range(0, n_samples)] # number of bar stacks, which is n_features in df
        width = 0.25
        fig, ax = plt.subplots(figsize=(10, 5))
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
        fig.savefig("./output_duplex_plots/" + file_name.rsplit(".", 1)[0] +"_" + nuc + "_1st_pos_.png")
        # plt.show() # comment to save figure, otherwise it'll save blank file
        plt.close(fig)  # comment the line to show the figure in the jupyter or wherever


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


def main():
    start_time = time()
    file_counter = 0
    sheet_names = ["nuc_A_percentage", 
                   "nuc_T_percentage", 
                   "nuc_G_percentage", 
                   "nuc_C_percentage" ]
    
    
    
    
    if os.path.exists("./output_apobec"):
        input_files = get_input_files_names("./output_apobec")
    
    if not os.path.exists("./output_duplex_plots"):
        os.mkdir("./output_duplex_plots")
    
    print("""
               ---------------
               job started at {0} ...
               ---------------
               """.format(get_current_time()))
    
    num_files = len(input_files)
    progress_bar = IntProgress(min=0, max=num_files, bar_style='success')
    display(progress_bar)
    
    for f in input_files:
        
        df_perc_container = []
        
        try:
            for sheet in sheet_names:
                df = pd.read_excel("./output_apobec/%s" % f, sheet_name=sheet)
                
                # the following two lines needed to assign new row index
                # and delete 'ref_duplex' column
                # as we want to delete index created by 'read_excel' method
                df.index = df["ref_duplex"]
                df = df.drop(["ref_duplex"], axis=1)
                df_perc_container.append(df)
                
        except PermissionError as permerr:
                print("""
                      warning: It seems that excel spreadsheet '{0}' is open in excel.
                      To process the corresponding file, rerun the script
                      now this file is skipped. error: {1}
                      """.format(f, permerr))
                
                progress_bar.value += 1
                continue
                
            
        
        try:
            largest_percent = get_largest_percent(df_perc_container)
            create_bar_chart(f, df_perc_container, largest_percent)
            
        except PermissionError as permerr:
                print("""
                      warning: It seems that excel spreadsheet '{0}' is open in excel.
                      To process the corresponding file, rerun the script
                      now this file is skipped. error: {1}
                      """.format(f, permerr))
        except ValueError as valerr:
                print(""""
                      plot can't be constructed. 
                      sequences in the input file may be too short. {0}
                      """.format(valerr))
        except IndexError as inderr:
            print("index error:%s" % inderr)
        
        
        progress_bar.value += 1
        file_counter += 1

    finish_time = time()
    total_time = finish_time - start_time
    show_report(total_time, file_counter)
        
 
    
if __name__ == "__main__":
    main()


