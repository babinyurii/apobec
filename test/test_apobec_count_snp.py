import sys
sys.path.append("..")


import os

from Bio.SeqIO.FastaIO import SimpleFastaParser  # low level fast fasta parser

from apobec.count_snp_duplex import (get_input_files_names, 
get_ref, 
get_duplex_posits_with_ref_nuc_at_first_pos,
get_duplex_posits_with_ref_nuc_at_second_pos,
create_df,
count_snp_with_ref_nuc_at_first_pos,
count_snp_with_ref_nuc_at_second_pos,
create_pivot_df)


'''
test_file_1 is :
------------------------
>ref_seq
TTCTTCTTCTCT
>test_seq_1
TTATTGTTCTTT
>test_seq_1
TTCTTGTTCTTT
>test_seq_1
TTCTTCTTCTTT 
'''  


def test_snp_df_nuc_C():
    print(os.listdir("./input_data"))

    input_files = get_input_files_names("./input_data")
    #print("input_files: ", input_files)

    ref_seq = get_ref("data_C.fasta")
    input_file = "data_C.fasta"
    #print(ref_seq)

    nuc = "C"
    snp_type = ["A", "T", "G"]

    duplex_posits_at_first = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, nuc)
    duplex_posits_at_second = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, nuc)

    
    df_first_pos = create_df(duplex_posits_at_first, snp_type)    
    df_second_pos = create_df(duplex_posits_at_second, snp_type)    
    df_snp_first_pos, coverage_first_pos, record_id_first_pos = count_snp_with_ref_nuc_at_first_pos(input_file, df_first_pos, nuc)
    df_snp_second_pos, coverage_second_pos, record_id_second_pos = count_snp_with_ref_nuc_at_second_pos(input_file, df_second_pos, nuc)
    
    # uncomment to have a look at dataframes
    #print(df_snp_first_pos)
    #print(df_snp_second_pos)

    #########################################################
    # at first position
    assert df_snp_first_pos.loc[0, "start_pos"] == 2
    assert df_snp_first_pos.loc[0, "stop_pos"] == 3
    assert df_snp_first_pos.loc[0, "ref_duplex"] == "CT"
    assert df_snp_first_pos.loc[0, "A"] == 1
    assert df_snp_first_pos.loc[0, "T"] == 0
    assert df_snp_first_pos.loc[0, "G"] == 0

    assert df_snp_first_pos.loc[1, "start_pos"] == 5
    assert df_snp_first_pos.loc[1, "stop_pos"] == 6
    assert df_snp_first_pos.loc[1, "ref_duplex"] == "CT"
    assert df_snp_first_pos.loc[1, "A"] == 0
    assert df_snp_first_pos.loc[1, "T"] == 0
    assert df_snp_first_pos.loc[1, "G"] == 2

    assert df_snp_first_pos.loc[2, "start_pos"] == 8
    assert df_snp_first_pos.loc[2, "stop_pos"] == 9
    assert df_snp_first_pos.loc[2, "ref_duplex"] == "CT"
    assert df_snp_first_pos.loc[2, "A"] == 0
    assert df_snp_first_pos.loc[2, "T"] == 3
    assert df_snp_first_pos.loc[2, "G"] == 0

    assert df_snp_first_pos.loc[3, "start_pos"] == 10
    assert df_snp_first_pos.loc[3, "stop_pos"] == 11
    assert df_snp_first_pos.loc[3, "ref_duplex"] == "CT"
    assert df_snp_first_pos.loc[3, "A"] == 0
    assert df_snp_first_pos.loc[3, "T"] == 0
    assert df_snp_first_pos.loc[3, "G"] == 0
    
    ######################################################
    # at second position
    assert df_snp_second_pos.loc[0, "start_pos"] == 1
    assert df_snp_second_pos.loc[0, "stop_pos"] == 2
    assert df_snp_second_pos.loc[0, "ref_duplex"] == "TC"
    assert df_snp_second_pos.loc[0, "A"] == 1
    assert df_snp_second_pos.loc[0, "T"] == 0
    assert df_snp_second_pos.loc[0, "G"] == 0

    assert df_snp_second_pos.loc[1, "start_pos"] == 4
    assert df_snp_second_pos.loc[1, "stop_pos"] == 5
    assert df_snp_second_pos.loc[1, "ref_duplex"] == "TC"
    assert df_snp_second_pos.loc[1, "A"] == 0
    assert df_snp_second_pos.loc[1, "T"] == 0
    assert df_snp_second_pos.loc[1, "G"] == 2

    assert df_snp_second_pos.loc[2, "start_pos"] == 7
    assert df_snp_second_pos.loc[2, "stop_pos"] == 8
    assert df_snp_second_pos.loc[2, "ref_duplex"] == "TC"
    assert df_snp_second_pos.loc[2, "A"] == 0
    assert df_snp_second_pos.loc[2, "T"] == 3
    assert df_snp_second_pos.loc[2, "G"] == 0

    assert df_snp_second_pos.loc[3, "start_pos"] == 9
    assert df_snp_second_pos.loc[3, "stop_pos"] == 10
    assert df_snp_second_pos.loc[3, "ref_duplex"] == "TC"
    assert df_snp_second_pos.loc[3, "A"] == 0
    assert df_snp_second_pos.loc[3, "T"] == 0
    assert df_snp_second_pos.loc[3, "G"] == 0



def test_snp_df_nuc_T():
    #print(os.listdir("./input_data"))

    input_files = get_input_files_names("./input_data")
    #print("input_files: ", input_files)

    ref_seq = get_ref("data_T.fasta")
    input_file = "data_T.fasta"
    #print(ref_seq)

    nuc = "T"
    snp_type = ["C", "A", "G"]

    duplex_posits_at_first = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, nuc)
    duplex_posits_at_second = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, nuc)

    
    df_first_pos = create_df(duplex_posits_at_first, snp_type)    
    df_second_pos = create_df(duplex_posits_at_second, snp_type)    
    df_snp_first_pos, coverage_first_pos, record_id_first_pos = count_snp_with_ref_nuc_at_first_pos(input_file, df_first_pos, nuc)
    df_snp_second_pos, coverage_second_pos, record_id_second_pos = count_snp_with_ref_nuc_at_second_pos(input_file, df_second_pos, nuc)
    
    # uncomment to have a look at dataframes
    #print(df_snp_first_pos)
    #print(df_snp_second_pos)
    #########################################################
    # at first position
    assert df_snp_first_pos.loc[0, "start_pos"] == 2
    assert df_snp_first_pos.loc[0, "stop_pos"] == 3
    assert df_snp_first_pos.loc[0, "ref_duplex"] == "TA"
    assert df_snp_first_pos.loc[0, "A"] == 1
    assert df_snp_first_pos.loc[0, "C"] == 0
    assert df_snp_first_pos.loc[0, "G"] == 0

    assert df_snp_first_pos.loc[1, "start_pos"] == 5
    assert df_snp_first_pos.loc[1, "stop_pos"] == 6
    assert df_snp_first_pos.loc[1, "ref_duplex"] == "TA"
    assert df_snp_first_pos.loc[1, "A"] == 0
    assert df_snp_first_pos.loc[1, "C"] == 0
    assert df_snp_first_pos.loc[1, "G"] == 2

    assert df_snp_first_pos.loc[2, "start_pos"] == 8
    assert df_snp_first_pos.loc[2, "stop_pos"] == 9
    assert df_snp_first_pos.loc[2, "ref_duplex"] == "TA"
    assert df_snp_first_pos.loc[2, "A"] == 0
    assert df_snp_first_pos.loc[2, "C"] == 3
    assert df_snp_first_pos.loc[2, "G"] == 0

    assert df_snp_first_pos.loc[3, "start_pos"] == 10
    assert df_snp_first_pos.loc[3, "stop_pos"] == 11
    assert df_snp_first_pos.loc[3, "ref_duplex"] == "TA"
    assert df_snp_first_pos.loc[3, "A"] == 0
    assert df_snp_first_pos.loc[3, "C"] == 0
    assert df_snp_first_pos.loc[3, "G"] == 0
    
    ######################################################
    # at second position
    assert df_snp_second_pos.loc[0, "start_pos"] == 1
    assert df_snp_second_pos.loc[0, "stop_pos"] == 2
    assert df_snp_second_pos.loc[0, "ref_duplex"] == "AT"
    assert df_snp_second_pos.loc[0, "A"] == 1
    assert df_snp_second_pos.loc[0, "C"] == 0
    assert df_snp_second_pos.loc[0, "G"] == 0

    assert df_snp_second_pos.loc[1, "start_pos"] == 4
    assert df_snp_second_pos.loc[1, "stop_pos"] == 5
    assert df_snp_second_pos.loc[1, "ref_duplex"] == "AT"
    assert df_snp_second_pos.loc[1, "A"] == 0
    assert df_snp_second_pos.loc[1, "C"] == 0
    assert df_snp_second_pos.loc[1, "G"] == 2

    assert df_snp_second_pos.loc[2, "start_pos"] == 7
    assert df_snp_second_pos.loc[2, "stop_pos"] == 8
    assert df_snp_second_pos.loc[2, "ref_duplex"] == "AT"
    assert df_snp_second_pos.loc[2, "A"] == 0
    assert df_snp_second_pos.loc[2, "C"] == 3
    assert df_snp_second_pos.loc[2, "G"] == 0

    assert df_snp_second_pos.loc[3, "start_pos"] == 9
    assert df_snp_second_pos.loc[3, "stop_pos"] == 10
    assert df_snp_second_pos.loc[3, "ref_duplex"] == "AT"
    assert df_snp_second_pos.loc[3, "A"] == 0
    assert df_snp_second_pos.loc[3, "C"] == 0
    assert df_snp_second_pos.loc[3, "G"] == 0

def test_snp_df_nuc_G():
    #print(os.listdir("./input_data"))

    input_files = get_input_files_names("./input_data")
    #print("input_files: ", input_files)

    ref_seq = get_ref("data_G.fasta")
    input_file = "data_G.fasta"
    #print(ref_seq)

    nuc = "G"
    snp_type = ["A", "T", "C"]

    duplex_posits_at_first = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, nuc)
    duplex_posits_at_second = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, nuc)

    
    df_first_pos = create_df(duplex_posits_at_first, snp_type)    
    df_second_pos = create_df(duplex_posits_at_second, snp_type)    
    df_snp_first_pos, coverage_first_pos, record_id_first_pos = count_snp_with_ref_nuc_at_first_pos(input_file, df_first_pos, nuc)
    df_snp_second_pos, coverage_second_pos, record_id_second_pos = count_snp_with_ref_nuc_at_second_pos(input_file, df_second_pos, nuc)
    
    # uncomment to have a look at dataframes
    #print(df_snp_first_pos)
    #print(df_snp_second_pos)
    #########################################################
    # at first position
    assert df_snp_first_pos.loc[0, "start_pos"] == 2
    assert df_snp_first_pos.loc[0, "stop_pos"] == 3
    assert df_snp_first_pos.loc[0, "ref_duplex"] == "GT"
    assert df_snp_first_pos.loc[0, "A"] == 1
    assert df_snp_first_pos.loc[0, "C"] == 0
    assert df_snp_first_pos.loc[0, "T"] == 0

    assert df_snp_first_pos.loc[1, "start_pos"] == 5
    assert df_snp_first_pos.loc[1, "stop_pos"] == 6
    assert df_snp_first_pos.loc[1, "ref_duplex"] == "GT"
    assert df_snp_first_pos.loc[1, "A"] == 0
    assert df_snp_first_pos.loc[1, "C"] == 2
    assert df_snp_first_pos.loc[1, "T"] == 0

    assert df_snp_first_pos.loc[2, "start_pos"] == 8
    assert df_snp_first_pos.loc[2, "stop_pos"] == 9
    assert df_snp_first_pos.loc[2, "ref_duplex"] == "GT"
    assert df_snp_first_pos.loc[2, "A"] == 0
    assert df_snp_first_pos.loc[2, "C"] == 0
    assert df_snp_first_pos.loc[2, "T"] == 3

    assert df_snp_first_pos.loc[3, "start_pos"] == 10
    assert df_snp_first_pos.loc[3, "stop_pos"] == 11
    assert df_snp_first_pos.loc[3, "ref_duplex"] == "GT"
    assert df_snp_first_pos.loc[3, "A"] == 0
    assert df_snp_first_pos.loc[3, "C"] == 0
    assert df_snp_first_pos.loc[3, "T"] == 0
    
    ######################################################
    # at second position
    assert df_snp_second_pos.loc[0, "start_pos"] == 1
    assert df_snp_second_pos.loc[0, "stop_pos"] == 2
    assert df_snp_second_pos.loc[0, "ref_duplex"] == "TG"
    assert df_snp_second_pos.loc[0, "A"] == 1
    assert df_snp_second_pos.loc[0, "C"] == 0
    assert df_snp_second_pos.loc[0, "T"] == 0

    assert df_snp_second_pos.loc[1, "start_pos"] == 4
    assert df_snp_second_pos.loc[1, "stop_pos"] == 5
    assert df_snp_second_pos.loc[1, "ref_duplex"] == "TG"
    assert df_snp_second_pos.loc[1, "A"] == 0
    assert df_snp_second_pos.loc[1, "C"] == 2
    assert df_snp_second_pos.loc[1, "T"] == 0

    assert df_snp_second_pos.loc[2, "start_pos"] == 7
    assert df_snp_second_pos.loc[2, "stop_pos"] == 8
    assert df_snp_second_pos.loc[2, "ref_duplex"] == "TG"
    assert df_snp_second_pos.loc[2, "A"] == 0
    assert df_snp_second_pos.loc[2, "C"] == 0
    assert df_snp_second_pos.loc[2, "T"] == 3

    assert df_snp_second_pos.loc[3, "start_pos"] == 9
    assert df_snp_second_pos.loc[3, "stop_pos"] == 10
    assert df_snp_second_pos.loc[3, "ref_duplex"] == "TG"
    assert df_snp_second_pos.loc[3, "A"] == 0
    assert df_snp_second_pos.loc[3, "C"] == 0
    assert df_snp_second_pos.loc[3, "T"] == 0

def test_snp_df_nuc_A():
    #print(os.listdir("./input_data"))

    input_files = get_input_files_names("./input_data")
    #print("input_files: ", input_files)

    ref_seq = get_ref("data_A.fasta")
    input_file = "data_A.fasta"
    #print(ref_seq)

    nuc = "A"
    snp_type = ["T", "G", "C"]

    duplex_posits_at_first = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, nuc)
    duplex_posits_at_second = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, nuc)

    
    df_first_pos = create_df(duplex_posits_at_first, snp_type)    
    df_second_pos = create_df(duplex_posits_at_second, snp_type)    
    df_snp_first_pos, coverage_first_pos, record_id_first_pos = count_snp_with_ref_nuc_at_first_pos(input_file, df_first_pos, nuc)
    df_snp_second_pos, coverage_second_pos, record_id_second_pos = count_snp_with_ref_nuc_at_second_pos(input_file, df_second_pos, nuc)
    
    # uncomment to have a look at dataframes
    #print(df_snp_first_pos)
    #print(df_snp_second_pos)
    #########################################################
    # at first position
    assert df_snp_first_pos.loc[0, "start_pos"] == 2
    assert df_snp_first_pos.loc[0, "stop_pos"] == 3
    assert df_snp_first_pos.loc[0, "ref_duplex"] == "AT"
    assert df_snp_first_pos.loc[0, "G"] == 0
    assert df_snp_first_pos.loc[0, "C"] == 0
    assert df_snp_first_pos.loc[0, "T"] == 1

    assert df_snp_first_pos.loc[1, "start_pos"] == 5
    assert df_snp_first_pos.loc[1, "stop_pos"] == 6
    assert df_snp_first_pos.loc[1, "ref_duplex"] == "AT"
    assert df_snp_first_pos.loc[1, "G"] == 0
    assert df_snp_first_pos.loc[1, "C"] == 2
    assert df_snp_first_pos.loc[1, "T"] == 0

    assert df_snp_first_pos.loc[2, "start_pos"] == 8
    assert df_snp_first_pos.loc[2, "stop_pos"] == 9
    assert df_snp_first_pos.loc[2, "ref_duplex"] == "AT"
    assert df_snp_first_pos.loc[2, "G"] == 3
    assert df_snp_first_pos.loc[2, "C"] == 0
    assert df_snp_first_pos.loc[2, "T"] == 0

    assert df_snp_first_pos.loc[3, "start_pos"] == 10
    assert df_snp_first_pos.loc[3, "stop_pos"] == 11
    assert df_snp_first_pos.loc[3, "ref_duplex"] == "AT"
    assert df_snp_first_pos.loc[3, "G"] == 0
    assert df_snp_first_pos.loc[3, "C"] == 0
    assert df_snp_first_pos.loc[3, "T"] == 0
    
    ######################################################
    # at second position
    assert df_snp_second_pos.loc[0, "start_pos"] == 1
    assert df_snp_second_pos.loc[0, "stop_pos"] == 2
    assert df_snp_second_pos.loc[0, "ref_duplex"] == "TA"
    assert df_snp_second_pos.loc[0, "G"] == 0
    assert df_snp_second_pos.loc[0, "C"] == 0
    assert df_snp_second_pos.loc[0, "T"] == 1

    assert df_snp_second_pos.loc[1, "start_pos"] == 4
    assert df_snp_second_pos.loc[1, "stop_pos"] == 5
    assert df_snp_second_pos.loc[1, "ref_duplex"] == "TA"
    assert df_snp_second_pos.loc[1, "G"] == 0
    assert df_snp_second_pos.loc[1, "C"] == 2
    assert df_snp_second_pos.loc[1, "T"] == 0

    assert df_snp_second_pos.loc[2, "start_pos"] == 7
    assert df_snp_second_pos.loc[2, "stop_pos"] == 8
    assert df_snp_second_pos.loc[2, "ref_duplex"] == "TA"
    assert df_snp_second_pos.loc[2, "G"] == 3
    assert df_snp_second_pos.loc[2, "C"] == 0
    assert df_snp_second_pos.loc[2, "T"] == 0

    assert df_snp_second_pos.loc[3, "start_pos"] == 9
    assert df_snp_second_pos.loc[3, "stop_pos"] == 10
    assert df_snp_second_pos.loc[3, "ref_duplex"] == "TA"
    assert df_snp_second_pos.loc[3, "G"] == 0
    assert df_snp_second_pos.loc[3, "C"] == 0
    assert df_snp_second_pos.loc[3, "T"] == 0

def test_create_pivot_df():
    print(os.listdir("./input_data"))

    input_files = get_input_files_names("./input_data")
    #print("input_files: ", input_files)

    ref_seq = get_ref("data_C_1.fasta")
    input_file = "data_C_1.fasta"
    #print(ref_seq)

    nuc = "C"
    snp_type = ["A", "T", "G"]

    duplex_posits_at_first = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, nuc)
    duplex_posits_at_second = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, nuc)

    
    df_first_pos = create_df(duplex_posits_at_first, snp_type)    
    df_second_pos = create_df(duplex_posits_at_second, snp_type)    
    df_snp_first_pos, coverage_first_pos, record_id_first_pos = count_snp_with_ref_nuc_at_first_pos(input_file, df_first_pos, nuc)
    df_snp_second_pos, coverage_second_pos, record_id_second_pos = count_snp_with_ref_nuc_at_second_pos(input_file, df_second_pos, nuc)

    pivot_df_C_first_pos = create_pivot_df(df_snp_first_pos)
    pivot_df_C_second_pos = create_pivot_df(df_snp_second_pos)

    print(pivot_df_C_first_pos)
    print(pivot_df_C_second_pos)

    assert pivot_df_C_first_pos.loc["CT", "A"] == 4
    assert pivot_df_C_first_pos.loc["CT", "G"] == 2
    assert pivot_df_C_first_pos.loc["CT", "T"] == 3

    assert pivot_df_C_second_pos.loc["TC", "A"] == 4
    assert pivot_df_C_second_pos.loc["TC", "G"] == 2
    assert pivot_df_C_second_pos.loc["TC", "T"] == 3
  

    


  
