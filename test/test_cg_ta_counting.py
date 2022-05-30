import sys
sys.path.append("..")


import os

from Bio.SeqIO.FastaIO import SimpleFastaParser  # low level fast fasta parser

from apobec.count_cg_to_ta import (get_ref,
                                    get_ref_name,
                                    get_input_files_names,
                                    get_cg_posits_in_ref,
                                    get_num_of_reads_in_file,
                                    count_num_of_c_and_g_in_ref,
                                    count_total_c_and_g_in_reads,
                                    count_c_t_change,
                                    count_g_a_change,
                                    calculate_percent_c_t_change,
                                    calculate_percent_g_a_change,
                                    calculate_c_t_and_g_a_change_percent)






def test_get_num_of_reads_in_file():
    """
    """
    in_file = "data_count_reads_one_read.fasta"
    num_reads_in_file = get_num_of_reads_in_file(in_file)
    assert num_reads_in_file == 0

    in_file = "data_count_reads_ten_reads.fasta"
    num_reads_in_file = get_num_of_reads_in_file(in_file)
    assert num_reads_in_file == 9



def test_count_num_of_c_and_g_in_ref():
    """
    """
    in_file = "data_c_pos.fasta"
    ref_seq = get_ref(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)
    num_c_in_ref, num_g_in_ref = count_num_of_c_and_g_in_ref(cg_posits)
    assert num_c_in_ref == 4
    assert num_g_in_ref == 0

    in_file = "data_g_pos.fasta"
    ref_seq = get_ref(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)
    num_c_in_ref, num_g_in_ref = count_num_of_c_and_g_in_ref(cg_posits)
    assert num_c_in_ref == 0
    assert num_g_in_ref == 4

    in_file = "data_cg_pos.fasta"
    ref_seq = get_ref(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)
    num_c_in_ref, num_g_in_ref = count_num_of_c_and_g_in_ref(cg_posits)
    assert num_c_in_ref == 4
    assert num_g_in_ref == 4



def test_count_total_c_and_g_in_reads():
    """
    """
    in_file = "data_c_pos.fasta"
    ref_seq = get_ref(in_file)
    ref_name = get_ref_name(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)
    num_reads_in_file = get_num_of_reads_in_file(in_file)
    num_c_in_ref, num_g_in_ref = count_num_of_c_and_g_in_ref(cg_posits)
    total_num_c_in_reads, total_num_g_in_reads = count_total_c_and_g_in_reads(ref_name, cg_posits, in_file)
    assert total_num_c_in_reads == 12
    assert total_num_g_in_reads == 0

    in_file = "data_g_pos.fasta"
    ref_seq = get_ref(in_file)
    ref_name = get_ref_name(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)
    num_reads_in_file = get_num_of_reads_in_file(in_file)
    num_c_in_ref, num_g_in_ref = count_num_of_c_and_g_in_ref(cg_posits)
    total_num_c_in_reads, total_num_g_in_reads = count_total_c_and_g_in_reads(ref_name, cg_posits, in_file)
    assert total_num_c_in_reads == 0
    assert total_num_g_in_reads == 12



    in_file = "data_cg_pos.fasta"
    ref_seq = get_ref(in_file)
    ref_name = get_ref_name(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)
    num_reads_in_file = get_num_of_reads_in_file(in_file)
    num_c_in_ref, num_g_in_ref = count_num_of_c_and_g_in_ref(cg_posits)
    total_num_c_in_reads, total_num_g_in_reads = count_total_c_and_g_in_reads(ref_name, cg_posits, in_file)
    assert total_num_c_in_reads == 12
    assert total_num_g_in_reads == 12

    
def test_count_g_a_change(): 

    """
    """
    in_file = "data_c_t_change.fasta"
    ref_seq = get_ref(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)
    num_c_t_change_in_reads = count_c_t_change(in_file, cg_posits)
    num_g_a_change_in_reads = count_g_a_change(in_file, cg_posits)
    assert num_c_t_change_in_reads == 4
    assert num_g_a_change_in_reads == 0


    in_file = "data_g_a_change.fasta"
    ref_seq = get_ref(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)
    num_c_t_change_in_reads = count_c_t_change(in_file, cg_posits)
    num_g_a_change_in_reads = count_g_a_change(in_file, cg_posits)
    assert num_c_t_change_in_reads == 0
    assert num_g_a_change_in_reads == 4


    in_file = "data_cg_ta_change.fasta"
    ref_seq = get_ref(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)
    num_c_t_change_in_reads = count_c_t_change(in_file, cg_posits)
    num_g_a_change_in_reads = count_g_a_change(in_file, cg_posits)
    assert num_c_t_change_in_reads == 4
    assert num_g_a_change_in_reads == 4



def test_calculate_percent_c_t_change():
    """
    """
    in_file = "data_c_t_10_percent_change.fasta" #
    ref_seq = get_ref(in_file)
    ref_name = get_ref_name(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)

    num_reads_in_file = get_num_of_reads_in_file(in_file)
    num_c_in_ref, num_g_in_ref = count_num_of_c_and_g_in_ref(cg_posits)
    total_num_c_in_reads, total_num_g_in_reads = count_total_c_and_g_in_reads(ref_name, cg_posits, in_file)

    num_c_t_change_in_reads = count_c_t_change(in_file, cg_posits)
    num_g_a_change_in_reads = count_g_a_change(in_file, cg_posits)

    percent_c_t_change = calculate_percent_c_t_change(num_c_t_change_in_reads, total_num_c_in_reads,num_c_in_ref)
    percent_g_a_change = calculate_percent_g_a_change(num_g_a_change_in_reads, total_num_g_in_reads, num_g_in_ref)

    assert percent_c_t_change == 10
    assert percent_g_a_change == "no G in reference"



def test_calculate_percent_g_a_change():
    """
    """
    in_file = "data_g_a_10_percent_change.fasta" #
    ref_seq = get_ref(in_file)
    ref_name = get_ref_name(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)

    num_reads_in_file = get_num_of_reads_in_file(in_file)
    num_c_in_ref, num_g_in_ref = count_num_of_c_and_g_in_ref(cg_posits)
    total_num_c_in_reads, total_num_g_in_reads = count_total_c_and_g_in_reads(ref_name, cg_posits, in_file)
    num_c_t_change_in_reads = count_c_t_change(in_file, cg_posits)
    num_g_a_change_in_reads = count_g_a_change(in_file, cg_posits)

    percent_c_t_change = calculate_percent_c_t_change(num_c_t_change_in_reads, total_num_c_in_reads,num_c_in_ref)
    percent_g_a_change = calculate_percent_g_a_change(num_g_a_change_in_reads, total_num_g_in_reads, num_g_in_ref)

    assert percent_c_t_change == "no C in reference"
    assert percent_g_a_change == 10

def test_calculate_percent_c_t_g_a_change():
    """
    """
    in_file = "data_c_t_11_g_a_17_percent_change.fasta" #
    ref_seq = get_ref(in_file)
    ref_name = get_ref_name(in_file)
    cg_posits = get_cg_posits_in_ref(ref_seq)

    num_reads_in_file = get_num_of_reads_in_file(in_file)
    num_c_in_ref, num_g_in_ref = count_num_of_c_and_g_in_ref(cg_posits)
    total_num_c_in_reads, total_num_g_in_reads = count_total_c_and_g_in_reads(ref_name, cg_posits, in_file)
    num_c_t_change_in_reads = count_c_t_change(in_file, cg_posits)
    num_g_a_change_in_reads = count_g_a_change(in_file, cg_posits)

    percent_c_t_change = calculate_percent_c_t_change(num_c_t_change_in_reads, total_num_c_in_reads, num_c_in_ref)
    percent_g_a_change = calculate_percent_g_a_change(num_g_a_change_in_reads, total_num_g_in_reads, num_g_in_ref)

    assert percent_c_t_change == 11
    assert percent_g_a_change == 17



def test_calculate_total_percent_c_t_g_a_change_10_15_15():
    """
    """
    test_percent_files = ["data_ct_10_ga_10.fasta",
                            "data_ct_20_ga_10.fasta",
                            "data_ct_10_ga_20.fasta",
                            "data_ct_1_ga_10.fasta"]
    results = [10, 15, 15, 5.5]
    counter = 0
    for file in test_percent_files:
        in_file = file #
        ref_seq = get_ref(in_file)
        ref_name = get_ref_name(in_file)
        cg_posits = get_cg_posits_in_ref(ref_seq)

        num_reads_in_file = get_num_of_reads_in_file(in_file)
        num_c_in_ref, num_g_in_ref = count_num_of_c_and_g_in_ref(cg_posits)
        total_num_c_in_reads, total_num_g_in_reads = count_total_c_and_g_in_reads(ref_name, cg_posits, in_file)
        num_c_t_change_in_reads = count_c_t_change(in_file, cg_posits)
        num_g_a_change_in_reads = count_g_a_change(in_file, cg_posits)

        percent_c_t_change = calculate_percent_c_t_change(num_c_t_change_in_reads, total_num_c_in_reads, num_c_in_ref)
        percent_g_a_change = calculate_percent_g_a_change(num_g_a_change_in_reads, total_num_g_in_reads, num_g_in_ref)
        
        c_t_and_g_a_change_percent = calculate_c_t_and_g_a_change_percent(total_num_c_in_reads,
                                                                            total_num_g_in_reads,
                                                                        num_c_t_change_in_reads,
                                                                        num_g_a_change_in_reads)
        assert c_t_and_g_a_change_percent == results[counter]
        counter += 1
                                                
    


    

    



               