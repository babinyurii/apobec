import sys
sys.path.append("..")

from apobec.count_snp_duplex import get_duplex_posits_with_ref_nuc_at_first_pos, get_duplex_posits_with_ref_nuc_at_second_pos
 
"""

"""    
def test_get_duplex_posits_T():
    """
    FIRST and second position test
    """

    ref_seq = "ATAAATAATA"

    d = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, "T")
    assert d["start_pos"][0] == 1
    assert d["stop_pos"][0] == 2
    assert d["ref_duplex"][0] == "TA"
    assert d["start_pos"][1] == 5
    assert d["stop_pos"][1] == 6
    assert d["ref_duplex"][1] == "TA"
    assert d["start_pos"][2] == 8
    assert d["stop_pos"][2] == 9
    assert d["ref_duplex"][2] == "TA"

    d = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, "T")
    assert d["start_pos"][0] == 0
    assert d["stop_pos"][0] == 1
    assert d["ref_duplex"][0] == "AT"
    assert d["start_pos"][1] == 4
    assert d["stop_pos"][1] == 5
    assert d["ref_duplex"][1] == "AT"
    assert d["start_pos"][2] == 7
    assert d["stop_pos"][2] == 8
    assert d["ref_duplex"][2] == "AT"

def test_get_duplex_posits_A():
    """
    FIRST and second position test
    """
    ref_seq = "TATTTATTAT"

    d = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, "A")
    assert d["start_pos"][0] == 1
    assert d["stop_pos"][0] == 2
    assert d["ref_duplex"][0] == "AT"
    assert d["start_pos"][1] == 5
    assert d["stop_pos"][1] == 6
    assert d["ref_duplex"][1] == "AT"
    assert d["start_pos"][2] == 8
    assert d["stop_pos"][2] == 9
    assert d["ref_duplex"][2] == "AT"

    d = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, "A")
    assert d["start_pos"][0] == 0
    assert d["stop_pos"][0] == 1
    assert d["ref_duplex"][0] == "TA"
    assert d["start_pos"][1] == 4
    assert d["stop_pos"][1] == 5
    assert d["ref_duplex"][1] == "TA"
    assert d["start_pos"][2] == 7
    assert d["stop_pos"][2] == 8
    assert d["ref_duplex"][2] == "TA"

def test_get_duplex_posits_G():
    """
    FIRST and second position test
    """
    ref_seq = "TGTTTGTTGT"

    d = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, "G")
    assert d["start_pos"][0] == 1
    assert d["stop_pos"][0] == 2
    assert d["ref_duplex"][0] == "GT"
    assert d["start_pos"][1] == 5
    assert d["stop_pos"][1] == 6
    assert d["ref_duplex"][1] == "GT"
    assert d["start_pos"][2] == 8
    assert d["stop_pos"][2] == 9
    assert d["ref_duplex"][2] == "GT"

    d = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, "G")
    assert d["start_pos"][0] == 0
    assert d["stop_pos"][0] == 1
    assert d["ref_duplex"][0] == "TG"
    assert d["start_pos"][1] == 4
    assert d["stop_pos"][1] == 5
    assert d["ref_duplex"][1] == "TG"
    assert d["start_pos"][2] == 7
    assert d["stop_pos"][2] == 8
    assert d["ref_duplex"][2] == "TG"

def test_get_duplex_posits_C():
    """
    FIRST and second position test
    """
    ref_seq = "TCTTTCTTCT"

    d = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, "C")
    assert d["start_pos"][0] == 1
    assert d["stop_pos"][0] == 2
    assert d["ref_duplex"][0] == "CT"
    assert d["start_pos"][1] == 5
    assert d["stop_pos"][1] == 6
    assert d["ref_duplex"][1] == "CT"
    assert d["start_pos"][2] == 8
    assert d["stop_pos"][2] == 9
    assert d["ref_duplex"][2] == "CT"

    d = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, "C")
    assert d["start_pos"][0] == 0
    assert d["stop_pos"][0] == 1
    assert d["ref_duplex"][0] == "TC"
    assert d["start_pos"][1] == 4
    assert d["stop_pos"][1] == 5
    assert d["ref_duplex"][1] == "TC"
    assert d["start_pos"][2] == 7
    assert d["stop_pos"][2] == 8
    assert d["ref_duplex"][2] == "TC"

def test_get_duplex_posits_at_first_pos_at_the_start_and_at_the_end():
    """
    """
    ref_seq = "CTTTTTTTCT"

    d = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, "C")
    assert d["start_pos"][0] == 0
    assert d["stop_pos"][0] == 1
    assert d["ref_duplex"][0] == "CT"
    assert d["start_pos"][1] == 8
    assert d["stop_pos"][1] == 9
    assert d["ref_duplex"][1] == "CT"

    # nucs to find is at the end without duplexes: like false positives
    ref_seq = "TTTTTTTTTC"

    d = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, "C")
    assert d["start_pos"] == []
    assert d["stop_pos"] == []
    assert d["ref_duplex"] == []
    

def test_get_duplex_posits_at_second_pos_at_the_start_and_at_the_end():
    """
    """
    ref_seq = "TCTTTTTTTC"

    d = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, "C")
    assert d["start_pos"][0] == 0
    assert d["stop_pos"][0] == 1
    assert d["ref_duplex"][0] == "TC"
    assert d["start_pos"][1] == 8
    assert d["stop_pos"][1] == 9
    assert d["ref_duplex"][1] == "TC"

    # nucs to find is at the start without duplexes: like false positives
    ref_seq = "CTTTTTTTT"

    d = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, "C")
    assert d["start_pos"] == []
    assert d["stop_pos"] == []
    assert d["ref_duplex"] == []
    
  

def test_get_duplex_posits_at_first_pos_with_gaps():
    """
    """
    ref_seq = "C-TTTC---TTTC--T"

    d = get_duplex_posits_with_ref_nuc_at_first_pos(ref_seq, "C")
    assert d["start_pos"][0] == 0
    assert d["stop_pos"][0] == 2
    assert d["ref_duplex"][0] == "CT"
    assert d["start_pos"][1] == 5
    assert d["stop_pos"][1] == 9
    assert d["ref_duplex"][1] == "CT"
    assert d["start_pos"][2] == 12
    assert d["stop_pos"][2] == 15
    assert d["ref_duplex"][2] == "CT"

def test_get_duplex_posits_at_second_pos_with_gaps():
    """
    """
    ref_seq = "T-CTTT---CTTT--C"

    d = get_duplex_posits_with_ref_nuc_at_second_pos(ref_seq, "C")
    assert d["start_pos"][0] == 0
    assert d["stop_pos"][0] == 2
    assert d["ref_duplex"][0] == "TC"
    assert d["start_pos"][1] == 5
    assert d["stop_pos"][1] == 9
    assert d["ref_duplex"][1] == "TC"
    assert d["start_pos"][2] == 12
    assert d["stop_pos"][2] == 15
    assert d["ref_duplex"][2] == "TC"

   
    