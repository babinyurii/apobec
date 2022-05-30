import sys
sys.path.append("..")

from apobec.count_cg_to_ta import get_ref, get_cg_posits_in_ref
 
"""

"""    
def test_get_c_positions():
    """
    """

    ref_seq = get_ref("data_c_pos.fasta")



    d = get_cg_posits_in_ref(ref_seq)

    assert d["C_pos"][0] == 0
    assert d["C_pos"][1] == 3
    assert d["C_pos"][2] == 6
    assert d["C_pos"][3] == 9
    assert d["G_pos"] == []


def test_get_g_positions():
    """
    """

    ref_seq = get_ref("data_g_pos.fasta")

    d = get_cg_posits_in_ref(ref_seq)

    assert d["G_pos"][0] == 0
    assert d["G_pos"][1] == 3
    assert d["G_pos"][2] == 6
    assert d["G_pos"][3] == 9
    assert d["C_pos"] == []

def test_get_cg_positions():
    """
    """

    ref_seq = get_ref("data_cg_pos.fasta")

    d = get_cg_posits_in_ref(ref_seq)

    assert d["C_pos"][0] == 0
    assert d["C_pos"][1] == 2
    assert d["C_pos"][2] == 4
    assert d["C_pos"][3] == 6
    assert d["G_pos"][0] == 1
    assert d["G_pos"][1] == 3
    assert d["G_pos"][2] == 5
    assert d["G_pos"][3] == 7
    
    




    

   
    