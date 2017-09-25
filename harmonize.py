# Input files:
# glgc (http://csg.sph.umich.edu/abecasis/public/lipids2013/​)
# test_data
# snp_info (for test_data)

# Output file:
# corrected_data

import pandas as pd

g = pd.read_table("glgc", sep="\s+", index_col="rsid", header=0)
g.loc[:, "A2"] = g.A2.str.upper()

td = pd.read_table("test_data.txt", sep=" ")

si = pd.read_table("snp_info.txt", sep=" ", header=0)
si = si.set_index("rs", append=False)

r = pd.read_table("corrected_data.txt", sep=" ")

j = g.join(si)["A1 A2 REF ALT AF MAF FreqA1".split()]

switch_alleles = {"G": "C", "C": "G", "A": "T", "T": "A"}

j.insert(1, "A1Reverse", j.A1.replace(switch_alleles))

# different: rs11206510_dose, rs2479394_dose, rs2479409_dose

correct_strand = (j.A1 == j.ALT) | (j.A1Reverse == j.ALT)
palindromic = (j.A1 == j.A2.replace(switch_alleles))
# palindromic_and_inferrable = j[palindromic & ]
