import pandas as pd

aa = pd.read_csv("../data/amino_acids_with_weight.csv")

combination_list = []
aa_names_list = []
for mass1, aa1, i in zip(aa["mono mass"],aa["letter"], range(len(aa))):
    combination_list.append([])
    aa_names_list.append(aa1)
    for mass2, aa2 in zip(aa["mono mass"],aa["letter"]):
        # combination_list[i].append([mass1 + mass2,(aa1 + "+" + aa2)])
        combination_list[i].append(mass1 + mass2)
# merge combo and aa into a single dataframe maybe this format:
"""
    full    letter  short   comp    mono    mass    G           A           S           P   V
0   glycine G       gly     c2h3no  57      57      g+g mass    g+a mass    g+s mass
"""
# index = aa_names_list
combo = pd.DataFrame(combination_list, columns = aa_names_list)
aa_combo = pd.concat([aa, combo], axis=1)
aa_combo.to_csv("../data/single_double_amino_acids.csv", index=False)
