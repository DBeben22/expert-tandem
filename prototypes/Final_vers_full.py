# basic dependencies and useful math/organization
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import random
from scipy.stats import chi2

# to read mzML files
from pyteomics import mzml

# to visualize MS/MS and obtain ground truth
from pyteomics import pylab_aux as pa
from pyteomics import mass

# to find peaks
from scipy.signal import find_peaks

pep_loc = int(input("Index of peptide location: "))

# file that contains the associations between tandem spectra name and amino acids sequence
df = pd.read_csv("../data/psm_a.tsv", sep="\t")
#%%
# file that contains the amino acids and their weights
aa = pd.read_csv("../data/single_double_amino_acids.csv")
# beware that leucine (L) and isoleucine (I) weight the same and are indistinguishable
#%%
mz_path = '../data/2015-05-19_MRC5_a.mzML'
#%%
# this creates a dictionary of matches between the mzMl file and the psm_a.tsv
# this should only return MS/MS where we have the Peptide sequence already identified.
i = 0
matches = []
plotting_dict = None
with mzml.MzML(mz_path) as reader:
    for spectrum in reader:
        # looking for the first match between the 2 files
        for name in df["Spectrum"]:
            if name in spectrum.get('spectrum title') :
                # Extract relevant information
                matches.append([spectrum,df[df["Spectrum"] == name]["Peptide"].sum(),
                               df[df["Spectrum"] == name]["Hyperscore"].sum()])
m_df = (pd.DataFrame(matches, columns=["Spectrum","Peptide","Hyperscore"])
    .sort_values(by = ["Hyperscore"], axis=0, ascending= False)
    .reset_index(drop=True))

# collecting the data in variables
peptide = m_df["Peptide"].iloc[pep_loc]
print("Peptide sequence:", peptide)
mz_array = m_df["Spectrum"].iloc[pep_loc]['m/z array']
intensity_array = m_df["Spectrum"].iloc[pep_loc]['intensity array']
precursor_mz = m_df["Spectrum"].iloc[pep_loc]['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
plotting_dict = {'m/z array': mz_array, 'intensity array': intensity_array}
# I am subtracting 1.00782 because that is the mass of the extra proton that gives the precursor a double charge
# every other y ion has just one extra proton
original_peptide_mass = precursor_mz*2 - 1.00782

# converting to dataframe for simplicity?
s1 = pd.DataFrame(plotting_dict)

sp_height = intensity_array.max()*0.05
sp_dist = None
sp_prom = None
# avoiding early thresholding peter and kristoffer
# 5% of max height should be the cutoff threshold
# i.e.: based on how long we expect the amino acid sequence to be
# what determines wheter a peak is useful, is there someway we can find out?
#

sci_peak, _ = find_peaks(s1["intensity array"],
                      height=sp_height,
                      distance=sp_dist,
                      prominence=sp_prom
                      )


combo_df = aa
#%%
# the best threshold here would be 0.02 but I want to see how well the program works with 0.5 because I used 2.5 before
aa_thres = 0.05
def find_closest_aa(value, thres = aa_thres):
    """ Find the closest amino acid to a given mass.
    this code takes a dataframe with single amino acids and combinations of amino acids in this format

        full    letter  short   comp    mono    mass    G           A           S           P   V etc.
    0   glycine G       gly     c2h3no  57      57      g+g mass    g+a mass    g+s mass
    1   alanine A       ala     c3h5no  71      71      a+g mass    a+a mass    a+s mass
    etc.

    using this dataframe and a given mass, it will find the closest amino acid
     or combination of 2 amino acids to the given mass.
    It will only return amino acids that are within the threshold.
    It returns a list of lists, where each list contains the letter of the amino acid, its mass and the error.
    """
    if math.isnan(value) == True:
        return None
    single_df = combo_df.iloc[:,:-22]
    double_df = combo_df.set_index(["letter"]).iloc[
                :,[i for i in range(-22, -0)]]
    closeness_list = []
    # find the closest single amino acids
    loop_for_single = True
    while loop_for_single:
        aam_array = np.asarray(single_df["mono mass"])
        idx = (np.abs(aam_array - value)).argmin()
        error = np.abs(aam_array[idx] - value)
        if error > thres:
            loop_for_single = False
        else:
            name_idx = single_df["letter"].iloc[idx]
            closeness_list.append([name_idx, aam_array[idx], error])
            single_df = single_df.drop(single_df.index[idx])
    # find closest combination of amino acids
    loop_for_combo = True
    while loop_for_combo:
        error = (np.abs(double_df - value)).min().min()
        r, c = np.where(double_df == error + value)
        # if error wouldve been negative np.where will not find r, c
        # and pass empty arrays creating error
        if r.size == 0 :
            # print("boink")
            r, c = np.where(double_df == value - error)
        if error > thres:
            loop_for_combo = False
        else:
            name_idx = double_df.index[r[0]]+ double_df.columns[c[0]]
            closeness_list.append([name_idx, double_df.iloc[r[0],c[0]], error])
            double_df.iloc[r[0],c[0]] = None
    # print("closest aa is: ", name_idx, " ,with mass: ",
    # aam_array[idx], "Da. With an error of: ", error, "Da.")
    closeness_list.sort(key=lambda x: x[2])
    if closeness_list:
        # return closeness_list[0]
        return closeness_list[0][0]
    return None

mz_peaks = mz_array[sci_peak]
intensity_peaks = intensity_array[sci_peak]
# the precursor mass is not in the original mz_array or intensity_array dataset
# but since we know the precursor mass from the MS1 adding it makes our search easier
mz_peaks = np.append(mz_peaks,original_peptide_mass)
# I made the intensity of the manually added precursor max+1 so its the higest in the ranking
# this is because I think its always going to be part of the protein
intensity_peaks = np.append(intensity_peaks, intensity_peaks.max()+1)
# then another addition because at the C-terminus there is a water molecule and an extra proton
# it can be difficult to identify the last amino acid without using a little trick
# I will add another artificial peak at 19 Da to which the program can find a connected amino acid.
mz_peaks = np.append(19, mz_peaks)
# I had to append something at the 0th index and I choose the max intensity as well
# this time -0.5 instead of +1 because the +1 became the new max peak
# because the 19 Daltons should also always be there
intensity_peaks = np.append(intensity_peaks.max()-0.5,intensity_peaks)
mz_round = np.round(mz_peaks, decimals= 6)

# functions for automatically generating the groundtruth m/z values of peak from the known peptide sequence

def get_theoretical_y_ions(sequence):
    y_ions = []
    for i in range(1, len(sequence)+1):
        frag = sequence[-i:]
        mz = mass.fast_mass(frag, ion_type='y')
        y_ions.append(round(mz, 4)+1) # +1 because of extra proton charge
    return np.flip(y_ions).tolist()

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
#%%
real_sequence = []
real_peaks = []

for i in peptide:
    real_sequence.append(i)

for i in get_theoretical_y_ions(peptide):
    nearest = find_nearest(mz_round, i)
    real_peaks.append(nearest)


real_peptide = [real_sequence, real_peaks + [19]]

def rename(column):
    new = mz_round[int(column)]
    return new
#%%
# this code builds a distance matrix for the given mz_array
def build_dm():
    distance_matrix = []
    for column, i in zip(mz_round, range(len(mz_round))):
        distance_matrix.append([])
        for row in mz_round:
            distance_matrix[i].append(np.abs(row - column))

    df_dm = pd.DataFrame(distance_matrix)
    df_dm.columns = df_dm.columns.map(str)

    # for i in range(len(mz_round)):
    #     for j in range(i):
    #         df_dm.iat[i,j] = None
    for i in range(len(mz_round)):
        for j in range(len(mz_round)):
            if j >= i:
                df_dm.iat[i,j] = None
            if df_dm.iat[i,j] < 55 or df_dm.iat[i,j] > 480:
                df_dm.iat[i,j] = None

    # if j <= i:
    #     df_dm.iat[i,j] = None
    # drop columns with no values (all NaN)
    df_dm.dropna(axis = 1, how="all", inplace = True)

    print("starting mapping")
    df_aa = df_dm.map(find_closest_aa)
    print("mapping complete")


    #pavel pevzner
    df_aa.rename(columns=rename, inplace=True)
    df_aa.rename(index=rename, inplace=True)

    df_aa.to_csv(f"distance_matrix_{peptide}.csv")
    df_aa = pd.read_csv(f"distance_matrix_{peptide}.csv").set_index("Unnamed: 0")
    # suboptimal process but sadly I dont have the time to optimize this as I'd like.
    # I had this because I used the to_csv and read_csv functions when testing
    # and now I realized the csv thats being read is different to the one being written.
    # Because dictionaries get turned to strings.

    return df_aa

#%%
df_aa = build_dm()
#%% md
## function to find all possible paths through this distance matrix

#%%
# pathfinding algorithm based on depth first search
def dfs(cur_peak, aa_path, peak_path, all_paths, mz_round, df_aa):
    """
    Recursively explores all valid peptide paths using depth-first search.

    Parameters:
    - cur_peak: current m/z value we're at
    - aa_path: list of amino acids collected so far
    - peak_path: list of m/z values used in the path
    - all_paths: master list of all complete peptide paths found
    - mz_round: full list/array of m/z peaks
    - df_aa: distance matrix with AA matches (rows and cols = m/z values)
    """
    try:
        peak_pos = np.where(mz_round == cur_peak)[0][0]
    except IndexError:
        return  # skip if cur_peak not in mz_round

    found = False

    for col in df_aa.columns:
        aa = df_aa[col].iloc[peak_pos]
        # print("aa",aa)
        if isinstance(aa, str):

            aa_list = aa.split('+')  # support combos like "G+P"
            # print("aa_list",aa_list)
            for aa_letter in aa_list:
                next_peak = float(col)
                if next_peak in peak_path:
                    continue  # avoid cycles

                # recurse to explore next peak
                dfs(
                    next_peak,
                    aa_path + [aa_letter],
                    peak_path + [next_peak],
                    all_paths,
                    mz_round,
                    df_aa
                )
                found = True

    if not found:
        # Dead end: store complete path
        all_paths.append((aa_path, peak_path))

all_paths = []
start_peak = float(mz_round[-1])  # typically highest m/z
dfs(start_peak, [], [start_peak], all_paths, mz_round, df_aa)

pep_df = pd.DataFrame(all_paths, columns = ["sequence","peaks"])
pep_df['len_pep'] = pep_df["peaks"].str.len()
pep_df.sort_values(by='len_pep', ascending=False)

ranked_intensity = np.flip(np.sort(intensity_peaks))

def intensity_rank(peaks_list):
    rank_sum = 0
    for i in peaks_list:
        intensity = intensity_peaks[np.where(mz_round == i)[0][0]]
        # print(intensity, np.where(ranked_intensity == intensity)[0][0])
        # I added the + 1 so that the best rank is 1
        # I did this because the ChatGPT code did this too and I didnt want to mess it up lol
        rank_sum += np.where(ranked_intensity == intensity)[0][0] + 1
    return rank_sum
#%%
def simulate_rank_sum_distribution(N, T, num_samples=1000000):
    rank_sums = []
    ranks = list(range(1, N+1))
    for _ in range(num_samples):
        sample = random.sample(ranks, T)
        rank_sums.append(sum(sample))
    return rank_sums

def empirical_p_value(observed_sum, distribution):
    # Probability of getting an equal or smaller sum by chance
    count = sum(1 for s in distribution if s <= observed_sum)
    return count / len(distribution)
#%%
# creating random rank distributions for each observed peptide length in the dataset

pep_lens_in_list = range(pep_df["len_pep"].min(),pep_df["len_pep"].max()+1)
rank_distributions = {}
for T in pep_lens_in_list:
    rank_distributions[T] = simulate_rank_sum_distribution(len(intensity_peaks), T)
#%%
# applying the intensity rank functions to the entire dataframe

pep_df["rank_sum"] = pep_df["peaks"].apply(intensity_rank)
pep_df["rank_p"] = pep_df.apply(
lambda row: empirical_p_value(row["rank_sum"], rank_distributions[row["len_pep"]]),
axis=1
)
#%% md
## 2nd m/z fidelity
#%%
def find_SSE(sequence, peaks):
    # find the m/z fidelity from the first peak to the last
    # important here our first peak is the first y-ion so its already missing the first peptide
    estimation = []
    distance = 0
    for i in range(1,len(sequence)+1):
        if len(sequence[-i]) == 1:
            distance += aa[aa["letter"]==sequence[-i]]["mono mass"].sum()
        if len(sequence[-i]) == 2:
            distance += aa[aa["letter"]==sequence[-i][0]]["mono mass"].sum()
            distance += aa[aa["letter"]==sequence[-i][1]]["mono mass"].sum()
        estimation.append(peaks[-(i+1)] - distance)
    mean = np.mean(estimation)
    estimation = np.array(estimation)
    SSE = np.sum((estimation - mean) ** 2)

    return SSE


# apply SSE to dataframe
pep_df["SSE"] = pep_df.apply(
lambda row: find_SSE(row["sequence"], row["peaks"]),
axis=1
)

def create_SSE_dist(n):
    SSE_distributions = {}
    dist = np.random.uniform(10-aa_thres,10+aa_thres,1000)
    for T in pep_lens_in_list:
        single_dist = []
        for i in range(n):
            estimation = []
            for i in range(T):
                estimation.append(random.choice(dist))
            mean = np.mean(estimation)
            estimation = np.array(estimation)
            SSE = np.sum((estimation - mean) ** 2)
            single_dist.append(SSE)
        SSE_distributions[T] = single_dist
    return SSE_distributions

SSE_len = 100000
SSE_distributions = create_SSE_dist(SSE_len)
def find_mz_fid_p(row):
    count = sum(1 for s in SSE_distributions[row["len_pep"]] if s <= row["SSE"])
    return (count / SSE_len) + 0.000001 # to avoid math errors when value is really low

# apply mz fidelity to dataframe

pep_df["mz_fid_p"] = pep_df.apply(
lambda row: find_mz_fid_p(row),
axis=1
)

def fisher_combined_p(p1, p2):
    """
    Combine two p-values using Fisher's method.
    """
    chi_stat = -2 * (math.log(p1) + math.log(p2))
    combined_p = chi2.sf(chi_stat, df=4)  # df = 2 * number of p-values
    return combined_p

# apply fisher combined to dataframe

pep_df["combined_p"] = pep_df.apply(
lambda row: fisher_combined_p(row["rank_p"],row["mz_fid_p"]),
axis=1
)

# I assume the best fit peptide is in the top 20

top_20 = pep_df.sort_values(by=["combined_p"], ascending=True).reset_index()[0:20]

def num_doubles(str_list):
    counter = 0
    for i in str_list:
        if len(i) == 2:
            counter += 1
    return counter


# read double amino acids in a sequence

top_20["doubles"] = top_20.apply(
lambda row: num_doubles(row["sequence"]),
axis=1
)

# join the sequence list into a string

top_20["sequence"] = top_20.apply(
lambda row: "".join(row["sequence"]),
axis=1
)

counted_final = {i:top_20["sequence"].tolist().count(i) for i in top_20["sequence"].tolist()}
Best_fit_peptide = (top_20[top_20["sequence"] ==
       max(counted_final, key=counted_final.get)]
        .sort_values(by=["doubles"], ascending=True)
        .reset_index().iloc[0]
)

def I_to_L(sequence):
    for i in range(len(sequence)):
        if sequence[i] == "I":
            sequence[i] = "L"
    return sequence

real_peptide[0] = I_to_L(real_peptide[0])

def equality_check(a, b):
    equality_rank = 0
    if len(a) == len(b):
        for x, y in zip(a, b):
            if x == y:
                equality_rank += 1
            if x == y:
                pass
    else:
        print("peptides dont have same len")
    return equality_rank/len(a)


equality_score = equality_check(Best_fit_peptide["sequence"], "".join(real_peptide[0]))

print("equality score", equality_score)
print("real peptide",real_peptide)
print("Top 20:\n", top_20)