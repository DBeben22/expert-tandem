

def find_first_path(mz_round, intensity_peaks, df_aa):
    # This is written to work form the highest peak down to the lowest similar to how its done by nick wedd.
    # however my distance matrix removed all points below the diagonal meaning anytime I identify a hit of an amino acid it looks like follows
    # hit a row: x and column: y. where y > x and since my program iterates over rows and then columns it sets the column to be the new row.
    # which means it starts back at the top and repeats itself in a pattern. Ill have to change the orientation of the matrix.
    path_one = []
    cur_peak = float(mz_round[-1])
    #
    while cur_peak > float(mz_round[0]):
        found = False
        for i in reversed(df_aa.columns):
            # print("conk", type(aa))
            peak_pos = np.where(mz_round == cur_peak)[0][0]
            aa_inf = df_aa[i].iloc[peak_pos]#.values[0]
            aa_int = intensity_peaks[peak_pos]
            print("row:",mz_round[peak_pos],"column: ", i,"aa:",aa_inf,"int:",aa_int)
            #time.sleep(0.1)
            # choosing the first possible connection always as of now
            if isinstance(aa_inf, str) :
                print("\n boingoloingo \n \n")
                path_one.append(aa_inf)
                # ast.literal_eval is a security risk but I dont know a better way
                # Problem is that dictionaries get saved as string when .to_csv is applied
                found = True
                cur_peak = float(i)
        if not found:
            print("jump to next row")
            cur_peak = float(mz_round[peak_pos-1])
            found = False

    print(path_one)


def find_closest_aa(value, thres = 5):
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
    # TODO remove useless calculations from this like closeness_list having values that I never use
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
            name_idx = double_df.index[r[0]]+ "+" + double_df.columns[c[0]]
            closeness_list.append([name_idx, double_df.iloc[r[0],c[0]], error])
            double_df.iloc[r[0],c[0]] = None
    # print("closest aa is: ", name_idx, " ,with mass: ",
    # aam_array[idx], "Da. With an error of: ", error, "Da.")
    closeness_list.sort(key=lambda x: x[2])
    if closeness_list:
        # return closeness_list[0]
        return closeness_list[0][0]
    return None
