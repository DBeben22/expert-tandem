import numpy as np
import pandas as pd
from Final_vers_function import create_fit_df


m_df = pd.read_csv("../data/PSM_df.csv")
r_list = np.random.randint(0,len(m_df),80)
for index in r_list:
    create_fit_df(index)

with open("../data/r_list_2.txt", "w") as output:
    output.write(str(r_list))
