import pandas as pd
import numpy as np

def read_pickle_file(file):
    pickle_data = pd.read_pickle(file)
    return pickle_data

def write_pickle_df(df, file):
    df.to_pickle(file)

def read_npy(file):
    return np.load(file).item()
