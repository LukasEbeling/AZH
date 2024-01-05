#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python

from dnn_config import CACHE, TARGETS, LOAD, OBSERVABLES

import pandas as pd
import numpy as np

class DataLoader():

    ### load data set for given signal and variation ###
    def __init__(self, signal: str, k: int, variation: str = 'nominal', shuffeled: bool = False):

        self.signal = signal
        self.variation = variation

        data = self.load()
        if shuffeled: data = self.shuffel(data)

        self.full_set = data
        self.evalu_set = data[(data['kfold'] - k)%5 == 0]
        self.valid_set = data[(data['kfold'] - k)%5 == 1]
        self.train_set = data[(data['kfold'] - k)%5 >= 2]
    
    ### reset k-fold value and split data into train, evalu, and validation set ###
    def reset(k: int):
        self.evalu_set = data[(data['kfold'] - k)%5 == 0]
        self.valid_set = data[(data['kfold'] - k)%5 == 1]
        self.train_set = data[(data['kfold'] - k)%5 >= 2]

    ### load all needed data as data frame ###
    def load(self):
        
        df_list = []

        processes = list(TARGETS.keys()) + [self.signal]

        for process in processes:
            print(process)
            df = self.load_process(process)
            df = self.reweight(df)
            df = self.calculate(df)
            df_list.append(df)

        data = pd.concat(df_list)
        data = self.normalize(data)
        return data

    ### load data of a single process as pandas data frame ###
    def load_process(self, process: str):
        process_dict = {"process": process}

        for obs in LOAD: 
            path = CACHE + f'/mc_UL17_{process}_{self.variation}_{obs}.parquet'
            data = self.load_parquet(path)
            process_dict[obs] = data

        if 'AZH' in process: targets = [1,0,0,0,0]
        else: targets = TARGETS[process]
        for i,t in enumerate(targets): process_dict[f"t{i}"] = t
        return pd.DataFrame(process_dict)

    ### load parquet file (with headline 'foo') ###
    def load_parquet(self, path):
        df = pd.read_parquet(path, columns=['foo'])
        arr = df.to_numpy()
        return(np.array([subarray[0] for subarray in arr]))

    ### up-scale weight of rare processes (eg signal) ###
    def reweight(self, df):
        if 'event_weight' not in df.columns: return false
        df = df
        df['weight'] = df['event_weight']/sum(df['event_weight'])
        return df

    ### calculate additional variables ###
    def calculate(self, df): 
        df = df

        if 'delta_phi_12' in OBSERVABLES:
            df["delta_phi_12"] = abs(np.array(df["jetsAk4CHS.m_phi_1"])-np.array(df["jetsAk4CHS.m_phi_2"]))

        if 'delta_phi_13' in OBSERVABLES:  
            df["delta_phi_13"] = abs(np.array(df["jetsAk4CHS.m_phi_1"])-np.array(df["jetsAk4CHS.m_phi_3"]))
        
        if 'delta_phi_23' in OBSERVABLES:
            df["delta_phi_23"] = abs(np.array(df["jetsAk4CHS.m_phi_2"])-np.array(df["jetsAk4CHS.m_phi_3"]))

        for i in [1,2,3]:
            if f"bscore{i}" not in OBSERVABLES: continue
            df[f"bscore{i}"] = df[f"jetsAk4CHS.m_btag_DeepFlavour_probbb_{i}"]
            df[f"bscore{i}"] += df[f"jetsAk4CHS.m_btag_DeepFlavour_probb_{i}"]
            df[f"bscore{i}"] += df[f"jetsAk4CHS.m_btag_DeepFlavour_problepb_{i}"]

        return df

    ### normalize values of observables to ensure comparability ###
    def normalize(self, df):
        for obs in OBSERVABLES:
            #df[obs] = df[obs] + min(df[obs])
            #df[obs] = df[obs]/np.mean(df[obs])
            df[obs] = (df[obs] - np.mean(df[obs]))/np.std(df[obs])

        return df
    
    ### shuffel data frame ###
    def shuffel(self, df):
        return df.sample(frac=1).reset_index(drop=True)

### split data (df format) into observables (x), targets (y) and weights (w) ###
def split_data(dataset):
    x = np.array([np.array(dataset[o]) for o in OBSERVABLES]).T
    y = np.array([np.array(dataset[t]) for t in ["t0","t1","t2","t3","t4"]]).T
    #w = dataset["event_weight"]    
    w = dataset["weight"]    
    return x,y,w   

    



    