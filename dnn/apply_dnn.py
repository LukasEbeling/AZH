#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import argparse
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from keras.models import load_model
from load_data import DataLoader, split_data
from dnn_config import TARGETS, CACHE

### use trained dnn to evaluate data ###
### eg ./apply_dnn
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("signal", default="", help="eg AZH_1000_400")
    parser.add_argument("variation", default="nominal", help="JES and JER variation")    

    args = parser.parse_args()
    signal = args.signal
    variation = args.variation

    loader = DataLoader(signal,1,variation)
    data = loader.full_set 

    positive_class = 0

    for k in [1,2,3,4,5]:
        evalu = data[data['kfold'] == k].copy()
        x,y,foo = split_data(evalu) 
        w = data["event_weight"]  
        model = load_model(f"dnn_{signal}_kfold{k}.keras".replace("AZH_",""))
        output = model.predict(x)
        score = [arr[positive_class]/sum(arr) for arr in output]
        node = np.argmax(output, axis=1)

        mask = (data['kfold'] == k)
        data.loc[mask, 'score'] = score
        data.loc[mask, 'node'] = node

    processes =  list(TARGETS.keys()) + [signal]
    
    for process in processes:
        mask = (data['process'] == process)
        pa_table = pa.table({"foo": np.array(data.loc[mask,'score'])})
        pq.write_table(pa_table, CACHE+f'/mc_UL17_{process}_{variation}_score.parquet')
        pa_table = pa.table({"foo": np.array(data.loc[mask,'node'])})
        pq.write_table(pa_table, CACHE+f'/mc_UL17_{process}_{variation}_node.parquet')
        

        



