#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import argparse
import os
import json

from keras.models import Sequential
from keras.layers import Dense
from keras.optimizers import Adam

from load_data import DataLoader, split_data
from dnn_config import NODES, EPOCHS, BATCHSIZE, OBSERVABLES

### define dnn structure with 3 hidden layers ###
def build_model():
    opt = Adam(learning_rate=0.0005)
    model = Sequential()
    model.add(Dense(NODES, input_dim=len(OBSERVABLES), activation='relu'))
    model.add(Dense(NODES, activation='relu'))
    model.add(Dense(NODES, activation='relu'))
    model.add(Dense(5, activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer=opt, weighted_metrics=['accuracy'])
    model.summary()
    return model

### train dnn ###
def train_model(model, train, valid):
    x_train, y_train, w_train = split_data(train)
    x_valid, y_valid, w_valid = split_data(valid)

    history = model.fit(
        x_train, 
        y_train, 
        sample_weight=w_train,
        epochs=EPOCHS, 
        batch_size=BATCHSIZE, 
        validation_data=(x_valid, y_valid, w_valid)
    )

    return model, history.history

### save dnn ###
def export(model, name):
    os.system(f"rm {name}.keras")
    model.save(f"{name}.keras")

### save training history ###
def save(hist, name):
    with open(f"{name}.json", "w") as file:
        json.dump(hist, file)
 

### main function - load data, build model, train dnn, export ###
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("signal", default="", help="eg AZH_1000_400")
    parser.add_argument("kfold", default=1, help="enter kfold value (eq 1,2,3...)")

    args = parser.parse_args()
    signal = args.signal
    k = int(args.kfold)

    loader = DataLoader(signal, k, shuffeled=True)
    print(loader.full_set)
    train = loader.train_set
    valid = loader.valid_set

    dnn = build_model()
    dnn, history = train_model(dnn, train, valid)
    name = f'dnn_{signal}_kfold{k}'.replace("AZH_","")
    export(dnn, name)
    save(history, name)