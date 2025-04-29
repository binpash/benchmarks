#!/usr/bin/env python3

from sklearn.model_selection import train_test_split
from sklearn import datasets
import pickle
import os


raw_data = datasets.fetch_covtype(data_home="inputs", download_if_missing=False)
    
data = train_test_split(raw_data.data, 
                        raw_data.target, 
                        test_size=0.2, 
                        random_state=0)
filenames = ['X_train', 'X_test', 'y_train', 'y_test']
tmp = os.environ.get('TMP')
filepath = os.path.join(tmp, 'model.obj')
for datum, name in zip(data, filenames):
    filepath = os.path.join(tmp, f'{name}.obj')
    with open(filepath, 'w+b') as file:
        pickle.dump(datum, file)
