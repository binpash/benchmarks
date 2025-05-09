from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
from sklearn import datasets
import pickle
import pandas as pd
import numpy as np
import os

X, y = datasets.fetch_kddcup99(data_home=f"{os.environ.get('DATA')}", percent10=False, return_X_y=True, as_frame=True, download_if_missing=True)
X = pd.DataFrame(X).drop(columns=["protocol_type", "service", "flag"]).astype(float)
X[X.columns] = MinMaxScaler().fit_transform(X[X.columns])
X = X.to_numpy()
y = LabelEncoder().fit_transform(y).astype(np.int32)

data = train_test_split(X, 
                        y,
                        test_size=0.2, 
                        random_state=0)
filenames = ['X_train', 'X_test', 'y_train', 'y_test']
for datum, name in zip(data, filenames):
    with open(f'{os.environ.get("TMP","./tmp")}/{name}.obj', 'w+b') as file:
        pickle.dump(datum, file)
