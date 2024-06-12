from sklearn.model_selection import train_test_split
from sklearn import datasets
import pickle
import numpy as np

raw_data = datasets.fetch_rcv1(data_home="inputs", download_if_missing=False)
    
data = train_test_split(raw_data.data, 
                        np.argmax(raw_data.target.toarray()), 
                        test_size=0.2, 
                        random_state=0)
filenames = ['X_train', 'X_test', 'y_train', 'y_test']
for datum, name in zip(data, filenames):
    with open(f'./tmp/{name}.obj', 'w+b') as file:
        pickle.dump(datum, file)
