from sklearn.model_selection import train_test_split
from sklearn import datasets
import pickle
import os

tmp = os.environ.get('TMP')
raw_data = datasets.fetch_covtype(data_home=f"{os.getenv('TMP')}", download_if_missing=True)
    
data = train_test_split(raw_data.data, 
                        raw_data.target, 
                        test_size=0.2, 
                        random_state=0)
filenames = ['X_train', 'X_test', 'y_train', 'y_test']
filepath = os.path.join(tmp, 'model.obj')
for datum, name in zip(data, filenames):
    filepath = os.path.join(tmp, f'{name}.obj')
    with open(filepath, 'w+b') as file:
        pickle.dump(datum, file)
