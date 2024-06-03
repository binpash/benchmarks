from sklearn.model_selection import train_test_split
from sklearn import datasets
import pickle


raw_data = datasets.fetch_covtype(data_home="inputs", download_if_missing=False)
    
data = train_test_split(raw_data.data, 
                        raw_data.target, 
                        test_size=0.2, 
                        random_state=0)
filenames = ['X_train', 'X_test', 'y_train', 'y_test']
for datum, name in zip(data, filenames):
    with open(f'./tmp/{name}.obj', 'w+b') as file:
        pickle.dump(datum, file)
