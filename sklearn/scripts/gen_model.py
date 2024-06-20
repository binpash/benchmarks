from sklearn.linear_model import LogisticRegression
import pickle
import sys
import os

reg = LogisticRegression(max_iter=int(sys.argv[1]), 
                         solver='newton-cholesky', 
                         multi_class='ovr')
with open(f'{os.environ.get("TMP","./tmp")}/model.obj', 'w+b') as file:
    pickle.dump(reg, file)
