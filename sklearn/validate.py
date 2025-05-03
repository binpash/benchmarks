#!/usr/bin/env python

from sklearn.model_selection import train_test_split
from sklearn import datasets
from sklearn.linear_model import LogisticRegression
import pickle
import numpy as np

max_iter = 100

dataset = datasets.fetch_covtype(data_home="inputs", download_if_missing=False)

X_train, X_test, y_train, y_test = train_test_split(dataset.data, dataset.target, 
                                                    test_size=0.2, 
                                                    random_state=0)

control_model = LogisticRegression(max_iter=max_iter, 
                                   solver='newton-cholesky', 
                                   multi_class='ovr')
control_model.fit(X_train, y_train)
control_score = control_model.score(X_test, y_test)

with open('tmp/trained_model.obj', 'rb') as file:
    experiment_model = pickle.load(file)
experiment_score = experiment_model.score(X_test, y_test)

assert experiment_score == control_score
assert np.array_equal(control_model.coef_, experiment_model.coef_)