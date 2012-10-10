#!/usr/bin/python
# -*- coding: utf-8 -*-
#


import numpy as np
import scipy as sp
import pylab as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
import interval
import itertools
import networkx as NX
import ffnet  # ffnet, mlgraph, savenet, loadnet, exportnet
from numpy.random import RandomState
from itertools import cycle
from sklearn import linear_model as lm  # LogisticRegression, SGDClassifier
from sklearn import svm  # SVC, NuSVC, LinearSVC, SVR
from sklearn import neighbors as ngb  # NeighborsClassifier, NeighborsRegressor
from sklearn import mixture as mix  # GMM
from sklearn import cluster as clus  # KMeans, MeanShift, AffinityPropagation, SpectralClustering


def ffnetwork(inputs=1, layers=1, outputs=1):
    """
    Define a feedforward neural network

    Parameters
    ----------
    inputs : integer
        default = 1
        defines the length of input vector
    layers : integer
        default = 1
        defines the number of layers
    outputs : integer
        default = 1
        defines the length of output vector

    Returns
    -------
    net : Feed-forward neural network
    """
    conec = ffnet.mlgraph((inputs, layers, outputs), biases=False)
    net = ffnet.ffnet(conec)
    '''NX.draw_graphviz(net.graph, prog='dot')
    plt.show()'''
    return net


def knn_learn(nneighbors=1, data_train=np.array([]), target_train=np.array([]), data_test=np.array([])):
    """
    estimate position using the K neareast neighbors (KNN) technique

    Parameters
    ----------

    nneighbors : int
        default = 1
    data_train : numpy.ndarray
        default = array([])
    target_train : numpy.ndarray
        default = array([])
    data_test : numpy.ndarray
        default = array([])

    Returns
    -------
    targets : numpy.ndarray

    """

    clf = ngb.NeighborsClassifier(nneighbors).fit(data_train, target_train)
    targets = clf.predict(data_test)
    return targets


def ffann_learn(layers=1, data_train=np.array([]), target_train=np.array([]), data_test=np.array([])):
    """
    estimate position using the feed-forward ANN technique

    Parameters
    ----------

    layers : int
        default = 1
    data_train : numpy.ndarray
        default = array([])
    target_train : numpy.ndarray
        default = array([])
    data_test : numpy.ndarray
        default = array([])

    Returns
    -------
    targets : numpy.ndarray

    """

    net = ffnetwork(np.shape(data_train)[1], layers, np.shape(data_train)[1])
    target_train = target_train * 1.0
    net.train_tnc(data_train, target_train)
    targets = net(data_test)
    return targets

def svm_learn(kernel='linear', data_train=np.array([]), target_train=np.array([]), data_test=np.array([])):
    """
    estimate position using the svm techniques

    Parameters
    ----------

    kernel : string
        default = 'linear'
        kernel can be 'linear', 'poly', 'rbf', 'sigmoid'
    data_train : numpy.ndarray
        default = array([])
    target_train : numpy.ndarray
        default = array([])
    data_test : numpy.ndarray
        default = array([])

    Returns
    -------
    targets : numpy.ndarray

    """

    clf = svm.SVC(kernel).fit(data_train, target_train)
    targets = clf.predict(data_test)
    return targets


def lreg_learn(data_train=np.array([]), target_train=np.array([]), data_test=np.array([])):
    """
    estimate position using logistic regression

    Parameters
    ----------

    data_train : numpy.ndarray
        default = array([])
    target_train : numpy.ndarray
        default = array([])
    data_test : numpy.ndarray
        default = array([])

    Returns
    -------
    targets : numpy.ndarray

    """

    clf = lm.LogisticRegression().fit(data_train, target_train)
    targets = clf.predict(data_test)
    return targets
