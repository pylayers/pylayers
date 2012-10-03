from pylayers.location.learning.loclearn import *

data_train = np.array([[1.0,2.3,3.5],[5.2,6.3,7.2],[5.4,6.8,9.5]])
target_train = np.array([[2.27],[6.23],[7.23]])
data_test = np.array([[2.0,3.3,3.5],[4.2,6.3,8.2],[9.4,6.8,5.5]])

print 'KNN'
print knn_learn(1, data_train, target_train, data_test)
print 'FFANN'
print ffann_learn(1, data_train, target_train, data_test)
print 'SVM Linear'
print svm_learn('linear', data_train, target_train, data_test)
print 'SVM Poly'
print svm_learn('poly', data_train, target_train, data_test)
print 'SVM rbf'
print svm_learn('rbf', data_train, target_train, data_test)
print 'SVM sigmoid'
print svm_learn('sigmoid', data_train, target_train, data_test)
print 'Logistic regression'
print lreg_learn(data_train, target_train, data_test)
