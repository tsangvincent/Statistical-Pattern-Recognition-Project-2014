# Statistical-Pattern-Recognition-Project-2015

The dataset ”Data 20” from ”http://wwwf.imperial.ac.uk/ eakc07/S7/data20.dat”
contains 1103 observations with 28 features and their corresponding class indicators.

The set of features to be measured is very large, which requires advanced methodology to build the prediction model. Therefore, it could potentially be time consuming and computational expensive. 

I have also loaded the packages that might be useful for my classifications for later sections. The term ”Curse of Dimensionality which is first mentioned by Bellman has particularly stressed on the problem when dealing with high dimensional data requiring large sets of training data.

(A). Data classifiers:
I have chosen 4 data classification methods based on the training set data to build the corresponding prediction model. The 4 chosen classifiers are:
1. K-nearest Neighbour (KNN)
2. Quadratic Discriminant Analysis (QDA)
3. Classification and Regression Tree (Tree-based) 
4. Multilayer Perception (MLP)

(B). Feature Selections:
There are 2 approaches:
1. Bottom-up (Sequential forward selection): Start from a null feature set. We are required to find out the individual best-performing feature, which has the lowest error rate. Add feature one at a time and find out the best- performing combination of features. Stop when the best performance of the new combination worsen off.
2. Top-down (Sequential backward selection): Start from a full feature set with n features. We are required to delete one feature at a time and find out the best-performing combinations of n-1 features, which has the lowest error rate. Stop when the best performance of the new combination worsen off.

(C). Performance Assessment:
10-fold cross-validation
