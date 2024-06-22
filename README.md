# K-prototypes-analysis
K-prototype clustering method is an extension to k-means method that accounts for both categorical and numerical data. The parameter lambda defines the trade-off between Euclidean distance of numeric variables and simply matching coefficient between categorical variables. Feature Importance can be used as an additional guide to tune this parameter. Briefly, the algorithm finds the variables which drive the cluster assignment using Gower distance and scores them according to their relevance (Huang, 1997). If X and Y are observations in a mixed dataset. X represents a vector with numerical variables (x1, x2, …xk) and Y represents a vector with categorical variables (y1, y2, ….., yp). The dissimilarity (D) between these two observations can be measured as following:

D(X,Y)=∑_(i=1)^k▒〖(x_i-y_i )²+y ∑_(i=k+1)^p▒(x_i,y_i ) 〗        (1)


The first term in the equation estimates the Euclidean distance between the observations X and Y considering ith numerical variables and the second term estimates the dissimilarity between X and Y considering ith categorical variables and the accuracy of the model was calculated using actual vs predicted values of the validation set (Huang, 1997).The procedure of k-prototypes clustering algorithm is as follows: 
	Choose the number of ’k’ clusters.
	Create the space represented by the clustered objects which are represented in the initial group of centroids.
	Measure the distance between new objectives to each group centroid, and then aggregate the point to a closest center.
	Recalculate the group centroid based on classified objectives. 
	Repeat Steps 2 and 3 until the centroids does not change. 

