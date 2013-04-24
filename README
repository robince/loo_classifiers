loo\_classifiers
================

This is a Matlab package to efficiently implement leave-one-out cross-validation for a number of basic supervised learning algorithms.

I was surprised I couldn't find anything online to do this, so perhaps this will be useful to other people.

Since this was developed specifically for my work, studying neural coding, it makes some strong assumptions - namely balanced data (same number of trials for every class) and assumes uniform priors over the classes.

All the algorithms take as the first argument a 3d data array with dimensions `(Nftr, Ncls, Ntrl)`, and returns two arguments - the confusion matrix `(Ncls, Ncls)` and the information in the confusion matrix. To calculate the information value you will need the [InfoToolbox](http://www.ibtb.org/) package installed. Some classifiers have extra parameters which should be passed as additional input arguments.

Where possible, results were tested against Matlab `classify` from the Statistics Toolbox.

Algorithms
----------

The currently implemented algorithms are:

* `diag_linear` : MVN clusters, pooled diagonal covariance matrix
* `linear` : MVN clusters, pooled covariance matrix
* `diag_quadratic` : MVN clusters, diagonal covariance matrices
* `quadratic` : MVN clusters
* `nearest_mean` : MVN clusters, pooled diagonal covariance, equal variances (template matching)
* `poisson_bayes` : Independent poisson features for each cluster. *Count data only*
* `multinom_bayes` : Independent multinomial features for each cluster. *Count data only*
* `knn` : k-Nearest neighbour classifier.

Contact
-------

Bugs / comments / ideas / criticism to robince at gmail dt com

This is unpublished research code so comes with no guarentees. Please feel free to try it out, but if it is useful to you and forms a key part of an analysis please consider contacting me. Perhaps I could help in some way with feedback or modifications or optimisations specific to your case. I don't have any citation yet to provide for this but will update here if/when I do. 

License
-------

This project is licensed under version 3 of the GNU General Public License. For the exact terms please see the [LICENSE file](https://github.com/robince/loo_classifiers/blob/master/LICENSE).

ToDo
----

- update to use MatlabAPI\_lite
- change interface to more conventional seperate feature and class indicator variables + allow non-uniform classes
