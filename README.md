perturb R package
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
The perturb R package evaluates collinearity by adding small random values (perturbations) to specified variables in a model, re-estimating the model, and repeating this process e.g. 100 times. The model estimates can then be evaluated to determine how much these random perturbations have affected the results. The perturbances are applied to the untransformed variables or to the variables included in an interaction effect, making it possible to evaluate how small changes are propagated to derived model terms. This principle can be extended to categorical variables where a probability of misclassification is specified. The perturb package can also be used in generalized linear models such as logistic regression.

Version 3 of perturb can also be used with mixed models estimated by the `nlme` or `lme4` packages. This can take into account that a variable for a random slope is correlated with other variables in the mixed model. Or the impact of changes to a levels variable such as the site where data was collected. The site variable often contains small categories, applying random reclassifications can indicate how robust the results are.
