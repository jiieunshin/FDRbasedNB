### FDR-based categorical variables selection in Naive Bayes classification
---

## Usage Example
```
# fix a seed and generate an example dataset (see in paper)
set.seed(123)
tr = sim_data(N_tr = 200, p = 20)
te = sim_data(N_te - 1000, p = 20)

# fitting for training data
fit = nbayes(tr$X, tr$Y)

# predict test data
fit_pred = predict_nbayes(te$X, te$Y, fit$prob, fit$prior)$miss

# compair the result between a true Y and the predicted Y
table(te$Y, fit_pred$miss)
```
paper URL (only Korean): https://doi.org/10.7465/jkdi.2021.32.6.1329
