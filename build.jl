include("csv2table.jl")

maketable("results\\all_accuracy","rcccccc","Dataset & Logit Train & Logit Test & MCMC Train & MCMC Test & MFVB Train & MFVB Test",
  "all_acc","Training and test set accuracy for logistic regression, MCMC logistic regression, and MFVB logistic regression on all datasets.
  The highest out-of-sample accuracy score is highlighted for each dataset.  \\ *** : The training dataset was completely separable;
  therefore the posterior covariance matrix was zero and the precision matrix was undefined.  
  Thus, we did not obtain MCMC predictive probabilities for this dataset.",
  headerPresent=true)

