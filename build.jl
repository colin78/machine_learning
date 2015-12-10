include("csv2table.jl")

maketable("results\\all_accuracy","rcccccc","Dataset & Logit Train & Logit Test & MCMC Train & MCMC Test & MFVB Train & MFVB Test",
  "all_acc","Training and test set accuracy for logistic regression, MCMC logistic regression, and MFVB logistic regression on all datasets.",
  headerPresent=true)

