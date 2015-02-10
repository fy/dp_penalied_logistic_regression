run_cvx <- function(XX, yy, lambda, alpha, noise, noise_scale){
  ## Save data to Matlab data file
  writeMat('cvx_data.mat', XX=XX, yy=yy, lambda=lambda, alpha_var=alpha,
           noise=noise, noise_scale=noise_scale);
  warnings()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Use MATLAB server
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  evaluate(matlab, "opt_elastic_net_for_R()")

  ## Load results from Matlab
  fit_beta = unlist(read.csv("cvx_results.csv", header=FALSE))

  return(list(estimate=fit_beta))
}
