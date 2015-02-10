## This function performs differentially private elastic-net penalized logistic 
## regression multiple times and saves the results to an RData file. It also
## performs non-private penalized logistic regression and saves the results to 
## an RData file.

###############################################################################
## load data
###############################################################################

source('./analyze_hapsample.R')



###############################################################################
## functions
###############################################################################
obj.f = function(beta, XX, yy, lambda, alpha, noise, noise_scale) {
  ## object function of penalized logistic regression
  ## (with intercept)
  NN = nrow(XX)
  ss = 0
  for (ii in 1:nrow(XX)) {
    ss = ss + yy[ii] * sum(XX[ii, ] * beta) - log(1 + exp(sum(XX[ii, ] *
                                                                beta)))
  }
  ## object function of elastic-net penalized logistic regression
  -1 / NN * ss + lambda / 2 * (1 - alpha) * sum(beta[-1]^2) + lambda * alpha *
    sum(sapply(beta[-1], abs)) + noise_scale * sum(noise[-1] * beta[-1])
}

get_design_matrix = function(MM) {
  ## get the design matrix of the pair-wise interaction regression model
  top_M_snp = valid_snps[order(valid_snp_chisq, decreasing=TRUE)[1:MM]]
  top_snp_data = training_data[, match(top_M_snp, names(training_data))]
  pairwise_interact_formula = sprintf("~(%s)^2",paste(names(top_snp_data),
                                                      collapse='+'))
  design_matrix = model.matrix(formula(pairwise_interact_formula),
                               top_snp_data)
  design_matrix = design_matrix[, -1] ## remove intercept
  design_matrix
}

get_kappa = function(MM, norm=2) {
  ## Get upper bound of the the p-norm of the gradient,
  ## and the gradient is of the form form xx^T.
  if (norm==2) {
    kappa = sqrt(8 * MM^2 - 4 * MM)
  }
  if (norm==1) {
    kappa = 2 * MM^2
  }
  kappa
}

rB2 = function(deg) {
  ## perturbation noise function (l2-norm)
  WW = rnorm(deg)
  YY = rchisq(1, 2*deg)
  YY * WW / sqrt(sum(WW^2))
}




###############################################################################
## simulation parameters
###############################################################################
NN = nrow(training_data)

nn_sim = 100 ## number of simulations
nn_lambda = 5 ## number of regularization parameters

MM_vec = c(5) ## top M SNPs
epsilon_vec = c(0.5, 1, 5, 10)
alpha_vec = c(0.1, 0.5, 0.9)

param_list = lapply(MM_vec, function(MM) {
  kappa = get_kappa(MM)
  phi = 2 * kappa
  res = lapply(epsilon_vec, function(epsilon) {
    lambda_convex_min = kappa^2 / (NN *  (exp(epsilon / 4) - 1))
    return(list(kappa=kappa,
                phi=phi,
                lambda_convex_min=lambda_convex_min,
                epsilon=epsilon,
                MM=MM))
  })
  names(res) = sapply(epsilon_vec, paste)
  return(res)
})
names(param_list) = sapply(MM_vec, paste)

param_dtf = c()
for (MM_res in param_list) {
  for (epsilon_res in MM_res) {
    param_dtf = rbind(param_dtf, c(kappa=epsilon_res$kappa,
                        phi=epsilon_res$phi,
                        lambda_convex_min=epsilon_res$lambda_convex_min,
                        epsilon=epsilon_res$epsilon,
                        MM=epsilon_res$MM))
  }
}
param_dtf= as.data.frame(param_dtf)


###############################################################################
## Run simulation on a separate dataset to obtain the **optimal regularization
## constants**. At the moment, use the validation dataset. In the future,
## should use a new one.
###############################################################################
require(glmnet)

glmnet_results = lapply(MM_vec, function(MM) {
  ## only use top MM SNPs
  design_matrix = get_design_matrix(MM=MM)
  res = lapply(alpha_vec, function(alpha) {
    glmnet(design_matrix, training_data$status, family="binomial",
                    standardize=FALSE, alpha=alpha)
  })
  names(res) = sapply(alpha_vec, paste)
  return(res)
})
names(glmnet_results) = sapply(MM_vec, paste)




###############################################################################
## Simulations.
###############################################################################

##########################################
## helper functions
##########################################

sim_result = function(MM, epsilon, alpha, lambda, phi, noisy=TRUE) {
  ## aggregate the simulation results
  design_matrix = get_design_matrix(MM)
  nn_var = 1 + ncol(design_matrix)  # add intercept
  if (noisy) {
    noise = rB2(nn_var)
    noise_scale = phi / (epsilon * nrow(design_matrix))
  } else  {
    noise_scale = 0
    noise = rep(0, nn_var)
  }
  XX = cbind(1, design_matrix)
  yy = as.numeric(training_data$status) - 1

  optim_result = run_cvx(XX, yy, lambda, alpha, noise, noise_scale)

  list(optim_result=optim_result,
       params=list(alpha=alpha,
                   lambda=lambda,
                   MM=MM,
                   epsilon=epsilon,
                   phi=phi,
                   noise_scale=noise_scale))
}


get_lambda_vec = function(all_lambda, lambda_convex_min, alpha, nn_lambda) {
  ## first test glmnet's default lambda list
  rescaled_lambda = all_lambda * (1 - alpha) / 2
  valid_lambda_vec = all_lambda[rescaled_lambda > lambda_convex_min]
  while (TRUE) {
    nn_valid_lambda = length(valid_lambda_vec)
    if (nn_valid_lambda > nn_lambda) {
      mult = floor((nn_valid_lambda - 1) / nn_lambda)
      lambda_vec = sort(valid_lambda_vec)[1 + c(0:(nn_lambda - 1)) * mult]
      break
    }
    if (nn_valid_lambda > 0) {
      lambda_vec = c(valid_lambda_vec,
                     max(valid_lambda_vec) * c(2: (1 + nn_lambda - nn_valid_lambda)))
      break
    }
    lambda_vec = lambda_convex_min * c(1:nn_lambda)
    break
  }
  lambda_vec
}


sim_wrapper = function(noisy=TRUE) {
  ## wraper function for the simulations
  res_MM_list = lapply(MM_vec, function(MM) {
    res_epsilon_list = lapply(epsilon_vec, function(epsilon) {
      lambda_convex_min = param_dtf[((param_dtf$MM==MM) & (param_dtf$epsilon==epsilon)), "lambda_convex_min"]
      phi = param_dtf[param_dtf$MM==MM, "phi"][1]
      res_alpha_list = lapply(alpha_vec, function(alpha) {
        ## first test glmnet's default lambda list
        glmnet_fit = glmnet_results[[paste(MM)]][[paste(alpha)]]
        lambda_vec = get_lambda_vec(glmnet_fit$lambda, lambda_convex_min, alpha, nn_lambda)
        res = lapply(lambda_vec, function(lambda) {
          print(sprintf("Processed MM=%s, epsilon=%s, alpha=%s, lambda=%s", MM, epsilon, alpha, lambda))
          sim_result(MM=MM, epsilon=epsilon, alpha=alpha, lambda=lambda, phi=phi, noisy=noisy)
        })
        names(res) = sapply(lambda_vec, paste)
        return(res)
      })
      names(res_alpha_list) = sapply(alpha_vec, paste)
      return(res_alpha_list)
    })
    names(res_epsilon_list) = sapply(epsilon_vec, paste)
    return(res_epsilon_list)
  })
  names(res_MM_list) = sapply(MM_vec, paste)
  return(res_MM_list)
}


##########################################
## get the simulation results when noise is added
## run optimization using cvx in matlab
##########################################


ptm = proc.time()
print(paste("Number of nlm fits = ",
            length(MM_vec) * length(epsilon_vec) * length(alpha_vec) * nn_lambda * nn_sim))

if (TRUE) {
  require(R.matlab)

  Matlab$startServer()
  ## Create a MATLAB client object used to communicate with MATLAB
  matlab <- Matlab()
  ## Connect to the MATLAB server.
  isOpen <- open(matlab)
  ## Confirm that the MATLAB server is open, and running
  if (!isOpen) {
    throw("MATLAB server is not running: waited 30 seconds.")
  }

  source('./run_cvx.R')
  noisy_result = lapply(1:nn_sim, function(ee) {
    print(paste("Simulation iteration", ee))
    return(sim_wrapper(noisy=TRUE))
  })

  ## When done, close the MATLAB client, which will also shutdown
  ## the MATLAB server and the connection to it.
  close(matlab)
}
print(floor((proc.time() - ptm) / 60))

## save the results
save(noisy_result, file="./noisy_result.RData")

##########################################
## get the simulation results when NO noise is added
##########################################
## when we require that the smallest candidate lambda to depend on epsilon

## use glmnet to find the estimates
sim_wrapper.glmnet = function() {
  res_MM_list = lapply(MM_vec, function(MM) {
    design_matrix = get_design_matrix(MM=MM)
    res_epsilon_list = lapply(epsilon_vec, function(epsilon) {
      lambda_convex_min = param_dtf[((param_dtf$MM==MM) & (param_dtf$epsilon==epsilon)), "lambda_convex_min"]
      phi = param_dtf[param_dtf$MM==MM, "phi"][1]
      res_alpha_list = lapply(alpha_vec, function(alpha) {
        ## first test glmnet's default lambda list
        glmnet_fit = glmnet_results[[paste(MM)]][[paste(alpha)]]
        lambda_vec = get_lambda_vec(glmnet_fit$lambda, lambda_convex_min, alpha, nn_lambda)
        res = lapply(lambda_vec, function(lambda) {
          print(sprintf("Processed MM=%s, epsilon=%s, alpha=%s, lambda=%s", MM, epsilon, alpha, lambda))
          fit = glmnet(design_matrix, training_data$status, family="binomial",
                       standardize=FALSE, alpha=alpha, lambda=lambda)
          list(optim_result=list(estimate=c(as.vector(fit$a0),
                                            as.vector(fit$beta))),
               params=list(alpha=alpha,
                           lambda=lambda,
                           MM=MM,
                           epsilon=epsilon,
                           phi='NA',
                           noise_scale=0))
        })
        names(res) = sapply(lambda_vec, paste)
        return(res)
      })
      names(res_alpha_list) = sapply(alpha_vec, paste)
      return(res_alpha_list)
    })
    names(res_epsilon_list) = sapply(epsilon_vec, paste)
    return(res_epsilon_list)
  })
  names(res_MM_list) = sapply(MM_vec, paste)
  return(res_MM_list)
}

##
ptm = proc.time()

print(paste("Number of nlm fits = ",
            length(MM_vec) * length(epsilon_vec) * length(alpha_vec) * nn_lambda ))
no_noise_result = sim_wrapper.glmnet()
print(floor((proc.time() - ptm) / 60))
## save the results
save(no_noise_result, file="./no_noise_result.RData")


