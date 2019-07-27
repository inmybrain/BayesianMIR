##------- Source from BMIR_Fun.v6.0.R: do not edit by hand


#' @name getHypara
#' @title Generate hyperparameters
#' @description Generate hyperparameters for the Monte Carlo Markov Chain sampling.
#' @param tidydata A list created from \code{\link{Tidy_dataset}}.
#' @return A list with components:
#' \item{hp_mu_coef}{For an intercept, the sample mean of labels is given, while \code{0} is given for the others.}
#' \item{hp_a}{\code{=0.1}. Vague informative prior.}
#' \item{hp_b}{\code{=0.1}. Vague informative prior.}
#' \item{hp_g_coef}{\code{1/n} where \code{n} is the sample size. The unit information prior.}
#' \item{hp_pi}{Uniform distribution over a bag. If a bag \code{i} has \code{m_i} instances, then the prior probability is \code{1/m_i}.}
#' \item{hp_inv_Sig_coef}{Identity matrix.}
#' \item{hp_sig2_y}{The sample variance of labels.}
#' @examples
#' bags <- replicate(3, matrix(rnorm(10), 2), simplify = FALSE)
#' tidydata <- Tidy_dataset(label = rnorm(3), feature_inst = bags)
#' getHypara(tidydata)
#' @export
getHypara <- function(tidydata){
  res <- list(
    hp_mu_coef = c(mean(tidydata$label), rep(0, tidydata$nfeature_inst)),
    hp_a = 0.1,
    hp_b = 0.1,
    hp_g_coef = 1 / tidydata$nsample,
    hp_pi = rep(1 / tidydata$ninst, tidydata$ninst),
    hp_inv_Sig_coef = diag(1 + tidydata$nfeature_inst),
    hp_sig2_y = var(tidydata$label)
  )
  res$hp_pi_list <- split(res$hp_pi, tidydata$membership)
  return(res)
}

#' @name getInit
#' @title Generate initial values
#' @description Generate initial values for parameters
#' @param hypara A list created from \code{\link{getHypara}}.
#' @return A list with components:
#' \item{coef}{Random variables from MVN(\code{hypara$hp_mu_coef}, \code{10I}).}
#' \item{sig2_error}{Random variable from Gamma(\code{hypara$hp_a}, \code{hypara$hp_b}).}
#' \item{delta}{Binary variables. Within a bag \code{i}, each random vector is generated from \code{rmultinom(1, hp_pi_list[[i]])}.}
#' @examples
#' bags <- replicate(3, matrix(rnorm(10), 2), simplify = FALSE)
#' tidydata <- Tidy_dataset(label = rnorm(3), feature_inst = bags)
#' getInit(getHypara(tidydata))
#' @export
getInit <- function(hypara){
  list(
    ## coef
    coef = unlist(lapply(hypara$hp_mu_coef, function(mu) rnorm(1, mu, 10))),
    
    ## sig2_error
    sig2_error = rgamma(1, shape = hypara$hp_a, rate = hypara$hp_b),
    
    ## delta
    delta = as.numeric(unlist(lapply(hypara$hp_pi_list, function(x) rmultinom(1, 1, x))))
  )
}

#' @name getPar
#' @title Generate all input parameters
#' @description Generate all input parameters required in \code{\link{Run1Gibbs_cpp}}.
#' @param tidydata A list created from \code{\link{Tidy_dataset}}.
#' @return A list with components:
#' \item{n}{\code{nsample} in \code{\link{Tidy_dataset}}}
#' \item{p}{\code{nfeature} in \code{\link{Tidy_dataset}}}
#' \item{ninst}{\code{ninst} in \code{\link{Tidy_dataset}}}
#' \item{m}{Total number of instances}
#' \item{membership}{{\code{membership} in \code{\link{Tidy_dataset}}}}
#' \item{Y}{{\code{label} in \code{\link{Tidy_dataset}}}}
#' \item{X1}{Covariate matrix of all instances.}
#' \item{hp_sig2_y}{Same as in \code{\link{getHypara}}}
#' \item{hp_mu_coef}{Same as in \code{\link{getHypara}}}
#' \item{hp_a}{Same as in \code{\link{getHypara}}}
#' \item{hp_b}{Same as in \code{\link{getHypara}}}
#' \item{hp_pi}{Same as in \code{\link{getHypara}}}
#' \item{hp_inv_Sig_coef}{Same as in \code{\link{getHypara}}}
#' \item{hp_g_coef}{Same as in \code{\link{getHypara}}}
#' \item{coef}{Same as in \code{\link{getInit}}}
#' \item{sig2_error}{Same as in \code{\link{getInit}}}
#' \item{delta}{Same as in \code{\link{getInit}}}
#' @examples
#' getPar(Tidy_dataset(label = rnorm(3), feature_inst = replicate(3, matrix(rnorm(10), 2), simplify = FALSE)))
#' @export
getPar <- function(tidydata){
  
  list_hypara <- getHypara(tidydata)
  list_init <- getInit(list_hypara)
  
  list(n = tidydata$nsample,
       p = tidydata$nfeature_inst,
       ninst = tidydata$ninst,
       m = sum(tidydata$ninst),
       membership = tidydata$membership,
       Y = tidydata$label,
       X1 = model.matrix(~ 1 +., data = Reduce(rbind, tidydata$feature_inst)),
       
       hp_sig2_y = list_hypara$hp_sig2_y,                    
       hp_mu_coef = list_hypara$hp_mu_coef,
       hp_a = list_hypara$hp_a,
       hp_b = list_hypara$hp_b,
       hp_pi = list_hypara$hp_pi,
       hp_inv_Sig_coef = list_hypara$hp_inv_Sig_coef,
       hp_g_coef = list_hypara$hp_g_coef,
       
       coef = list_init$coef,
       sig2_error = list_init$sig2_error,
       delta = list_init$delta
  )
}


#' @name BMIR_sampler
#' @title The Bayesian Muliple Instance Regression model
#' @description Generate the Monte Carlo Markov Chain samples from the specified Bayesian model.
#' @param ntotal Total number of samplings. Default is \code{100000}.
#' @param nwarm Number of iterations used in burn-in steps. Default is \code{ntotal / 2}.
#' @param nthin The thinning interval. Default is \code{1}.
#' @param nchain The number of the Markov chains. Default is $1$.
#' @param tidydata A list created from \code{\link{Tidy_dataset}}.
#' @param return_delta Logical. Whether to return binary indicators.
#' @return A list with components:
#' \item{mcmclist}{A list of class \code{mcmc} from \code{coda}.}
#' \item{pip}{A matrix with each column having the posterior inclusion probability of each instance, which indicates the probability of being a primary instance in a bag. There are \code{nchain} columns.}
#' @import Rcpp Matrix coda reshape2
#' @examples
#' ## Basic set-up
#' nsample <- 100
#' ninst <- 5
#' nprime <- 1 # do not change
#' nfeature <- 1
#' 
#' ## Generate synthetic data
#' set.seed(6)
#' npoint <- (ninst - nprime) * nfeature
#' bag <- list()
#' for(i in 1:nsample){
#'   prime_inst <- matrix(runif(nfeature, -5, 5), ncol = nfeature)
#'   nonprime_inst <- matrix(c(rnorm(npoint, -5, 1),
#'                             rnorm(npoint, 5, 1)),
#'                           nrow = 2, byrow = TRUE)
#'   nonprime_inst <- nonprime_inst[cbind(sample(1:2, npoint, replace = TRUE), 1:npoint)]
#'   nonprime_inst <- matrix(nonprime_inst, ncol = nfeature, byrow = FALSE)
#'   bag[[i]] <- as.data.frame(rbind(prime_inst, nonprime_inst))
#' }
#' beta_true <- rep(2, nfeature + 1)
#' label <- unlist(lapply(bag, function(feature){
#'   beta_true[1] + as.matrix(feature[1,,drop = FALSE]) %*% beta_true[-1]
#' }))
#' label <- label + rnorm(nsample, mean = 0, sd = 1)
#' 
#' ## Tidy data
#' ### Try MIScatterPlot function to visualize the data
#' tidydata <- Tidy_dataset(label = label,
#'                          feature_inst = bag)
#'                          
#' ## BMIR model fitting
#' BMIR_fit <- BMIR_sampler(ntotal = 1000, tidydata = tidydata)
#' @export
BMIR_sampler <- function(ntotal = 100000, 
                         nwarm = ntotal / 2, 
                         nthin = 1,
                         nchain = 1,
                         tidydata,
                         return_delta = FALSE
){
  parlist <- getPar(tidydata)
  
  cat("=============================================================\n")
  cat(sprintf("Multiple Instance Bayesian Regression\n"))
  res_mcmc <- vector("list", nchain)
  pip <- c()
  for(nc in 1:nchain){
    s_time <- Sys.time()
    parlist <- getPar(tidydata)
    res_mcmc[[nc]] <- BMIR_cpp(ntotal = ntotal,
                               nwarm = nwarm,
                               nthin = nthin,
                               
                               ninst = parlist$ninst,
                               Y = parlist$Y,
                               X1 = parlist$X1,
                               
                               hp_mu_coef = parlist$hp_mu_coef,
                               hp_a = parlist$hp_a,
                               hp_b = parlist$hp_b,
                               hp_g_coef = parlist$hp_g_coef, 
                               hp_sig2_y = parlist$hp_sig2_y,
                               hp_pi = parlist$hp_pi,
                               hp_inv_Sig_coef = parlist$hp_inv_Sig_coef,
                               
                               coef = parlist$coef,
                               sig2_error = parlist$sig2_error,
                               delta = parlist$delta,
                               return_delta = return_delta
    )
    cat(sprintf("Elapsed time for chain%d=%.3f mins: MCMC sampling is done!\n", nc, difftime(Sys.time(), s_time, units = "mins")))
    
    ## save pip
    pip <- cbind(pip, res_mcmc[[nc]]$pip)
    
    ## save MCMC samples for coef, sig2_error
    colnames(res_mcmc[[nc]]$coef) <- paste0("coef", seq_along(parlist$coef))
    colnames(res_mcmc[[nc]]$sig2_error) <- "sig2_error"
    
    if(!return_delta){
      res_mcmc[[nc]]$delta <- NULL
    }
    res_mcmc[[nc]] <- rbind(unlist(parlist[c("coef", "sig2_error")]), # initial
                            cbind(res_mcmc[[nc]]$coef, res_mcmc[[nc]]$sig2_error) # MCMC samples
    )
    rownames(res_mcmc[[nc]]) <- 1:nrow(res_mcmc[[nc]])
    res_mcmc[[nc]] <- mcmc(data = res_mcmc[[nc]], start = 1, end = nrow(res_mcmc[[nc]]), thin = nthin)
  }
  res_mcmc <- mcmc.list(res_mcmc)
  colnames(pip) <- names(res_mcmc) <- paste0("Chain", 1:nchain)
  
  res <- list(
    mcmclist = res_mcmc,
    pip = pip
  )
  return(res)
}


