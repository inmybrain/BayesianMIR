##------- Source from BMIR_Fun.v6.0.R: do not edit by hand


#' @name predict.BMIR
#' @title Prediction for the Bayesian Muliple Instance Regression Model
#' @description Predict primary instances and labels in new samples (bags).
#' @param BMIRchain One of the Markov chains \code{mcmclist} from \code{\link{BMIR_sampler}}.
#' @param pip One of columns of \code{pip} from \code{\link{BMIR_sampler}}.
#' @param tidydata A list created from \code{\link{Tidy_dataset}} used in model training.
#' @param newtidydata A list created from \code{\link{Tidy_dataset}} only using \code{feature_inst}.
#' @param k Numeric. The number of primary instances in each bag. Default is \code{1}.
#' @return A list with components:
#' \item{newtidydata}{\code{newtidydata} appended with predicted labels (in \code{label}), predicted posterior inclusion probability (in \code{predpip}), and  predicted binary indicators (in \code{predind})}.
#' @import randomForest
#' @examples
#' ## Basic set-up
#' nsample <- 200 # 100 for training, 100 for test
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
#' tidydata <- Tidy_dataset(label = label[1:100],
#'                          feature_inst = bag[1:100])
#' newtidydata <- Tidy_dataset(feature_inst = bag[-(1:100)])
#' 
#' ## BMIR model fitting
#' BMIR_fit <- BMIR_sampler(ntotal = 1000, tidydata = tidydata)
#' 
#' ## Prediction in new bags
#' pred_fit <- predict.BMIR(BMIRchain = BMIR_fit$mcmclist$Chain1, 
#'                          pip = BMIR_fit$pip[,1], 
#'                          tidydata = tidydata, 
#'                          newtidydata = newtidydata, 
#'                          k = 1)
#' ggplot2::ggplot(data = data.frame(pred = pred_fit$newtidydata$label, 
#'                          true = label[-(1:100)]), 
#'        mapping = ggplot2::aes(x = pred, y = true)) + 
#'   ggplot2::geom_point() + ggplot2::geom_abline(intercept = 0, slope = 1, color = "red")
#' @export
predict.BMIR <- function(BMIRchain, pip, tidydata, newtidydata, k = 1){
  # BMIRchain = BMIR_fit$mcmclist$Chain1
  # pip = BMIR_fit$pip[,1]
  # tidydata = tidydata
  # newtidydata = newtidydata
  # k = 1
  
  if(k <= 0){
    stop("k should be strictly bigger than 0.\n")
  }
  if(any(k > newtidydata$ninst)){
    warning("k is larger than bag size of some bags.\nSuch bags will use all instances in prediction.\n")
    # stop("Please choose smaller k (<= bag sizes).\n") 
  }
  if(is.list(pip)){
    pip <- unlist(pip)
  }
  traindata <- data.frame(Reduce(rbind, tidydata$feature_inst),
                          membership = tidydata$membership,
                          pip = pip
  )
  traindata$pip <- logit(traindata$pip - rep(1 / tidydata$ninst, tidydata$ninst))
  colnames(traindata) <- c(paste0("X", 1:tidydata$nfeature_inst), "membership", "pip")
  newdata <- data.frame(X = Reduce(rbind, newtidydata$feature_inst),
                        membership = newtidydata$membership,
                        ninst = rep(newtidydata$ninst, newtidydata$ninst)
  )
  colnames(newdata) <- c(paste0("X", 1:newtidydata$nfeature_inst), "membership", "ninst")
  
  ## Prediction of primary instances in new bags
  fit <- randomForest(x = traindata[,paste0("X", 1:tidydata$nfeature_inst), drop = FALSE], 
                      y = traindata$pip, 
                      xtest = newdata[,paste0("X", 1:newtidydata$nfeature_inst), drop = FALSE])
  # tidydata$fittedpip <- split(fit$predicted, newdata$membership)
  # tidydata$fittedind <- lapply(tidydata$fittedpip, function(x) rank(x, ties.method = "min") <= k)
  
  newtidydata$predpip <- split(exp(fit$test$predicted) / (1 + exp(fit$test$predicted)) + 1 / newdata$ninst,
                               newdata$membership) # logit back-transform
  if(0 < k & k < 1){
    newtidydata$predind <- lapply(newtidydata$predpip, function(x) rank(-x, ties.method = "min") <= max(length(x) * k, 1)) # rank in decreasing order
  } else{
    newtidydata$predind <- lapply(newtidydata$predpip, function(x) rank(-x, ties.method = "min") <= k) # rank in decreasing order
  }
  
  
  ## Bag-level prediction
  coef <- colMeans(BMIRchain)[paste0("coef", 1:(tidydata$nfeature + 1))]
  newtidydata$label <- unlist(Map(function(feature, id, pip){
    ypred <- coef[1] + as.matrix(feature[id,,drop = FALSE]) %*% coef[-1]
    return(sum(ypred * pip[id] / sum(pip[id])))
  }, newtidydata$feature_inst, newtidydata$predind, newtidydata$predpip))
  # newtidydata$label <- unlist(Map(function(feature, id, pip){
  #   ypred <- coef[1] + as.matrix(feature[id,,drop = FALSE]) %*% coef[-1]
  #   return(sum(ypred * pip[id]))
  # }, newtidydata$feature_inst, newtidydata$predind, newtidydata$predpip))
  
  res <- list(
    # tidydata = tidydata,
    newtidydata = newtidydata
  )
  return(res)
}

logit <- function(x, tol = 1e-8){
  x[x < tol] <- tol
  x[x > 1 - tol] <- 1 - tol
  return(log(x / (1 - x)))
}
