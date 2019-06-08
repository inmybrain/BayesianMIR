## Set up data for BMIR fitting

#' @name Tidy_dataset
#' @title Tidy a dataset for the multiple instance learning
#' @description Generate a list of components frquently used in the multiple instance learning
#' @param label A vector of response variables (of length r). For unlabeled bags, it is set to NULL.
#' @param feature_inst Either a list of length r or a matrix. If a list is given, each component contains covariate information of multiple instances in a matrix form (of same column size, but possibly different row size). If a matrix is given, \code{membership} that indicate group indices of all instances should be specified.
#' @param membership A vector of integers that indicate group indices of all instances.
#' @return A list with components:
#' \item{label}{Same as in \code{label}.}
#' \item{feature_inst}{Same as in \code{feature_inst}.}
#' \item{nsample}{The number of bags (or samples).}
#' \item{nfeature_inst}{The dimension of covariates.}
#' \item{ninst}{A vector of length r, each entry of which is the number of instances in a bag.}
#' \item{membership}{A vector of integers that indicate group indices of all instances.}
#' @examples
#' ## labeled bags
#' Tidy_dataset(label = rnorm(3), feature_inst = replicate(3, matrix(rnorm(10), 2), simplify = FALSE))
#' 
#' ## unlabeled bags
#' Tidy_dataset(feature_inst = replicate(3, matrix(rnorm(10), 2), simplify = FALSE))
#' @export
Tidy_dataset <- function(label = NULL, feature_inst, membership = NULL){
  ## Transform a matrix or data frame input to a list  
  if(class(feature_inst) == "list"){
    # okay
  } else if(any(class(feature_inst) %in% c("data.frame", "matrix"))){
    if(is.null(membership)){
      stop("\"membership\" should be specified.\n")
    }
    feature_inst <- split(feature_inst, membership)
  } else{
    stop("\"feature_inst\" should be either a list or a matrix (or a data.frame).\n")
  }
  # if(length(label) == 0 & length(feature_inst) == 0){
  #   return(NULL)
  # }
  
  ## Check dimension of feature_inst
  p <- unique(unlist(lapply(feature_inst, ncol)))
  if(length(p) > 1){
    stop("\"feature_inst\" should have an equal number of features.\n")
  }
  
  ## List elements are converted to a data.frame
  if(any(!unlist(lapply(feature_inst, is.data.frame)))){
    feature_inst <- lapply(feature_inst, function(x){
      x <- data.frame(x)
      colnames(x) <- paste0("X", 1:p)
      return(x)
    })
  }
  
  ## Check the number of labels and bags
  if(!is.null(label)){
    if(length(label) != length(feature_inst)){
      stop("Number of bags and labels are not equal.\n")
    }
  }
  
  ## Return a list
  res <- list(label = label,
              feature_inst = feature_inst,
              nsample = length(feature_inst),
              nfeature_inst = p,
              ninst = unlist(lapply(feature_inst, nrow))
  )
  res$membership <- rep(seq_len(res$nsample), res$ninst)
  return(res)
}