##------- Source from BMIR_Fun.v6.0.R: do not edit by hand


#' @name MISummarize
#' @title Summarize Bags with Multiple Instances
#' @description Summarize basic informations about instances in bags
#' @param tidydata A list created from \code{\link{Tidy_dataset}}.
#' @param plot Logical. Draw a scatter plot using \code{\link{MIScatterPlot}}. Default is FALSE.
#' @param bag_size Maximum size of bags to be plotted. This parameter is for reducing the object size.
#' @return No return value.
#' @seealso MIScatterPlot
#' @examples 
#' ## Generate bags
#' tidydata <- Tidy_dataset(label = rnorm(3), feature_inst = replicate(3, matrix(rnorm(10), 2), simplify = FALSE))
#' MISummarize(tidydata)
#' @export
MISummarize <- function(tidydata, plot = FALSE, bag_size = 5){
  # tidydata = tidydata_real
  # plot = T
  # bag_size = 10
  
  cat(sprintf("Number of bags : %d\n", tidydata$nsample))
  cat(sprintf("Number of features : %d(", tidydata$nfeature_inst))
  cat(sprintf("%s,", colnames(tidydata$feature_inst[[1]])[1:min(10, tidydata$nfeature_inst)]))
  cat("...)\n")
  cat(sprintf("Number of instances in bags : \n"))
  print(summary(tidydata$ninst))
  if(plot){
    print(MIScatterPlot(tidydata = tidydata, bag_size = bag_size))
  }
  invisible()
}

#' @name MIScatterPlot
#' @title Scatter Plot for Multiple Instance Regression
#' @description Visualize multiple instances and/or identification of primary instances.
#' @param tidydata A list created from \code{\link{Tidy_dataset}}.
#' @param bag_size Maximum size of bags to be plotted. This parameter is for reducing the object size.
#' @param true_primary A list of logical vectors. Each logical vector denotes which instances are primary in a bag. In general, it is not available in real data.
#' @param pred_primary A list of logical vectors. Each logical vector is a fitted or predicted outcome that indicates which instances are primary in a bag.
#' @return Return a ggplot object created from \code{ggplot}.
#' @import reshape2 ggplot2
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
#' tidydata <- Tidy_dataset(label = label,
#'                          feature_inst = bag)
#'                          
#' ## BMIR model fitting
#' BMIR_fit <- BMIR_sampler(ntotal = 1000, tidydata = tidydata)
#'                          
#' ## Visualize fitted results
#' tp <- lapply(1:tidydata$nsample, function(x) rep(c(TRUE, FALSE), c(1, ninst - 1)))
#' pp <- lapply(split(BMIR_fit$pip[,1], tidydata$membership), function(x) rank(-x, "min") <= 1)
#' MIScatterPlot(tidydata = tidydata, 
#'               bag_size = 5,
#'               true_primary = tp, 
#'               pred_primary = pp)
#' @export
MIScatterPlot <- function(tidydata, bag_size,
                          true_primary = NULL, pred_primary = NULL){
  # tidydata = tidydata_real
  # bag_size = 10
  id_fix_true <- vector("list", tidydata$nsample)
  if(!is.null(true_primary)){
    id_fix_true <- list()
    for(i in 1:tidydata$nsample){
      id_fix_true[[i]] <- which(true_primary[[i]])
    }
    ggdata_true_primary <- Map(function(label, bag, id) {
      suppressWarnings(
        data.frame(label = rep(label, length(id)),
                   bag[id, , drop = FALSE])
      )
    }, tidydata$label, tidydata$feature_inst, id_fix_true)
    ggdata_true_primary <- Reduce(rbind, ggdata_true_primary)
    ggdata_true_primary$primary <- 1
    ggdata_true_primary$data <- 1
  }
  
  id_fix_pred <- vector("list", tidydata$nsample)
  if(!is.null(pred_primary)){
    id_fix_pred <- list()
    for(i in 1:tidydata$nsample){
      id_fix_pred[[i]] <- which(pred_primary[[i]])
    }
    ggdata_pred_primary <- Map(function(label, bag, id) {
      suppressWarnings(
        data.frame(label = rep(label, length(id)),
                   bag[id, , drop = FALSE])
      )
    }, tidydata$label, tidydata$feature_inst, id_fix_pred)
    ggdata_pred_primary <- Reduce(rbind, ggdata_pred_primary)
    ggdata_pred_primary$primary <- 1
    ggdata_pred_primary$data <- 0
  }
  
  id_fix <- Map(c, id_fix_true, id_fix_pred)
  id_fix <- lapply(id_fix, unique)
  
  id_random <- list()
  for(i in 1:tidydata$nsample){
    if(bag_size >= tidydata$ninst[i]){
      id_random[[i]] <- setdiff(1:tidydata$ninst[i], id_fix[[i]])
    } else{
      if(bag_size > length(id_fix[[i]])){
        id_random[[i]] <- sample(setdiff(1:tidydata$ninst[i], id_fix[[i]]), 
                                 size = bag_size - length(id_fix[[i]]),
                                 replace = FALSE)  
      } else{
        id_random[[i]] <- numeric(0)
      }
    }
  }
  ggdata_random <- Map(function(label, bag, id) {
    suppressWarnings(
      data.frame(label = rep(label, length(id)),
                 bag[id, , drop = FALSE])
    )
  }, tidydata$label, tidydata$feature_inst, id_random)
  ggdata_random <- Reduce(rbind, ggdata_random)
  ggdata_random$primary <- 0
  ggdata_random$data <- 1
  
  ggdata <- c()
  for(i in grep("ggdata_", ls())){
    ggdata <- rbind(ggdata, get(ls()[i]))
  }
  ## 0.0 : non-primary, not data (fitted) --> must have no case (by construction)
  ## 0.1 : non-primary, data
  ## 1.0 :     primary, not data (fitted)
  ## 1.1 :     primary, data
  ggdata$intn <- factor(interaction(ggdata$primary, ggdata$data),
                        levels = c("0.1", "1.1", "1.0"),
                        labels = c("Non-prime", "Prime", "Fit"))
  ggdata$data <- NULL 
  ggdata <- melt(ggdata, id.vars = c("label", "intn", "primary")) # wide -> long
  gg_manual_intn <- data.frame(breaks = unique(ggdata$intn),
                               labels = c("Non-prime", "Prime", "Fit")[as.numeric(unique(ggdata$intn))],
                               color_value = gg_color_hue(3)[as.numeric(unique(ggdata$intn))],
                               shape_value = c(4, 4, 2)[as.numeric(unique(ggdata$intn))],
                               stringsAsFactors = FALSE
  )
  gg_manual_intn <- gg_manual_intn[order(gg_manual_intn$breaks),]
  ggbase <- ggplot(data = ggdata, mapping = aes(x = value, 
                                                y = label, 
                                                by = variable, 
                                                color = intn, 
                                                shape = intn, 
                                                alpha = as.factor(primary))) + 
    geom_point() + facet_wrap(facets = ~ variable, ncol = 3, scales = "free_x") + 
    scale_alpha_ordinal(range = c(0.3, 1), guide = FALSE) + 
    scale_colour_manual(values = gg_manual_intn$color_value,
                        labels = gg_manual_intn$labels,
                        breaks = gg_manual_intn$breaks,
                        name = "Type") +
    scale_shape_manual(values = gg_manual_intn$shape_value,
                       labels = gg_manual_intn$labels,
                       breaks = gg_manual_intn$breaks,
                       name = "Type")
  return(ggbase)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

