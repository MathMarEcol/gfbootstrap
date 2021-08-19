#' Generate a bootstrapped GF object
#'
#' The bootstrapped GF object is bootstrapped by fitting many GF models, each with a single tree.
#' So this model is really a list of GF models, with some parameters pre-specified.
#' Many of the options will be disabled if possible.
#' The models are not intended to be combined.
#' The bootstrapped GF model can use predict.gfbootstrap, which generates many points per input point, with a total of n_models * nrow(predict_points)
#'
#' Given a training set of site by species+environment, and a set of relevant GF fitting parameters, you get
#' a single object back with class gfbootstrap, which is internally a list of gf models.
#' Each gf model in the list is fitted with exactly 1 tree and a bootstrapped set of inputs.
#' The ntree parameter is not available, see nbootstraps instead.
#' All other parameters are passed on as is, although some defaults have been changed.
#'
#' Are some parameters meaningless, and worth discarding?
#'
#' gfbootstrap is `future`-aware, so setting up a parallel plan before calling gfbootstrap will speed up processing.
#' Not specifying a `future` plan will default to standard sequential processing.
#' See future.seed
#'
#' @param x data.frame, rows as observations, columns are either predictors or responses
#' @param predictor.vars columns in x that are predictors
#' @param response.vars columns in x that are responses
#' @param nboostraps number of bootstrapped models to fit
#' @param nsamples number of samples taken at each bootstrapped model fit
#'
#' see gradientForest::gradientForest() for mtry, transform, maxLevel, corr.threshold, compact, nbin, trance parameters
#'
#' @export
#'
#' @examples
#'
#' library(gradientForest)
#' data(CoMLsimulation)
#' predictor.vars <- colnames(Xsimulation)
#' response.vars <- colnames(Ysimulation)
#' demobootstrap <- bootstrapGradientForest(x = data.frame(Ysimulation,Xsimulation),
#' predictor.vars = predictor.vars,
#' response.vars = response.vars,
#' nbootstrap = 100, #small, should be 500 or 1000 in a real experiment
#' compact = T, nbin = 200,
#' transform = NULL,
#' corr.threshold = 0.5,
#' maxLevel = floor(log2(length(response.vars)*0.368/2)),
#' trace = TRUE
#' )
bootstrapGradientForest <- function(
                        x,
                        predictor.vars,
                        response.vars,
                        nbootstrap=10,
                        nsamples = nrow(x),
                        mtry=NULL,
                        transform=NULL,
                        maxLevel=0,
                        corr.threshold=0.5,
                        compact=FALSE,
                        nbin=101,
                        max_retries = 10,
                        trace=FALSE
) {
  if(!require(gradientForest)){
    stop("gfbootstrap requires the gradientForest package. Please install it from https://r-forge.r-project.org/R/?group_id=973")
  }

  ##gradientForest is loaded

  gf_bootstrap <- future.apply::future_lapply(1:nbootstrap, function(i,
                        x,
                        predictor.vars,
                        response.vars,
                        nbootstrap,
                        nsamples,
                        mtry,
                        transform,
                        maxLevel,
                        corr.threshold,
                        compact,
                        nbin,
                        max_retries,
                        trace
                                                                     ){

    ##Fit GF with a single tree, but otherwise identical calls.

    ##In bootstrapping, the input data are randomly resampled.
    ##GF internally bootstraps each tree as part of the random forest algorithm
    ##therefore, each loop will produce a bootstrapped GF model

    valid_model <- FALSE
    tries <- 1
    while(!valid_model & tries < max_retries){
      if(tries > 1) message(paste0("GF model failed to fit. Restarting. i: [", i, "], try: ", tries))
      gf_list <-  tryCatch({
        gradientForest::gradientForest(data = x,
                                                predictor.vars=predictor.vars,
                                                response.vars=response.vars,
                                                ntree = 1,
                                                mtry=mtry,
                                                transform=transform,
                                                maxLevel=maxLevel,
                                                corr.threshold=corr.threshold,
                                                compact=compact,
                                                nbin=nbin,
                                                trace=trace
                                                )
      }, error = function(e){
        message(paste0("GF model failed to fit, restarting: ", conditionMessage(e)))
        return(NULL)
    }, warning = function(e){
      message(paste0("GF model failed to fit, restarting: ", conditionMessage(e)))
      return(NULL)
    })

      if(is.null(gf_list)){
        tries <- tries + 1
        valid_model <- FALSE
      } else {
        valid_model <- TRUE
      }

    }
    if(length(unique(gf_list$res$var)) < length(predictor.vars)){
      message("Some vars not found, i: [", i, "], vars: [", unique((gf_list$res$var)),"]")
    }
    return(gf_list)

  }, future.seed = TRUE,
                        x = x,
                        predictor.vars = predictor.vars,
                        response.vars = response.vars,
                        nbootstrap = nbootstrap,
                        nsamples = nsamples,
                        mtry = mtry,
                        transform = transform,
                        maxLevel = maxLevel,
                        corr.threshold = corr.threshold,
                        compact = compact,
                        nbin = nbin,
                        max_retries = max_retries,
                        trace = trace
  )

  if (any(vapply(gf_bootstrap, is.null, logical(1)))){
    stop(paste0("[", sum(vapply(gf_bootstrap, is.null, logical(1))),
                "] GradientForest objects failed to fit even after [", max_retries, "] tries"))
  }

  ##The offsets are applied to the cumimp curves,
  ##but the cumimp curves are calculated
  ##on demand by gradientForest.
  ##HI will just keep the offsets beside the gf models,
  ## and apply the offsets on demand.
  out <- list(
    gf_list = gf_bootstrap,
    offsets = NULL
  )
  class(out) <- "bootstrapGradientForest"
  ##Calculate offsets
    ##the optimal offset requires knowing the distance between curves for every tree in gfbootstrap
    ##
  df_dist <- gfbootstrap:::gfbootstrap_dist(out,
                              x_samples = 100)
  out$offsets <- gfbootstrap:::gfbootstrap_offsets(df_dist)

  return(out)
}

#' gfbootstrap distance matrix
#'
#' Given a list of gradient forest objects,
#' gfbootstrap_dist calculates the area bewteen curves for
#' each pair of gf objects.
#'
#' gfbootstrap_dist() assumes that every GF model has the same
#' set of predictor variables.
#'
#' If an offsets vector is provided, it is used to calculate the result.
#'
#' @param gf_boot a bootstrapGradientForest object
#' @param x_samples integer number of points sampled along each cumulative importance curve, evenly spaced.
#' higher values will give more accurate answers, but consume more time and resources
#'
#' @return Returns a data.frame of areas between curves for each pair and each predictor.
#'
#' @export
gfbootstrap_dist <- function(
                             gf_boot,
                             x_samples = 100){

  if(!require(gradientForest)){
    stop("gfbootstrap_dist requires the gradientForest package. Please install it from https://r-forge.r-project.org/R/?group_id=973")
  }
  gf_list <- gf_boot$gf_list

  assertthat::assert_that(any("gradientForest" %in% class(gf_list[[1]]),
                              "combinedGradientForest" %in% class(gf_list[[1]])))
  pred_vars <- pred_names(gf_boot)
  P <- length(pred_vars)
  K <- length(gf_list)

  if(hasName(gf_boot, "offsets") & !is.null(gf_boot$offsets)){
    offsets <- gf_boot$offsets
    assertthat::assert_that(nrow(offsets) == K)
    assertthat::assert_that(ncol(offsets) == P)
  } else {
    offsets <- matrix(
      rep.int(0, K*P),
      nrow = K,
      ncol = P
    )
    offsets <- as.data.frame(offsets)
    names(offsets) <- pred_vars
  }

  gf_pred_cross <- expand.grid(seq_along(gf_list), pred_vars)
  names(gf_pred_cross) <- c("gf", "pred")
  ##precalculate scores
  ## need the min and max ci for each predictor
  ## using min and max, predict for all gf, using extrap = NA and the same x_samples points
  pred_ci_range <- future.apply::future_lapply(pred_vars, function(pred, gf_list){
    gf_ci_range <- do.call("rbind", future.apply::future_lapply(seq_along(gf_list), function(i, pred) {
      gf <- gf_list[[i]]
              ##calculate overlap
              tryCatch({
                 ret <-gradientForest::cumimp(gf, pred)$x
                if(length(ret) >= 2){
                  return(data.frame(xmin = min(ret), xmax = max(ret)))
                }
                return(data.frame(xmin = NA, xmax = NA))
              },
              error = function(e){
                if(grepl("Predictor [^[:space:]]* does not belong to", e$message)) {
                message(paste0("Tree [", i,
                               "], had no occurences of predictor [", pred,
                               "]"))
                return(data.frame(xmin = NA, xmax = NA))
                } else {
                  stop(e)
                }
              })
    }, pred = pred)
    )
    return(list(xmin = min(gf_ci_range$xmin, na.rm = TRUE), xmax = max(gf_ci_range$xmax, na.rm = TRUE)))
    }, gf_list = gf_list)
  names(pred_ci_range) <- pred_vars

  empty_preds <- vapply(pred_ci_range, function(x){!is.finite(x$xmin)},
    logical(1))

  if(any(empty_preds)){
    stop(paste0("Some predictors were not selected by any bootstap run.\n",
                "Increase bootstrap runs, include more samples, or drop weak predictors.\n",
                "The following predictors were not selected:\n",
                "====\n",
                paste(names(empty_preds)[empty_preds], collapse = "\n"),
                "\n====")
         )
  }

  ## generate predictions for each gf and var
  newdata_df <- do.call("cbind", future.apply::future_lapply(pred_vars, function(pred, pred_ci_range, x_samples){
    pred_range <- pred_ci_range[[pred]]
    points <- seq(pred_range$xmin, pred_range$xmax, length.out = x_samples)
    return(points)
  }, pred_ci_range = pred_ci_range, x_samples = x_samples)
  )
  colnames(newdata_df) <- pred_vars
  newdata_df <- as.data.frame(newdata_df)

  ## Generate predictions
  gf_predictions_list <- future.apply::future_lapply(gf_list, function(gf, newdata_df, pred_vars) {
    ##
    if(class(gf)[1] == "combinedGradientForest"){
      gf_preds <- names(gf$CU)
    } else if(class(gf)[1] == "gradientForest"){
      gf_preds <- as.character(unique(gf$res$var))
    }
    gf_predictions <-  predict(gf, newdata_df[ , gf_preds], extrap = NA)
    missing_preds <- setdiff(pred_vars, gf_preds)
    full_prediction <- as.data.frame(lapply(pred_vars, function(pred, missing_preds, gf_predictions) {
      if(pred %in% missing_preds) {
        return(stats::setNames(list(rep(NA, length.out = nrow(gf_predictions))),
                               nm = pred)
               )
      } else {
        return(stats::setNames(list(gf_predictions[, pred]),
                               nm = pred)
               )
       }
      } , missing_preds = missing_preds, gf_predictions))

    return(full_prediction)
    }, newdata_df = newdata_df, pred_vars = pred_vars)

  ##generate pairs
  d_ij <- expand.grid(i = seq.int(K), j = seq.int(K))

  d_ij_diag <- d_ij[d_ij$i < d_ij$j, ]

  ## For each predictor, generate a distance long-form distance matrix between gf objects
  d_ij_pred <- future.apply::future_lapply(pred_vars, function(pred, gf_predictions_list, d_ij_diag, offsets) {
    dist_est <- apply(d_ij_diag, 1,  function(ind, pred, gf_predictions_list, offsets) {
      ind <- as.data.frame(t(ind))
      i <- ind$i
      j <- ind$j

      gf_i <- gf_predictions_list[[i]][ , pred]
      gf_j <- gf_predictions_list[[j]][ , pred]

      d <-  mean((gf_i + offsets[pred][i,]) - (gf_j + offsets[pred][j,]), na.rm = TRUE)
      }, pred = pred, gf_predictions_list = gf_predictions_list, offsets = offsets)

    d_ij_dist <- data.frame(d_ij_diag, d = dist_est)
    d_ji <- data.frame(d_ij_dist$j,
                       d_ij_dist$i,
                       -d_ij_dist$d)
    names(d_ji) <- names(d_ij_dist)

    d_ij_full <- rbind(d_ij_dist, d_ji)
    return(d_ij_full)
  }, gf_predictions_list = gf_predictions_list, d_ij_diag = d_ij_diag, offsets = offsets)
  names(d_ij_pred) <- pred_vars
  return(d_ij_pred)

}


#' Calculate opimal offsets
#'
#' The distance values returned by `gfbootstrap_dist()`
#' are used to calculate the optimal offsets.
#' Given a predictor, an offset is calculated for each
#' gf cumulative importance curve,
#' such that the total area between the curves is
#' minimised.
#'
#' @param x the data.frame returned by gfbootstrap_dist()
#'
#' @return numeric offsets data.frame, bootstrap sample by predictor
#' 
#' @export
gfbootstrap_offsets <- function(x){

  ## see ../../align_gf_curves.pdf
  ## this is the final form of a method
  ## that minimizes mean squared error between
  ## curves.
  K <- max(x[[1]]$j)
  ##find optimal offsets for each predictor
  mat_A <- K*diag(K) - matrix(1, nrow = K, ncol = K)

  ##Set the offset of the first curve to be 0,
  ##so all offsets are relative to the first curve.
  ##This is neccessary to create a unique solution
  ##because the system has an infinite number of solutions
  ##eg. distances do not change if ALL curves are shifted by +10
  mat_A[1, ] <- 0
  mat_A[1,1] <- 1
  offsets_by <- lapply(x, function(d_pred, mat_A){
                            vec_b <- aggregate(d_pred$d, by = list(alpha_i = d_pred$i), function(x){sum(x, na.rm = TRUE)})
                            vec_b$x[1] <- 0
                            vec_b <- - vec_b
                            solve(mat_A, vec_b$x)
                          },
                          mat_A = mat_A
                    )
  offsets_new <- as.data.frame(lapply(offsets_by, as.vector))
  return(offsets_new)
}

#' plot cumimp for bootStrapGradientForest
#'
#' @param x a bootStrapGradientForest object
#' @param vars names of predictor variables to use
#' @param n_curves Number of curves to plot. If less than the number of bootstrap samples, then samples are chosen randomly
#' @param debug additional plots and messages if TRUE
#'
#' @return Returns a ggplot object, ready for plotting.
#'
#' @export
#'
gg_bootstrapGF <- function(x,
                           vars = names(importance(x$gf_list[[1]], type = "Weighted", sorted = TRUE)),
                           n_curves = 10,
                           debug = TRUE) {
  assertthat::assert_that("bootstrapGradientForest" %in% class(x))
  assertthat::assert_that(is.vector(vars))
  assertthat::assert_that(is.numeric(n_curves))
  assertthat::assert_that(length(n_curves) == 1)

  if(!require(gradientForest)){
    stop("gg_bootstrapGF requires the gradientForest package. Please install it from https://r-forge.r-project.org/R/?group_id=973")
  }

  if (n_curves >= length(x$gf_list)) {
    gf_ind <- seq_along(x$gf_list)
  } else {
    gf_ind <- sample.int(n = length(x$gf_list), size = n_curves)
  }
  ##for each predictor,
  ##get curves from a subset of the gf objects
  ##put into a table and plot
  curve_groups <- expand.grid(vars, gf_ind)
  curve_all <- do.call("rbind", future.apply::future_apply(curve_groups, 1, function(cg, x){
    pred <- as.character(cg[1])
    gf <- as.integer(cg[2])
    
    if (is.element(pred, levels(x$gf_list[[gf]]$res$var)) ){
      curve_data <- gradientForest::cumimp(x$gf_list[[gf]], pred)
    } else {
      message(paste0("gg_bootstrapgf: Tree [", gf,
                     "], had no occurences of predictor [", pred,
                     "]"))
      return(data.frame(pred = NULL, gf = NULL, x = NULL, y = NULL, stringsAsFactors = FALSE))
    }
    curve_data$y <- curve_data$y + x$offsets[gf, pred]
    return(data.frame(pred = pred, gf = gf, x = curve_data$x, y = curve_data$y, stringsAsFactors = FALSE))

  }, x = x)
  )
  if(debug){
    curve_no_offset <- do.call("rbind", future.apply::future_apply(curve_groups, 1, function(cg, x){
      pred <- as.character(cg[1])
      gf <- as.integer(cg[2])
      
      if (is.element(pred, levels(x$gf_list[[gf]]$res$var)) ){
        curve_data <- gradientForest::cumimp(x$gf_list[[gf]], pred)
      } else {
        message(paste0("gg_bootstrapgf: Tree [", gf,
                       "], had no occurences of predictor [", pred,
                       "]"))
        return(data.frame(pred = NULL, gf = NULL, x = NULL, y = NULL, stringsAsFactors = FALSE))
      }

      return(data.frame(pred = pred, gf = gf, x = curve_data$x, y = curve_data$y, stringsAsFactors = FALSE))

    }, x = x)
    )
  }
  if(debug){
    return(list(
      no_offset = ggplot2::ggplot(data = curve_no_offset, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
  ggplot2::facet_wrap(ggplot2::vars(pred), scales = "free_x"),
  no_offset_stretch =  ggplot2::ggplot(data = curve_no_offset, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(pred), scales = "free"),
  offset = ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(pred), scales = "free_x"),
  offset_stretch = ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(pred), scales = "free")
  ))
  }
    return(
      ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(ggplot2::vars(pred), scales = "free_x")
    )
}


#' Predict bootstrapGradientForest
#'
#' Generic S3 function for predicting bootstrapGradientForest
#' objects.
#'
#' The behaviour is very similar to predict.gradientForest,
#' but due to the way bootstrapGradientForest is intended
#' to be used, an extra parameter, `type` is added.
#' @param object a bootstrapGradientForest object
#' @param newdata data.frame of new observations to predict. If NULL, use observations from object.
#' @param type is a vector containing elements from
#' c("mean", "variance", "points", "weight")
#'
#' "mean" gives the unweighted mean of the points (default, for
#' consistency with other predict functions)
#'
#' "diagonal" gives the unweighted variance of the points,
#' but only allows diagonal variance, that is, the variance
#' of each predictor independently.
#'
#' "variance" gives the unweighted variance matrix of the points
#'
#' "points" returns all the points, so for each new x row,
#' k (number of bootstrapped models) rows will be returned
#'
#' "weight" gives the R^2 performance of each bootstrapped model.
#'
#' Due to the fact that many rows of output can be returned
#' for each row of input, x_row is included in the returned
#' data.frame to match inputs to outputs.
#'
#' @param extrap passed on to predict.gradientForest
#' possible valuse are NA, TRUE, FALSE, or a number in the range [0,1]

#' @param ... arguments passed to gradientForest::cumimp()
#' TODO: This code is rather spaghetti, eg. activity 4 depends on 2 and 3,
#' then activity 5 depends on 1 and 4
#'
#' @return A long form data.frame, with $type giving the prediction type.
#'
#' @export
predict.bootstrapGradientForest <- function(object,
                                            newdata,
                                            type = c("mean"),
                                            extrap=TRUE,
                                            ...){

  out <- gfbootstrap:::bootstrap_predict_common(object,
                                            newdata,
                                            type = type,
                                            extrap=extrap,
                                            ...)



  class(out) <- c("list", "predict.bootstrapGradientForest")
  return(out)

}


#' Get predictor names from a gfbootstrap object
#'
#' Gradient Forest will drop variables that are
#' never used in any splits. During bootstrapping,
#' it is likely that at least one model drops a predictor.
#'
#' pred_names() searches all models and finds the
#' full set of predictors, even if some models only
#' use a subset.
#'
#' @param obj gfbootstrap object
#'
#' @return vector of character strings, naming the predictors
#'
#' @export
pred_names <- function(obj) {
  UseMethod("pred_names")
}

#' Get predictor names from a gfbootstrap object
#'
#' Gradient Forest will drop variables that are
#' never used in any splits. During bootstrapping,
#' it is likely that at least one model drops a predictor.
#'
#' pred_names() searches all models and finds the
#' full set of predictors, even if some models only
#' use a subset.
#'
#' @param obj gfbootstrap object
#'
#' @return vector of character strings, naming the predictors
#'
#' @export
pred_names.bootstrapGradientForest <- function(obj) {
unique(do.call("c", future.apply::future_lapply(obj$gf_list,
                                                            function(x){
                                                              levels(x$res$var)
                                                            })
                           ))
}
#             (2) List it as "suggests" in your DESCRIPTION file and
#precede each use with "if(require(desiredR_ForgePackage))".  If
#"require" returns TRUE, you do what you want.  Else issue an error message.


#' Combine bootstrapped gradient forest objects
#'
#' As for gradientForest::combinedGradientForest(),
#' takes two or more bootstrapGradientForest objects
#' and combines them.
#'
#' Each bootStrapGradientForest object has many individual
#' GF objects, and one GF model from each bootstrapGradientForest
#' parameter is taken and used to create one sample in the
#' combined bootstrapped Gradient Forest model.
#'
#' Each bootstrap sample in the combined bootstrapped Gradient Forest
#' model is a combinedGradientForest object, with one GF model from
#' each of the individual bootStrapGradientForest objects.
#'
#' The key difference to combinedGradientForest
#' is that weights must be specified in this call,
#' not at the predict.combinedGradientForest() call, 
#' so offsets for each curve can be calculated now.
#'
#' @param n_samp number of bootstrap samples
#' @param x_samples number of points along each cumimp curve for
#' deciding the best offsets
#' @param nbin number of bins for the cumimp curves
#' @param method see gradientForest::combinedGradientForest
#' @param standardize see gradientForest::combinedGradientForest
#' @param weight see gradientForest::combinedGradientForest
#'
#' @return combinedBootstrapGF object
#'
#' @export
#'
#' @examples
#'
#' library(gradientForest)
#' data(CoMLsimulation)
#' predictor.vars <- colnames(Xsimulation)
#' response.vars <- colnames(Ysimulation)
#' demobootstrap1 <- bootstrapGradientForest(x = data.frame(Ysimulation,Xsimulation),
#' predictor.vars = predictor.vars,
#' response.vars = response.vars[1:6],
#' nbootstrap = 100, #small, should be 500 or 1000 in a real experiment
#' compact = T, nbin = 200,
#' transform = NULL,
#' corr.threshold = 0.5,
#' maxLevel = floor(log2(length(response.vars)*0.368/2)),
#' trace = TRUE
#' )
#' demobootstrap2 <- bootstrapGradientForest(x = data.frame(Ysimulation,Xsimulation),
#' predictor.vars = predictor.vars,
#' response.vars = response.vars[1:6+6],
#' nbootstrap = 100, #small, should be 500 or 1000 in a real experiment
#' compact = T, nbin = 200,
#' transform = NULL,
#' corr.threshold = 0.5,
#' maxLevel = floor(log2(length(response.vars)*0.368/2)),
#' trace = TRUE
#' )
#' democombinedbootstrap <- combinedBootstrapGF(demobootstrap1, demobootstrap2, n_samp = 100)
combinedBootstrapGF <- function(...,
                                n_samp,
                                x_samples = 100,
                                nbin = 101,
                                method = 2,
                                standardize = c("before", "after")[1],
                                weight=c("uniform","species","rsq.total","rsq.mean","site","site.species","site.rsq.total","site.rsq.mean")[3]
                                ) {
  ##... are all the GF objects

  std.options <- c("before","after")
  #TODO: convert to assertthat syntax
  if (is.na(std.option <- pmatch(standardize,std.options)))
    stop(paste('Unmatched standardize value "',standardize,'". Expecting "before" or "after"',sep=""))
  if (is.na(option <- pmatch(weight,c("uniform","species","rsq.total","rsq.mean","site","site.species","site.rsq.total","site.rsq.mean"))))
    stop(paste('Unmatched weight "',weight,'". Expecting one of "uniform", "species", "rsq.total", "rsq.mean", "site", "site.species", "site.rsq.total" or "site.rsq.mean"',sep=""))

  gf_list <- list(...)
  n_gf <- length(gf_list)
  if(!all(sapply(gf_list,inherits,"bootstrapGradientForest")))
    stop("Every argument must be a bootstrapgradientForest")

  if(is.null(gf_names <- names(gf_list)))
    gf_names <- paste("F",1:n_gf,sep="")
  if (any(empty <- gf_names==""))
    gf_names[empty] <- paste("F",1:n_gf,sep="")[empty]

  names(gf_list) <- gf_names

  n_gf_boot <- lapply(gf_list, function(gf){
    length(gf$gf_list)
  })

  ## For each "trial",
  ## randomly select a gf model from each gfbootstrap
  ## and create a combinedGradientForest
  ## from the selected gf models
  combin <- data.frame(lapply(n_gf_boot, function(n, n_samp){
    sample.int(n, n_samp, replace = TRUE)
  }, n_samp = n_samp) )

  gf_combine <- future.apply::future_apply(combin, 1, function(samp, gf_list, nbin, method, standardize){
    gf_samp <- lapply(names(gf_list), function(gf, gf_list, samp){
      gf_list[[gf]]$gf_list[[samp[gf]]]
    }, gf_list = gf_list, samp = samp)
    gf_samp$nbin <- nbin
    gf_samp$method <- method
    gf_samp$standardize <- standardize
    gf_combined <- do.call(gradientForest::combinedGradientForest, gf_samp)
    return(gf_combined)
  }, gf_list = gf_list, nbin = nbin, method = method, standardize = standardize)

  ##calculate the new offsets
  gf_combine <- list(gf_list = gf_combine)
  gf_combine$weight <- weight
  class(gf_combine) <- c("combinedBootstrapGF", "list")

  gf_combine_dist <- gfbootstrap::gfbootstrap_dist(gf_combine, x_samples = x_samples)

  gf_combine$offsets <- gfbootstrap::gfbootstrap_offsets(gf_combine_dist)

  return(gf_combine)
}




#' Get predictor names from a combined gfbootstrap object
#'
#' Gradient Forest will drop variables that are
#' never used in any splits. During bootstrapping,
#' it is likely that at least one model drops a predictor.
#'
#' pred_names() searches all models and finds the
#' full set of predictors, even if some models only
#' use a subset.
#'
#' @param obj combinedBootstrapGF object
#'
#' @return vector of character strings, naming the predictors
#'
#' @export
pred_names.combinedBootstrapGF <- function(obj) {
  unique(do.call("c", future.apply::future_lapply(obj$gf_list,
                                                  function(x){
                                                    ##now each is a combined GF
                                                    names(x$CU)
                                                  })
                 ))
}

#' plot cumimp for bootStrapGradientForest
#'
#' @param x a combinedBootstrapGradientForest object
#' @param vars names of predictor variables to use
#' @param n_curves Number of curves to plot. If less than the number of bootstrap samples, then samples are chosen randomly
#' @param debug additional plots and messages if TRUE
#'
#' @return Returns a ggplot object, ready for plotting.
#'
#' @export
gg_combined_bootstrapGF <- function(x,
                           vars = gfbootstrap::pred_names(x),#names(importance(x$gf_list[[1]], type = "Weighted", sorted = TRUE)),
                           n_curves = 10,
                           debug = TRUE) {
  assertthat::assert_that(inherits(x, "combinedBootstrapGF"))
  assertthat::assert_that(is.vector(vars))
  assertthat::assert_that(is.numeric(n_curves))
  assertthat::assert_that(length(n_curves) == 1)

  if (n_curves >= length(x$gf_list)) {
    gf_ind <- seq_along(x$gf_list)
  } else {
    gf_ind <- sample.int(n = length(x$gf_list), size = n_curves)
  }
  ##for each predictor,
  ##get curves from a subset of the gf objects
  ##put into a table and plot
  curve_groups <- expand.grid(vars, gf_ind)
  curve_all <- do.call("rbind", future.apply::future_apply(curve_groups, 1, function(cg, x){
    pred <- as.character(cg[1])
    gf <- as.integer(cg[2])
    
    if (is.element(pred, names(x$gf_list[[gf]]$CU)) ){
      curve_data <- gradientForest::cumimp(x$gf_list[[gf]], pred, weight = x$weight)
    } else {
      message(paste0("gg_combined_bootstrapgf: Tree [", gf,
                     "], had no occurences of predictor [", pred,
                     "]"))
      return(data.frame(pred = NULL, gf = NULL, x = NULL, y = NULL, stringsAsFactors = FALSE))
    }
    curve_data$y <- curve_data$y + x$offsets[gf, pred]
    return(data.frame(pred = pred, gf = gf, x = curve_data$x, y = curve_data$y, stringsAsFactors = FALSE))

  }, x = x)
  )
  if(debug){
    curve_no_offset <- do.call("rbind", future.apply::future_apply(curve_groups, 1, function(cg, x){
      pred <- as.character(cg[1])
      gf <- as.integer(cg[2])
      
      if (is.element(pred, names(x$gf_list[[gf]]$CU)) ){
        curve_data <- gradientForest::cumimp(x$gf_list[[gf]], pred, weight = x$weight)
      } else {
        message(paste0("gg_combined_bootstrapgf: Tree [", gf,
                       "], had no occurences of predictor [", pred,
                       "]"))
        return(data.frame(pred = NULL, gf = NULL, x = NULL, y = NULL, stringsAsFactors = FALSE))
      }

      return(data.frame(pred = pred, gf = gf, x = curve_data$x, y = curve_data$y, stringsAsFactors = FALSE))

    }, x = x)
    )
  }
  if(debug){
    return(list(
      no_offset = ggplot2::ggplot(data = curve_no_offset, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
  ggplot2::facet_wrap(ggplot2::vars(pred), scales = "free_x"),
  no_offset_stretch =  ggplot2::ggplot(data = curve_no_offset, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(pred), scales = "free"),
  offset = ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(pred), scales = "free_x"),
  offset_stretch = ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(pred), scales = "free")
  ))
  }
    return(
      ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(ggplot2::vars(pred), scales = "free_x")
    )
}

#' Predict combinedBootstrapGF
#'
#' Generic S3 function for predicting combinedBootstrapGF
#' objects.
#'
#' The behaviour is very similar to predict.gradientForest,
#' but due to the way bootstrapGradientForest is intended
#' to be used, an extra parameter, `type` is added.
#'
#' @param object a combinedBootstrapGF object
#' @param newdata data.frame of new observations to predict. If NULL, use observations
#' from the first GF object.
#' @param type is a vector containing elements from
#' c("mean", "variance", "points", "weight")
#'
#' "mean" gives the unweighted mean of the points (default, for
#' consistency with other predict functions)
#'
#' "diagonal" gives the unweighted variance of the points,
#' but only allows diagonal variance, that is, the variance
#' of each predictor independently.
#'
#' "variance" gives the unweighted variance matrix of the points
#'
#' "points" returns all the points, so for each new x row,
#' k (number of bootstrapped models) rows will be returned
#'
#' "weight" gives the R^2 performance of each bootstrapped model.
#'
#' Due to the fact that many rows of output can be returned
#' for each row of input, x_row is included in the returned
#' data.frame to match inputs to outputs.
#'
#' @param extrap passed on to predict.gradientForest
#' possible valuse are NA, TRUE, FALSE, or a number in the range [0,1]
#'
#' @param ... arguments passed to gradientForest::cumimp()
#'
#' @return A long form data.frame, with $type giving the prediction type.
#'
#' @export
#'
predict.combinedBootstrapGF <- function(object,
                                            newdata,
                                            type = c("mean"),
                                            extrap=TRUE,
                                            ...){

  out <- gfbootstrap:::bootstrap_predict_common(object,
                                            newdata,
                                            type = type,
                                            extrap=extrap,
                                            ...)



  class(out) <- c("list", "predict.combinedBootstrapGF")
  return(out)

}

#' Internal function
#'
#' predicting bootstrapGradientForest and combinedBootstrapGF
#' have a lot of overlap.
#'
#' This function pulls together common operations.
bootstrap_predict_common <- function(object,
                                            newdata,
                                            type = c("mean"),
                                            extrap=TRUE,
                                            ...) {

  assertthat::assert_that(length(extrap) == 1)

  if (missing(newdata)){
    newdata <- object$gf_list[[1]]$X[,gfbootstrap::pred_names(object)]
  }
  assertthat::assert_that(inherits(newdata,"data.frame"))

  newnames <- names(newdata)
  pred_vars <- gfbootstrap::pred_names(object)
  if(!all(ok <- newnames %in% pred_vars)) {
    badnames <- paste(newnames[!ok], collapse=", ")
    stop(paste("the following predictors are not in any of the bootstrapGradientForests:\n\t",badnames,sep=""))
  }

  ## get all predictions
  gf_predictions_list <- future.apply::future_lapply(object$gf_list, function(gf, newdata, pred_vars, extrap) {
    if(class(gf)[1] == "combinedGradientForest"){
      gf_preds <- names(gf$CU)
    } else if(class(gf)[1] == "gradientForest"){
      gf_preds <- as.character(unique(gf$res$var))
    }
    gf_predictions <-  predict(gf, newdata[ , gf_preds], extrap = extrap)
    missing_preds <- setdiff(pred_vars, gf_preds)
    full_prediction <- as.data.frame(lapply(pred_vars, function(pred, missing_preds, gf_predictions) {
      if(pred %in% missing_preds) {
        return(stats::setNames(list(rep(NA, length.out = nrow(gf_predictions))),
                               nm = pred)
               )
      } else {
        return(stats::setNames(list(gf_predictions[, pred]),
                               nm = pred)
               )
       }
      } , missing_preds = missing_preds, gf_predictions))

    return(full_prediction)
    }, newdata = newdata, pred_vars = pred_vars, extrap = extrap)

  newdata_long <- stats::reshape(
                           newdata,
                           idvar = "x_row",
                           varying = names(newdata),
                           times = names(newdata),
                           v.names = "x",
                           direction = "long",
                           timevar = "pred")
  gf_predictions_long <- do.call("rbind",
                                 future.apply::future_lapply(
                                                 seq_along(gf_predictions_list),
                                                 function(i, gf_predictions_list, object, newdata_long) {
                                                   gf_pred <- gf_predictions_list[[i]]
                                                   gf_pred <- gf_pred + object$offsets[rep(i, length.out = nrow(gf_pred)), ]
                                                   gf_pred_long <- stats::reshape(
                                                                         gf_pred,
                                                                         idvar = "x_row",
                                                                         varying = names(gf_pred),
                                                                         times = names(gf_pred),
                                                                         v.names = "y",
                                                                         direction = "long",
                                                                         timevar = "pred")
                                                   gf_pred_full <- merge(gf_pred_long, newdata_long, by = c("pred", "x_row"))
                                                   gf_pred_full <- data.frame(gf = i, gf_pred_full)
                                                   return(gf_pred_full)
                                                   }, gf_predictions_list = gf_predictions_list, object = object, newdata_long = newdata_long)
                                 )

  ##Now generate summary results
  #out <- data.frame(type = NA, pred = NA, x_row = NA, x = NA, y = NA, gf_model = NA)[numeric(0), ]
  out <- list()
  all_opts <- c("mean", "variance", "points", "weight")
  if ("mean" %in% type){
    out$mean <- do.call("rbind", future.apply::future_by(gf_predictions_long,
                            list(pred = gf_predictions_long$pred,
                                 x_row = gf_predictions_long$x_row),
                            function(x) {
                              data.frame(pred = unique(x$pred), x_row = unique(x$x_row), x = unique(x$x),
                                         y = mean(x$y))
                            }

                            ))

  }
  if ("diagonal" %in% type){
    out$diagonal <- do.call("rbind", future.apply::future_by(gf_predictions_long,
                            list(pred = gf_predictions_long$pred,
                                 x_row = gf_predictions_long$x_row),
                            function(x) {
                              data.frame(pred = unique(x$pred), x_row  = unique(x$x_row), x = unique(x$x),
                                         y = var(x$y))
                            }
      ))

  }
  if ("variance" %in% type){
    out$variance <- future.apply::future_by(gf_predictions_long,
                            list(x_row = gf_predictions_long$x_row),
                            function(x, pred_vars) {

                              site_wide <- stats::reshape(
                                                    x,
                                                    idvar = "gf",
                                                    timevar = "pred",
                                                    v.names = "y",
                                                    direction = "wide",
                                                    drop = c("x_row", "x"),
                              )
                              names(site_wide) <- sub("^y.", "", names(site_wide))
                              cov(site_wide[, pred_vars])
                            }, pred_vars = pred_vars, simplify = FALSE
      )

  }
  if ("points" %in% type){
    out$points <- gf_predictions_long
  }

  return(out)
}

#' Get variable importance from a bootstrap Gradient Forest
#'
#' Finds the overall importance of each predictor by averaging
#' the importance of each predictor from each bootstrap run.
#' Predictors with more importance are better predictors of the
#' response variables.
#'
#' returns a named vector of predictor importances, possibly sorted.
#' @param x A bootstrapGradientForest object
#' @param type What kind of importance to use. See importance.gradientForest in the gradientForest package for more details.
#' @param sort Should the predictors be sorted in order of importance?
#' @export
#' @importFrom extendedForest importance
importance.bootstrapGradientForest <- function(x,
                                           type = c("Accuracy","Impurity","Weighted","Raw","Species")[3],
                                           sort = TRUE) {

    if (!inherits(x,"bootstrapGradientForest"))
      stop(paste("'x' must be a bootstrapGradientForest object"))

  return(importance_bootstrap_common(x, type, sort))
}
#
#' Get variable importance from a combined bootstrapped Gradient Forest
#'
#' Finds the overall importance of each predictor by averaging
#' the importance of each predictor from each bootstrap run.
#' Predictors with more importance are better predictors of the
#' response variables.
#'
#' returns a named vector of predictor importances, possibly sorted.
#' @param x A combinedBootstrapGF object
#' @param type What kind of importance to use. See importance.combinedGradientForest in the gradientForest package for more details.
#' @param sort Should the predictors be sorted in order of importance?
#' @export
#' @importFrom extendedForest importance
importance.combinedBootstrapGF <- function(x,
                                           type = c("Weighted","Raw","Species")[1],
                                           sort = TRUE) {

    if (!inherits(x,"combinedBootstrapGF"))
      stop(paste("'x' must be a combinedBootstrapGF object"))

  return(importance_bootstrap_common(x, type, sort))
}

#' Internal function
#'
#' calculating importance for bootstrapGradientForest and
#' combinedBootstrapGF has a lot of overlap.
#'
#' This function pulls together common operations.
#' @importFrom extendedForest importance
importance_bootstrap_common <- function(x,
                                        type,
                                        sort
                                        ) {
  preds <- pred_names(x)
  ## Basic logic: Take average importance over all bootstrap samples
  ## Deal with incomplete names
  imp_list <- lapply(x$gf_list, function(gf, preds, type) {
    imp <- importance(gf, type = type, sort = FALSE)
    out <- imp[preds]
    return(out)
  }, preds = preds, type = type)

  ## Assumes that GF assigns importance of 0 to predictors
  ## that were supplied but not used, and an importance of NA
  ## to predictors that were not supplied
  imp_combined <- colMeans(do.call(rbind, imp_list), na.rm=TRUE)
  if (sort) {
    out <- sort(imp_combined,decreasing=TRUE)
  }else {
    out <- imp_combined
  }
  return(out)

}
