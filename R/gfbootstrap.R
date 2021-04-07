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
                        trace=FALSE
) {
  if(!require(gradientForest)){
    stop("gfbootstrap requires the gradientForest package. Please install it from https://r-forge.r-project.org/R/?group_id=973")
  }

  ##gradientForest is loaded

  gf_bootstrap <- future.apply::future_lapply(X = 1:nbootstrap, function(i){

    ##Fit GF with a single tree, but otherwise identical calls.

    ##In bootstrapping, the input data are randomly resampled.
    ##GF internally bootstraps each tree as part of the random forest algorithm
    ##therefore, each loop will produce a bootstrapped GF model

    valid_model <- FALSE
    tries <- 1
    while(!valid_model & tries < 10){
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

  }, future.seed = TRUE)

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
  df_dist <- gfbootstrap:::gfbootstrap_dist_fast(out,
                              x_samples = 100)
  out$offsets <- gfbootstrap:::gfbootstrap_offsets_fast(df_dist)

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
#' @param gf_list a bootstrapped list of Gradient Forest objects, created by bootStrapGradientForest()
#' @param offsets data.frame of sample by predictors, giving the offset for each curve. The default, NULL, sets offsets to 0.
#' @param x_samples integer number of points sampled along each cumulative importance curve, evenly spaced.#'
#' @return Returns a data.frame of areas between curves for each pair and each predictor.
#'
#' @export
gfbootstrap_dist <- function(gf_boot,
                             x_samples = 100){

  if(!require(gradientForest)){
    stop("gfbootstrap_dist requires the gradientForest package. Please install it from https://r-forge.r-project.org/R/?group_id=973")
  }

  gf_list <- gf_boot$gf_list

  assertthat::assert_that("gradientForest" %in% class(gf_list[[1]]))
  pred_vars <- pred_names(gf_boot)
  P <- length(pred_vars)
  K <- length(gf_list)

  if(!is.null(gf_boot$offsets)){
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

  ##generate pairs
  d_ij <- expand.grid(i = seq.int(K), j = seq.int(K))

  d_ij_diag <- d_ij[d_ij$i < d_ij$j, ]

  ##for each pair, find area between curves
  d_ij_pred <- do.call("rbind", future.apply::future_apply(d_ij_diag, 1, function(ind){
    ind <- as.data.frame(t(ind))
    i <- ind$i
    j <- ind$j
    assertthat::assert_that( "gradientForest" %in% class(gf_list[[i]]))
    assertthat::assert_that( "gradientForest" %in% class(gf_list[[j]]))
    tmp <- do.call("rbind",
            future.apply::future_lapply(pred_vars, function(pred, gf_list, offsets){
              ##calculate overlap
              if (is.element(pred, levels(gf_list[[i]]$res$var)) ){
                x_cumimp_i <- gradientForest::cumimp(gf_list[[i]], pred)$x
              } else {
                message(paste0("Tree [", i,
                               "], had no occurences of predictor [", pred,
                               "]"))
                return(data.frame(i, j, pred,  d=0))
              }
              if (is.element(pred, levels(gf_list[[j]]$res$var)) ){
                x_cumimp_j <- gradientForest::cumimp(gf_list[[j]], pred)$x
              } else {
                message(paste0("Tree [", j,
                               "], had no occurences of predictor [", pred,
                               "]"))
                return(data.frame(i, j, pred,  d=0))
              }

              ## message(names(gf_list[[i]]$X))
              ## message(names(gf_list[[j]]$X))
              x_range <- c(max(
                  min(x_cumimp_i),
                  min(x_cumimp_j)
                ),
                min(
                  max(x_cumimp_i),
                  max(x_cumimp_j)
                )
                )
              if(x_range[2] <= x_range[1]){

                message(paste0("Tree [", i,
                               "] and Tree [", j,
                               "], for predictor [", pred,
                               "] had no overlap in cumimp curves"))
                return(data.frame(i, j, pred,  d=0))

              }
              ##predict x_samples points within overlapping ranges
              x_sample_points <- data.frame(seq(x_range[1], x_range[2],
                                                       length.out = x_samples))
              names(x_sample_points) <- pred
              ## message("i")
              x_i <- gradientForest::predict.gradientForest(gf_list[[i]],
                                                            x_sample_points,
                                                            extrap = FALSE)
              ## message("i passed. j")
              x_j <- gradientForest::predict.gradientForest(gf_list[[j]],
                                                            x_sample_points,
                                                            extrap = FALSE)
              ## message("j passed")
              d <-  sum((x_i + offsets[pred][i,]) - (x_j + offsets[pred][j,]))
              d <- d/x_samples
              return(data.frame(i, j, pred, d))
            }, gf_list = gf_list, offsets = offsets)
            )
    return(tmp)
  })
  )

  ## fill in the lower triange, d_ji
  ##take d_ij, and negate $d, swap i, j
  d_ji <- data.frame(d_ij_pred$j,
                     d_ij_pred$i,
                     d_ij_pred$pred,
                     -d_ij_pred$d)
  names(d_ji) <- names(d_ij_pred)

  d_ij_full <- rbind(d_ij_pred, d_ji)
  return(d_ij_full)
}

gfbootstrap_dist_fast <- function(
                                  gf_boot,
                             offsets = NULL,
                             x_samples = 100){

  if(!require(gradientForest)){
    stop("gfbootstrap_dist requires the gradientForest package. Please install it from https://r-forge.r-project.org/R/?group_id=973")
  }
  gf_list <- gf_boot$gf_list

  assertthat::assert_that("gradientForest" %in% class(gf_list[[1]]))
  pred_vars <- pred_names(gf_boot)
  P <- length(pred_vars)
  K <- length(gf_list)

  if(!is.null(gf_boot$offsets)){
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
    gf_ci_range <- do.call("rbind", future.apply::future_lapply(gf_list, function(gf, pred) {
              ##calculate overlap
              if (is.element(pred, levels(gf$res$var)) ){
                x_cumimp_i <- gradientForest::cumimp(gf, pred)$x
                if(length(x_cumimp_i) >= 2){
                  return(data.frame(xmin = min(x_cumimp_i), xmax = max(x_cumimp_i)))
                }
                return(data.frame(xmin = NA, xmax = NA))
              } else {
                message(paste0("Tree [", i,
                               "], had no occurences of predictor [", pred,
                               "]"))
                return(data.frame(xmin = NA, xmax = NA))
              }
    }, pred = pred)
    )
    return(list(xmin = min(gf_ci_range$xmin, na.rm = TRUE), xmax = max(gf_ci_range$xmax, na.rm = TRUE)))
    }, gf_list = gf_list)
  names(pred_ci_range) <- pred_vars

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
    gf_preds <- as.character(unique(gf$res$var))
    gf_predictions <-  predict.gradientForest(gf, newdata_df[ , gf_preds], extrap = NA)
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

  K <- max(x$j)
  ##find optimal offsets for each predictor
  mat_A <- K*diag(K) - matrix(1, nrow = K, ncol = K)

  ##Set the offset of the first curve to be 0,
  ##so all offsets are relative to the first curve.
  ##This is neccessary to create a unique solution
  ##because the system has an infinite number of solutions
  ##eg. distances do not change if ALL curves are shifted by +10
  mat_A[1, ] <- 0
  mat_A[1,1] <- 1
  offsets_by <- future.apply::future_by(x, list(pred = x$pred),
                          function(d_pred, mat_A){
                            vec_b <- aggregate(d_pred$d, by = list(alpha_i = d_pred$i), function(x){sum(x)})
                            vec_b$x[1] <- 0
                            vec_b <- - vec_b
                            solve(mat_A, vec_b$x)
                          },
                          mat_A = mat_A
                    )

  offsets_new <- as.data.frame(lapply(offsets_by, as.vector))
  return(offsets_new)
}

gfbootstrap_offsets_fast <- function(x){

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
  offsets_by <- future.apply::future_lapply(x, function(d_pred, mat_A){
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
  ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
  no_offset_stretch =  ggplot2::ggplot(data = curve_no_offset, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(var), scales = "free"),
  offset = ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
  offset_stretch = ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(var), scales = "free")
  ))
  }
    return(
      ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x")
    )
}


#' Compression Function
#'
#' Copied from rphildyerphd package

#' Helper function to apply power compression to
#' GF predictions.
#'
#' The compressed value lies between a linear extrapolation (ext)
#' and capping at the maximum known value (cap)
#'
#' @description
#'
#' The logic is:
#'
#' We have two lines:
#' \deqn{y_{cap} = b}{y_cap = b}
#' \deqn{y_{ext} = ax + b}{y_ext = ax + b}
#'
#' where we set \eqn{x = 0} to be the point where extrapolation begins, so \eqn{y_{cap} = y_{ext}} at \eqn{x = 0}, and
#' \eqn{b} is the capped value.
#'
#' We want to define a third line, \eqn{z(x)}, st.
#'
#' \deqn{y_{cap} \le z(x) \le y_{ext} \all x \ge 0}{y_cap <= z(x) <= y_ext AA x >= 0}
#'
#' Two clear choices are sigmoid functions that asymptote, and power functions that do not asymptote.
#'
#' Here I will apply a power function to allow compositional turnover to grow indefinitely.
#'
#' Let \eqn{0 \le p \le 1}{0 <= p <= 1} be the power x is raised to. Then:
#'
#' \deqn{z = (x + c)^p + d}
#'
#' where \eqn{c,d} are constants.
#'
#' \eqn{x+c} for \eqn{p} between 0 and 1 grows faster than linear when \eqn{x+c} is close to 0, but it is also convex and monotonically increasing.
#' Therefore choosing \eqn{c,d} st. \eqn{y_{cap}(x)}{y_cap} is tangent to \eqn{z(x)} at \eqn{x = 0} will satisfy
#' \eqn{y_{cap} \le z(x) \le y_{ext} \all x \ge 0}{y_cap <= z(x) <= y_ext AA x >= 0}.
#'
#' At the tangent we know:
#'
#' \deqn{x = 0}
#' \deqn{y_{ext}(0) = z(0)}{y_ext(0) = z(0)}
#' \deqn{y_{ext}'(0) = z'(0)}{y'_ext(0) = z'(0)}
#'
#' \deqn{y_{ext}'(0) = z'(0)}{y'_ext(0) = z'(0)}
#' \deqn{a = p(x+c)^{p-1}}{a = p(x+c)^(p-1)}
#' \deqn{a = p(c)^{p-1} by x = 0}{a = p(c)^(p-1) by x = 0}
#' \deqn{c = \frac{a}{p}^\frac{1}{p-1}}{c = (a/p)^(1/(p-1))}
#'
#' \deqn{y_{ext}(0) = z(0)}{y_ext(0) = z(0)}
#' \deqn{b = c^p + d}{b = c^p + d}
#' \deqn{d = b - c^p}{d = b - c^p}
#' \deqn{d = b - \frac{a}{p}^\frac{p}{p-1}}{d = b - (a/p)^(p/(p-1))}
#'
#' therefore, to compress the extrapolation by power \eqn{p}, use:
#'
#' \deqn{z(x) = (x + \frac{a}{p}^\frac{1}{p-1})^p + b - \frac{a}{p}^\frac{p}{p-1}}{z(x) = (x + (a/p)^(1/(p-1)))^p + b - (a/p)^(p/(p-1))}
#'
#' When extrapolating into the negative domain, just take the negative of both y and x, then call this function, then negate x and y again.
#'
#' @param x numeric vector, all >= 0
#' @param p power, in range [0, 1]
#' @param a gradient of linear extrapolation
#' @param b cap value
#'
#' @return numeric vector of compressed values \eqn{b \le z(x) \le ax + b}{b <= z(x) <= ax + b}
#'
#' @examples
#'
#' a <- 1
#' p <- 0.25
#' b <- 1
#' x <- seq(0, 5, 0.1)
#'
#' testthat::expect_true(all(rphildyerphd:::compress_extrap_z(x, p, a, b) <= a*x+b))
#' testthat::expect_true(all(rphildyerphd:::compress_extrap_z(x, p, a, b) >= b))
#' testthat::expect_true(all(rphildyerphd:::compress_extrap_z(x, p = 0, a, b) == b))
#'
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x, p = 2, a, b), "p not less than or equal to 1")
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x, p, a = -1, b), "a not greater than or equal to 0")
#' #testthat::expect_true(all(rphildyerphd:::compress_extrap_z(x, p, a = 0, b) == b))
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x, p = -1, a, b), "p not greater than or equal to 0")
#' testthat::expect_error(rphildyerphd:::compress_extrap_z(x = seq(-1, 1, 0.1), p, a, b), "Elements 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 of x >= 0 are not true")
#'
#' z <- rphildyerphd:::compress_extrap_z(x, p, a, b)
#' pl_all <- data.frame(x = x, y = a*x+b, z = z, cap = b)
#' matplot(pl_all$x, pl_all[,c("y","z", "cap")], type = "l")
#'
compress_extrap_z <- function(x, p, a, b){

  assertthat::assert_that(p >=0, p<=1)
  assertthat::assert_that(a >=0)
  assertthat::assert_that(all(x >= 0))

  if(p == 1){
    z <- a * x + b
  } else if (a == 0 | p == 0){
    z <- b
  } else {
    z <- (x + (a / p) ^ (1/(p-1)) ) ^ p + b - (a / p) ^ (p / (p-1))
  }

  return(z)

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
#' "variance" gives the unweighted variance of the points
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

  assertthat::assert_that("bootstrapGradientForest" %in% class(object))
  assertthat::assert_that(is.logical( extrap) )
  if (missing(newdata))
    newdata <- object$gf_list[[1]]$X

  assertthat::assert_that(inherits(newdata,"data.frame"))

  newnames <- names(newdata)
  pred_vars <- gfbootstrap::pred_names(object)
  if(!all(ok <- newnames %in% pred_vars)) {
    badnames <- paste(newnames[!ok], collapse=", ")
    stop(paste("the following predictors are not in any of the bootstrapGradientForests:\n\t",badnames,sep=""))
  }


  ##Calculating means and variances requires all points anyway.


  predict_groups <- expand.grid(newnames, seq_along(object$gf_list))
  predict_all <- do.call("rbind", future.apply::future_apply(predict_groups, 1, function(cg, x){
    pred <- as.character(cg[1])
    gf <- as.integer(cg[2])
    if (is.element(pred, levels(x$gf_list[[gf]]$res$var)) ){
      curve_data <- gradientForest::cumimp(x$gf_list[[gf]], pred)
    } else {
      message(paste0("Tree [", gf,
                     "], had no occurences of predictor [", pred,
                     "]"))
      return(data.frame(var = NULL, gf = NULL, x_row = NULL, x = NULL, y = NULL, stringsAsFactors = FALSE))
    }
    if(length(curve_data$x) < 2){
      message(paste0("Not using Tree [", gf,
                     "], it had only one occurences of predictor [", pred,
                     "]: x = ", curve_data$x, "\n"))
      return(data.frame(var = NULL, gf = NULL, x_row = NULL, x = NULL, y = NULL, stringsAsFactors = FALSE))
    }
    ##predictor exists
    ## get the cumimp for the predictor: curve_data
    ## figure out the range of the original data in x and y
    x_model <- range(curve_data$x)
    y_model <- range(curve_data$y)
    ## figure out the range of the new data
    ##   for x, just look at it
    x_new <- range(newdata[,pred], na.rm = TRUE)
    ## determine the range of the y data by either
    ##   if "clip", just keep old y max and min
    if (extrap == "clip"){
      interp_method <- 2 #approxfun, set extremes to max y_model
      y_new <- y_model
    } else if (extrap == "na"){
      interp_method <- 1 #approxfun, set extremes to na
      y_new <- y_model
    } else {
      interp_method <- 2
    ##   linearly extrapolate, using only the min and max of the new and old xy data
      y_new <- y_model[1] + (x_new - x_model[1]) * diff(y_model)/diff(x_model)
    }
    ## add the new extremes to the cumimp vector
    is_prepend <- FALSE
    is_append <- FALSE
    if(extrap != "na") {
      if (x_new[1] < x_model[1]) {
        is_prepend <- TRUE
        curve_data$x <- c(x_new[1], curve_data$x)
        curve_data$y <- c(y_new[1], curve_data$y)
      }
      if (x_new[2] > x_model[2]) {
        is_append <- TRUE
        curve_data$x <- c(curve_data$x, x_new[2])
        curve_data$y <- c(curve_data$y, y_new[2])
      }
    }
    ## use approxfun to interpolate along the discretized cumimp vector
    tmp_func <- approxfun(curve_data, rule = interp_method)
    ## for all extrap != "na", get the linear interpolation and apply further processing
    predicted <- tmp_func(newdata[,pred])
    if (extrap == "compress") {
      tmp_y <- predicted
      tmp_x <- newdata[, pred]
      ## from rphildyerphd::gf_extrap_compress
      ## find gradient of slope
      grad <- diff(y_model)/diff(x_model)
      ## find range of cumimp: x_model, y_model
      ## find points above the upper limit
      is_upper <- tmp_x > x_model[2]
      ## pass the points above the upper limit into compress_extrap_z
      if(any(is_upper)) {
        predicted[is_upper] <- gfbootstrap:::compress_extrap_z(tmp_x[is_upper] - x_model[2], extrap_pow, grad, y_model[2])
      }
      ## repeat for points below the lower limit
      is_lower <- tmp_x < x_model[1]
      if(any(is_lower)) {
        predicted[is_lower] <- -gfbootstrap:::compress_extrap_z(-(tmp_x[is_lower] - x_model[1]), extrap_pow, grad, -y_model[1])
      }
    }
    ##apply offset
    predicted_offset <- predicted + x$offsets[gf, pred]
    ## repeat for all models in the set
    ## with all the results, calculate and return statistics.


    return(data.frame(var = pred, gf_model = gf, x_row = seq_along(newdata[,pred]), x = newdata[,pred], y = predicted_offset, stringsAsFactors = FALSE))

  }, x = object)
  )

  ##Now generate summary results
  ##tidyverse style
  out <- data.frame(type = NA, pred = NA, x_row = NA, x = NA, y = NA, gf_model = NA)[numeric(0), ]
  all_opts <- c("mean", "variance", "points", "weight")
  if ("mean" %in% type){
    out <- rbind(out, do.call("rbind", future.apply::future_by(gf_predictions_long,
                            list(pred = gf_predictions_long$pred,
                                 x_row = gf_predictions_long$x_row),
                            function(x) {
                              data.frame(type = "mean", pred = unique(x$pred), x_row = unique(x$x_row), x = unique(x$x),
                                         y = mean(x$y), gf_model = NA)
                            }

                            ))
    )
  }
  if ("variance" %in% type){
    out <- rbind(out, do.call("rbind", future.apply::future_by(gf_predictions_long,
                            list(pred = gf_predictions_long$pred,
                                 x_row = gf_predictions_long$x_row),
                            function(x) {
                              data.frame(type = "variance",pred = unique(x$pred), x_row  = unique(x$x_row), x = unique(x$x),
                                         y = var(x$y), gf_model = NA)
                            }
      ))
    )
  }
  if ("points" %in% type){
    out <- rbind(out, data.frame(type = "points", gf_predictions_long[,c("pred", "x_row", "x", "gf_model")],
                          y = gf_predictions_long$y)
    )
  }
  # if("weight" %in% type){
  #   out["weight"] <- object$gf_list[[1]]
  # }

  class(out) <- c("data.frame", "predict.bootstrapGradientForest")
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
    stop("Every argument must be a gradientForest")

  if(is.null(gf_names <- names(gf_list)))
    gf_names <- paste("F",1:n_gf,sep="")
  if (any(empty <- gf_names==""))
    gf_names[empty] <- paste("F",1:n_gf,sep="")[empty]

  names(gf_list) <- gf_names

  n_gf_boot <- lapply(gf_list, function(gf){
    length(gf$gf_list)
  })

  combin <- data.frame(lapply(n_gf_boot, function(n, n_samp){
    sample.int(n, n_samp, replace = TRUE)
  }, n_samp = n_samp) )

  gf_combine <- apply(combin, 1, function(samp, gf_list, nbin, method, standardize){
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

  gf_combine_dist <- gfbootstrap::combine_gfbootstrap_dist(gf_combine, x_samples = x_samples)

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
pred_names_combined <- function(obj) {
  unique(do.call("c", future.apply::future_lapply(obj$gf_list,
                                                  function(x){
                                                    ##now each is a combined GF
                                                    names(x$CU)
                                                  })
                 ))
}

#' Combined gf_dist
#'
#' Calculate the distance between bootstrapped
#' samples of combined gradientForest models.
#'
#' @param gf_combine combinedBootstrapGF object
#' @param x_samples number of points along each cumimp curve for
#' deciding the best offsets
#'
#' @return data.frame, bootstraps by predictors
#'
#' @export
combine_gfbootstrap_dist <- function(gf_combine,
                                     x_samples = 100){
  assertthat::assert_that(inherits(gf_combine, "combinedBootstrapGF"))

  pred_vars <- gfbootstrap::pred_names_combined(gf_combine) 
  P <- length(pred_vars)
  K <- length(gf_combine$gf_list)

  if (!is.null(gf_combine$offsets)){
      assertthat::assert_that(nrow(gf_combine$offsets) == K)
      assertthat::assert_that(ncol(gf_combine$offsets) == P)
  } else {
      gf_combine$offsets <- as.data.frame(matrix(
        rep.int(0, K * P),
        nrow = K,
        ncol = P))
      names(gf_combine$offsets) <- pred_vars
  }

  ##generate pairs
  d_ij <- expand.grid(i = seq.int(K), j = seq.int(K))

  d_ij_diag <- d_ij[d_ij$i < d_ij$j, ]

  ##for each pair, find area between curves
  d_ij_pred <- do.call("rbind", future.apply::future_apply(d_ij_diag, 1, function(ind){
    ind <- as.data.frame(t(ind))
    i <- ind$i
    j <- ind$j
    assertthat::assert_that(inherits(gf_combine$gf_list[[i]], "combinedGradientForest"))
    assertthat::assert_that(inherits(gf_combine$gf_list[[j]], "combinedGradientForest"))
    tmp <- do.call("rbind",
            future.apply::future_lapply(pred_vars, function(pred, gf_combine){
              gf_list <- gf_combine$gf_list
              ##calculate overlap
              if (is.element(pred, names(gf_list[[i]]$CU)) ){
                x_cumimp_i <- gradientForest::cumimp(gf_list[[i]], pred, weight = gf_combine$weight)$x
              } else {
                message(paste0("Combined Tree [", i,
                               "], had no occurences of predictor [", pred,
                               "]"))
                return(data.frame(i, j, pred,  d=0))
              }
              if (is.element(pred, names(gf_list[[i]]$CU)) ){
                x_cumimp_j <- gradientForest::cumimp(gf_list[[i]], pred, weight = gf_combine$weight)$x
              } else {
                message(paste0("Combined Tree [", j,
                               "], had no occurences of predictor [", pred,
                               "]"))
                return(data.frame(i, j, pred,  d=0))
              }

              ## message(names(gf_list[[i]]$X))
              ## message(names(gf_list[[j]]$X))
              x_range <- c(max(
                  min(x_cumimp_i),
                  min(x_cumimp_j)
                ),
                min(
                  max(x_cumimp_i),
                  max(x_cumimp_j)
                )
                )
              if(x_range[2] <= x_range[1]){

                message(paste0("Combined Tree [", i,
                               "] and Combined Tree [", j,
                               "], for predictor [", pred,
                               "] had no overlap in cumimp curves"))
                return(data.frame(i, j, pred,  d=0))

              }
              ##predict x_samples points within overlapping ranges
              x_sample_points <- data.frame(seq(x_range[1], x_range[2],
                                                       length.out = x_samples))
              names(x_sample_points) <- pred
              ## message("i")
              x_i <- gradientForest::predict.combinedGradientForest(gf_list[[i]],
                                                            x_sample_points,
                                                            extrap = FALSE,
                                                            weight = gf_combine$weight)
              ## message("i passed. j")
              x_j <- gradientForest::predict.combinedGradientForest(gf_list[[j]],
                                                            x_sample_points,
                                                            extrap = FALSE,
                                                            weight = gf_combine$weight)
              ## message("j passed")
              d <-  sum((x_i + gf_combine$offsets[pred][i,]) - (x_j + gf_combine$offsets[pred][j,]))
              d <- d/x_samples
              return(data.frame(i, j, pred, d))
            }, gf_combine = gf_combine)
            )
    return(tmp)
  })
  )

  ## fill in the lower triange, d_ji
  ##take d_ij, and negate $d, swap i, j
  d_ji <- data.frame(d_ij_pred$j,
                     d_ij_pred$i,
                     d_ij_pred$pred,
                     -d_ij_pred$d)
  names(d_ji) <- names(d_ij_pred)

  d_ij_full <- rbind(d_ij_pred, d_ji)
  return(d_ij_full)

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
                           vars = gfbootstrap::pred_names_combined(x),#names(importance(x$gf_list[[1]], type = "Weighted", sorted = TRUE)),
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
  ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
  no_offset_stretch =  ggplot2::ggplot(data = curve_no_offset, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(var), scales = "free"),
  offset = ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
  offset_stretch = ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(var), scales = "free")
  ))
  }
    return(
      ggplot2::ggplot(data = curve_all, mapping = ggplot2::aes(x = x, y = y, group = gf, color = as.factor(gf))) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x")
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
#' "variance" gives the unweighted variance of the points
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
#' @param extrap can be TRUE or FALSE
#'
#' If extrap = FALSE, return NA outside of the training values. extrap = FALSE only uses
#' models that observed the training data when calculating mean and variance
#'
#' If extrap = TRUE, then the curves are extrapolated if the new value lies outside
#' the training values.
#' @param extrap_pow is only used with extrap = TRUE, and sets the
#' compression power. 1/4 gives the 4th root. 1 produces linear extrapolation
#' and 0 clips the extrapolation to the min and max values of Y.
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
                                            extrap_pow = 1/4,
                                            ...){


  assertthat::assert_that(inherits(object, "combinedBootstrapGF"))
  assertthat::assert_that(extrap %in% c(TRUE, FALSE))
  assertthat::assert_that(length(extrap) == 1)
  if (missing(newdata)){
    newdata <- object$gf_list[[1]]$X[,gfbootstrap::pred_names_combined(object)]
  } else {
    newnames <- names(newdata)
    gf_preds <- gfbootstrap::pred_names_combined(object)
    if(!all(ok <- newnames %in% gf_preds)) {
      badnames <- paste(newnames[!ok], collapse=", ")
      stop(paste("the following predictors are not in any of the bootstrapGradientForests:\n\t",badnames,sep=""))
    }
  }
  assertthat::assert_that(inherits(newdata,"data.frame"))

  ##Calculating means and variances requires all points anyway.


  predict_groups <- expand.grid(newnames, seq_along(object$gf_list))
  predict_all <- do.call("rbind", future.apply::future_apply(predict_groups, 1, function(cg, x){
    pred <- as.character(cg[1])
    gf <- as.integer(cg[2])
    if (is.element(pred, names(x$gf_list[[gf]]$CU)) ){
      curve_data <- gradientForest::cumimp(x$gf_list[[gf]], pred, weight = x$weight)
    } else {
      message(paste0("Tree [", gf,
                     "], had no occurences of predictor [", pred,
                     "]"))
      return(data.frame(pred = NULL, gf = NULL, x_row = NULL, x = NULL, y = NULL, stringsAsFactors = FALSE))
    }
    if(length(curve_data$x) < 2){
      message(paste0("Not using Tree [", gf,
                     "], it had only one occurences of predictor [", pred,
                     "]: x = ", curve_data$x, "\n"))
      return(data.frame(pred = NULL, gf = NULL, x_row = NULL, x = NULL, y = NULL, stringsAsFactors = FALSE))
    }
    ##predictor exists
    ## get the cumimp for the predictor: curve_data
    ## figure out the range of the original data in x and y
    x_model <- range(curve_data$x)
    y_model <- range(curve_data$y)
    ## figure out the range of the new data
    ##   for x, just look at it
    x_new <- range(newdata[,pred], na.rm = TRUE)
    ## determine the range of the y data by either
    ##   if "clip", just keep old y max and min
    if (extrap) {
      interp_method <- 2 #approxfun, set extremes to max y_model
    } else {
      interp_method <- 1 #approxfun, set extremes to na
    }

    if (extrap_pow == 0 || !extrap){
      y_new <- y_model
    } else {
      ##   linearly extrapolate, using only the min and max of the new and old xy data
      y_new <- y_model[1] + (x_new - x_model[1]) * diff(y_model)/diff(x_model)
    }
    ## add the new extremes to the cumimp vector
    is_prepend <- FALSE
    is_append <- FALSE
    if(extrap) {
      if (x_new[1] < x_model[1]) {
        is_prepend <- TRUE
        curve_data$x <- c(x_new[1], curve_data$x)
        curve_data$y <- c(y_new[1], curve_data$y)
      }
      if (x_new[2] > x_model[2]) {
        is_append <- TRUE
        curve_data$x <- c(curve_data$x, x_new[2])
        curve_data$y <- c(curve_data$y, y_new[2])
      }
    }
    ## use approxfun to interpolate along the discretized cumimp vector
    tmp_func <- approxfun(curve_data, rule = interp_method)
    ## for all extrap != "na", get the linear interpolation and apply further processing
    predicted <- tmp_func(newdata[,pred])
    if (extrap) {
      tmp_y <- predicted
      tmp_x <- newdata[, pred]
      ## from rphildyerphd::gf_extrap_compress
      ## find gradient of slope
      grad <- diff(y_model)/diff(x_model)
      ## find range of cumimp: x_model, y_model
      ## find points above the upper limit
      is_upper <- tmp_x > x_model[2]
      ## pass the points above the upper limit into compress_extrap_z
      if(any(is_upper)) {
        predicted[is_upper] <- gfbootstrap:::compress_extrap_z(tmp_x[is_upper] - x_model[2], extrap_pow, grad, y_model[2])
      }
      ## repeat for points below the lower limit
      is_lower <- tmp_x < x_model[1]
      if(any(is_lower)) {
        predicted[is_lower] <- -gfbootstrap:::compress_extrap_z(-(tmp_x[is_lower] - x_model[1]), extrap_pow, grad, -y_model[1])
      }
    }
    ##apply offset
    predicted_offset <- predicted + x$offsets[gf, pred]
    ## repeat for all models in the set
    ## with all the results, calculate and return statistics.


    return(data.frame(pred = pred, gf_model = gf, x_row = seq_along(newdata[,pred]), x = newdata[,pred], y = predicted_offset, stringsAsFactors = FALSE))

  }, x = object)
  )

  ##Now generate summary results
  ##tidyverse style
  out <- data.frame(type = NA, pred = NA, x_row = NA, x = NA, y = NA, gf_model = NA)[numeric(0), ]
  all_opts <- c("mean", "variance", "points", "weight")
  if ("mean" %in% type){
    out <- rbind(out, do.call("rbind", future.apply::future_by(predict_all,
                            list(pred = predict_all$pred,
                                 x_row = predict_all$x_row),
                            function(x) {
                              data.frame(type = "mean", pred = unique(x$pred), x_row = unique(x$x_row), x = unique(x$x),
                                         y = mean(x$y), gf_model = NA)
                            }

                            ))
    )
  }
  if ("variance" %in% type){
    out <- rbind(out, do.call("rbind", future.apply::future_by(predict_all,
                            list(pred = predict_all$pred,
                                 x_row = predict_all$x_row),
                            function(x) {
                              data.frame(type = "variance",pred = unique(x$pred), x_row  = unique(x$x_row), x = unique(x$x),
                                         y = var(x$y), gf_model = NA)
                            }
      ))
    )
  }
  if ("points" %in% type){
    out <- rbind(out, data.frame(type = "points", predict_all[,c("pred", "x_row", "x", "gf_model")],
                          y = predict_all$y)
    )
  }
  # if("weight" %in% type){
  #   out["weight"] <- object$gf_list[[1]]
  # }

  class(out) <- c("data.frame", "predict.combinedBootstrapGF")
  return(out)

}

#' Huberts Gamma statistic
#'
#' Calculates Huberts Gamma statistic
#'
#' norm_z rescales both x and y
#' by subtracting the mean and dividing by the
#' standard deviation.
#' The norm_z variant is from Tseng and Kao 2003.
#'
#' The original Huberts Gamma statistic (Hubert and Levin 1976)
#' uses x and y as is.
#'
#' Zhao et al. 2006 proposes another variant, where
#' y is defined as the distance between the cluster
#' centers of site i and site j.
#' This variant is harder to calculate and gives
#' a knee rather than a peak.
#' This variant is not implemented.
#'
#' This version requires a full matrix for both sim_mat and
#' member_mat. Using other forms is much more complicated.
#'
#' @param sim_mat a square matrix equivalent, all values in the range 0 to 1
#' @param member_mat square matrix equivalent, 0 if a pair of sites are in different
#' clusters, 1 if the pair are in the same cluster.
#' @param norm_z whether to normalize the sim_mat and member_mat into z-scores. See Tseng and Kao 2003
#'
#' @return the Hubert Gamma score
#'
#' @export
hubert_gamma <- function(sim_mat, member_mat, norm_z = TRUE){
  assertthat::assert_that(assertthat::are_equal(dim(sim_mat), dim(member_mat)))

  if(norm_z){
    mean_s <- mean(sim_mat)
    sd_s <- sd(sim_mat)
    mean_m <- mean(member_mat)
    sd_m <- sd(member_mat)

    member_mat <- (member_mat - mean_m) / sd_m
    sim_mat <- (sim_mat - mean_s) / sd_s
  } else {
    mean_m <- 0
    sd_m <- 1
    mean_s <- 0
    sd_s <- 1
  }

  h_gamma <- (1/sum(upper.tri(sim_mat))) *
    sum((upper.tri(sim_mat)*sim_mat) *
          member_mat)
  return(h_gamma)
}


