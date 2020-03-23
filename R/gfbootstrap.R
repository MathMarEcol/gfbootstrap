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
#' @param x
#' @param predictor.vars
#' @param response.vars
#' @param nboostraps number of bootstrapped models to fit
#' @param nsamples number of samples taken at each bootstrapped model fit
#' @param mtry
#' @param transform
#' @param maxLevel
#' @param corr.threshold
#' @param compact
#' @param nbin
#' @param trace
#'
#' @import future.apply
#'
#' @examples
#'
#' library(gradientForest)
#' data(CoMLsimulation)
#' predictor.vars <- colnames(Xsimulation)
#' response.vars <- colnames(Ysimulation)
#' demobootstrap <- boostrapGradientForest(x = data.frame(Ysimulation,Xsimulation),
#' predictor.vars = predictor.vars,
#' response.vars = response.vars,
#' nbootstrap = 100, #small, should be 500 or 1000 in a real experiment
#' compact = T, nbin = 200,
#' transform = NULL,
#' corr.threshold = 0.5,
#' maxLevel = floor(log2(length(.y)*0.368/2)),
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
        message("GF model failed to fit, restarting")
        message(e)
        return(NULL)
    }, warning = function(e){
      message("GF model failed to fit, restarting")
      message(e)
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

  })

  ##Calculate offsets
    ##the optimal offset requires knowing the distance between curves for every tree in gfbootstrap
    ##
  df_dist <- gfbootstrap_dist(gf_bootstrap,
                              x_samples = 100)
  offsets <- gfbootstrap_offsets(df_dist)

  ##The offsets are applied to the cumimp curves,
  ##but the cumimp curves are calculated
  ##on demand by gradientForest.
  ##HI will just keep the offsets beside the gf models,
  ## and apply the offsets on demand.
  out <- list(
    gf_list = gf_bootstrap,
    offset = offsets
  )
  class(out) <- "bootstrapGradientForest"
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
#' If an offsets vector is provided, it is included in the result.
#'
#' x_samples integer number of points sampled along each cumulative importance curve, evenly spaced.
#'
#' Returns a data.frame of areas between curves for each pair and each predictor.
gfbootstrap_dist <- function(gf_list,
                             offsets = NULL,
                             x_samples = 100){

  assertthat::assert_that("gradientForest" %in% class(gf_list[[1]]))
  pred_vars <- names(gf_list[[1]]$X)
  P <- length(pred_vars)
  K <- length(gf_list)

  if(!is.null(offsets)){
    assertthat::assert_that(nrow(offsets) == K)
    assertthat::assert_that(ncol(offsets) == P)
  } else {
    offsets <- matrix(
      rep.int(0, K*P),
      nrow = K,
      ncol = P
    )
  }

  offsets <- as.data.frame(offsets)
  names(offsets) <- pred_vars
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

#' Calculate opimal offsets
#'
#' The distance values returned by `gfbootstrap_dist()`
#' are used to calculate the optimal offsets.
#' Given a predictor, an offset is calculated for each
#' gf cumulative importance curve,
#' such that the total area between the curves is
#' minimised.
#'
#' \param x the data.frame returned by gfbootstrap_dist()
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

#' Predict bootstrapped GF object
#'
#' Given a fitted gfbootstrap object, predict new points.
#'
#' Unlike most predict functions, you get multiple output rows per input row.
#' The id_cols parameter allows you to specify which columns should be used as identifiers for each input row.
#' The total number of output rows is gfbootstrap$nbootstraps * nrow(predict_points)
#'
#'


#' plot cumimp for bootStrapGradientForest
#'
#' Returns a ggplot object, ready for plotting
gg_bootstrapGF <- function(x,
                           vars = names(importance(x$gf_list[[1]], type = "Weighted", sorted = TRUE)),
                           n_curves = 10,
                           debug = TRUE) {
  assertthat::assert_that("bootstrapGradientForest" %in% class(x))
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
    
    if (is.element(pred, levels(x$gf_list[[gf]]$res$var)) ){
      curve_data <- gradientForest::cumimp(x$gf_list[[gf]], pred)
    } else {
      message(paste0("gg_bootstrapgf: Tree [", gf,
                     "], had no occurences of predictor [", pred,
                     "]"))
      return(data.frame(var = NULL, gf = NULL, x = NULL, y = NULL, stringsAsFactors = FALSE))
    }
    curve_data$y <- curve_data$y + x$offset[gf, pred]
    return(data.frame(var = pred, gf = gf, x = curve_data$x, y = curve_data$y, stringsAsFactors = FALSE))

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
        return(data.frame(var = NULL, gf = NULL, x = NULL, y = NULL, stringsAsFactors = FALSE))
      }

      return(data.frame(var = pred, gf = gf, x = curve_data$x, y = curve_data$y, stringsAsFactors = FALSE))

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
#'
#' `type` is a vector containing elements from
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
#' extrap can take on 4 values, c("clip", "linear", "compress", "na")
#'
#' "linear" extrapolate each model linearly using average gradient
#'
#' "clip" returns the maximum or minimum Y if X exceeds training values
#'
#' "compress" extrapolates, but applies compression so predicted
#' values approach an asymptote as X moves veary far from training values.
#' extrap_pow is only used with "compress", and sets the
#' compression power. 1/4 gives the 4th root, 1 is equivalent to "linear"
#' and 0 is equivalent to "clip"
#'
#' TODO: just use "compress" with edge cases for linear and clip
#' TODO: This code is rather spaghetti, eg. activity 4 depends on 2 and 3,
#' then activity 5 depends on 1 and 4
#'
#' "na" return Na outside of the training values. "na" only uses
#' models that observed the training data when calculating mean and variance
#' ... arguments passed to cumimp
#' If `type` contains more than one string, a list of
#' data.frames will be returned, otherwise just a data.frame.
predict.bootstrapGradientForest <- function(object,
                                            newdata,
                                            type = c("mean"),
                                            extrap="compress",
                                            extrap_pow = 1/4,
                                            ...){

  assertthat::assert_that("bootstrapGradientForest" %in% class(object))
  assertthat::assert_that(extrap %in% c("linear", "clip", "compress", "na"))
  assertthat::assert_that(length(extrap) == 1)
  if (missing(newdata))
    newdata <- object$gf_list[[1]]$X

  assertthat::assert_that(inherits(newdata,"data.frame"))

  newnames <- names(newdata)
  gf_preds <- pred_names(object)
  if(!all(ok <- newnames %in% gf_preds)) {
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
        predicted[is_upper] <- compress_extrap_z(tmp_x[is_upper] - x_model[2], extrap_pow, grad, y_model[2])
      }
      ## repeat for points below the lower limit
      is_lower <- tmp_x < x_model[1]
      if(any(is_lower)) {
        predicted[is_lower] <- -compress_extrap_z(-(tmp_x[is_lower] - x_model[1]), extrap_pow, grad, -y_model[1])
      }
    }
    ##apply offset
    predicted_offset <- predicted + x$offset[gf, pred]
    ## repeat for all models in the set
    ## with all the results, calculate and return statistics.


    return(data.frame(var = pred, gf_model = gf, x_row = seq_along(newdata[,pred]), x = newdata[,pred], y = predicted_offset, stringsAsFactors = FALSE))

  }, x = object)
  )

  ##Now generate summary results
  ##tidyverse style
  out <- data.frame(type = NA, var = NA, x_row = NA, x = NA, y = NA, gf_model = NA)[numeric(0), ]
  all_opts <- c("mean", "variance", "points", "weight")
  if ("mean" %in% type){
    out <- rbind(out, do.call("rbind", future.apply::future_by(predict_all,
                            list(var = predict_all$var,
                                 x_row = predict_all$x_row),
                            function(x) {
                              data.frame(type = "mean", var = unique(x$var), x_row = unique(x$x_row), x = unique(x$x),
                                         y = mean(x$y), gf_model = NA)
                            }

                            ))
    )
  }
  if ("variance" %in% type){
    out <- rbind(out, do.call("rbind", future.apply::future_by(predict_all,
                            list(var = predict_all$var,
                                 x_row = predict_all$x_row),
                            function(x) {
                              data.frame(type = "variance",var = unique(x$var), x_row  = unique(x$x_row), x = unique(x$x),
                                         y = var(x$y), gf_model = NA)
                            }
      ))
    )
  }
  if ("points" %in% type){
    out <- rbind(out, data.frame(type = "points", predict_all[,c("var", "x_row", "x", "gf_model")],
                          y = predict_all$y)
    )
  }
  # if("weight" %in% type){
  #   out["weight"] <- object$gf_list[[1]]
  # }

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
#' 
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

aff_func <- function(x, u, sim_mat){
  ##Takes the vector x of rows and the row u of sim_match
  ##and returns the affinity from u to x
  assertthat::assert_that(length(x) == nrow(sim_mat))
  assertthat::assert_that(is.logical(x))
  assertthat::assert_that(is.numeric(u))

  ##Simple, this is just elements of sim_mat
  return(sim_mat[u, ]*x)
}

cast_alg <- function(sim_mat, aff_thres){
##rather than wrestle with integer sets, I will use boolean vectors
  
  ##initialise, no clusters, and all elements belong to no clusters
  n <- nrow(sim_mat)
  clust <- list()
  spares <- rep(TRUE, n)
  
  clust_id <- 1
  ##now we keep working until spares is all false
  while(any(spares)) {
    debug_spares_loop <- sum(spares)
    debug_spares_expected <- debug_spares_loop
    new_clust <- rep(FALSE, n)

    aff_local <- rep(0, n)

    ##Do the next steps until stability is reached
    repeat{
      is_change <- FALSE
      
      ##find initial cluster point
      if(all(!new_clust)){
        maxima <- which(spares)
        new_elem <- maxima[sample.int(length(maxima), 1)] ##intention: take one of the max at random
        new_clust[new_elem] <- TRUE
        spares[new_elem] <- FALSE
        is_change <- TRUE
        debug_spares_expected <- debug_spares_expected - 1
        assertthat::assert_that(debug_spares_expected == sum(spares))
        ##update aff_local
        ## aff_local[new_clust | spares] <- aff_local[new_clust | spares] + aff_func(new_clust | spares, new_elem, sim_mat)
        aff_local <- aff_local + aff_func(new_clust | spares, new_elem, sim_mat)
      }

      ##addition stage
      ##affinity seems to be used to bias group size
      ##Strong affinity is needed to join a large group

      ##Short-circuit evaluation, won't try to find max of an empty set
      while(any(spares) && (max(aff_local[spares]) >= aff_thres*mean(aff_local[new_clust]))){
        maxima <- which(aff_local == max(aff_local[spares]))
        new_elem <- maxima[sample.int(length(maxima), 1)] ##intention: take one of the max at random
        new_clust[new_elem] <- TRUE
        spares[new_elem] <- FALSE
        debug_spares_expected <- debug_spares_expected - 1
        assertthat::assert_that(debug_spares_expected == sum(spares))
        is_change <- TRUE
        ##update aff_local
        ## aff_local[new_clust | spares] <- aff_local[new_clust | spares] + aff_func(new_clust | spares, new_elem, sim_mat)
        aff_local <- aff_local + aff_func(new_clust | spares, new_elem, sim_mat)
      }

      ##Removal stage
      while(any(new_clust) && (min(aff_local[new_clust]) < aff_thres*mean(aff_local[new_clust]))){
        minima <- which(aff_local == min(aff_local[new_clust]))
        new_elem <- minima[sample.int(length(minima), 1)] ##intention: take one of the max at random
        new_clust[new_elem] <- FALSE
        spares[new_elem] <- TRUE
        debug_spares_expected <- debug_spares_expected + 1
        assertthat::assert_that(debug_spares_expected == sum(spares))
        is_change <- TRUE
        ##update aff_local
        aff_local <- aff_local - aff_func(new_clust | spares, new_elem, sim_mat)
      }

      if(!is_change){
        break
      }
    }

    assertthat::assert_that(debug_spares_expected == (debug_spares_loop - sum(new_clust)))
    clust[[clust_id]] <- new_clust
    clust_id <- clust_id + 1


  }

  return(clust)
}


cast_stabilize <- function(cast_obj, aff_thres, sim_mat, max_iter = 20){
  iter <- 1

  while(iter <= max_iter){
    ##For each vertex, find affinity to other clusters
    ##Cast_obj is not huge, but constantly recreating it when I am just
    ##flipping bools seems wasteful
    ##This approach tests each site, and updates clustering at end
    updates <- lapply(seq_along(cast_obj[[1]]), function(u, cast_obj, sim_mat){
      clust_id <- which(sapply(cast_obj, function(clust, u){
        clust[u]
      }, u = u))
      assertthat::assert_that(length(clust_id) == 1)
      ##u belongs to clust_id
      clust_aff <- lapply(cast_obj, function(clust, u, sim_mat){
          ##get the affinity to each cluster
          if(any(clust)){
            return(mean(aff_func(clust, u, sim_mat)[clust]))
          } else {
            ##empty clusters have 0 affinity
            return(0)
          }
        }, u = u, sim_mat = sim_mat)
        if(which.max(clust_aff) != clust_id){
          return(data.frame(u = u, old = clust_id, new = which.max(clust_aff)))
        } else {
          return(NULL)
        }
      }, cast_obj = cast_obj, sim_mat = sim_mat)
    updates <- do.call("rbind", updates)
    ##Apply updates
    if(!is.null(updates)){
      message("iteration [", iter, "] reassigned [", nrow(updates),"] samples")
      for(upd in 1:nrow(updates) ){
        cast_obj[[updates[upd,"old"]]][updates[upd,"u"]] <- FALSE
        cast_obj[[updates[upd,"new"]]][updates[upd,"u"]] <- TRUE
      }
    } else {
      break
    }
    iter <- iter + 1
  }
  return(cast_obj)
}

##within cluster affinity
##returns a list
aff_clust_inner <- function(cast_obj, sim_mat){
  lapply(seq_along(cast_obj), function(clust, cast_obj, sim_mat){
    elems <- cast_obj[[clust]]
    elem_cor <- sim_mat[elems,elems]
    if(sum(elems) > 1){
      elems_mean_aff <- rowSums(elem_cor)/sum(elems)
      mean_aff <- mean(elems_mean_aff)
    } else {
      mean_aff <- elem_cor
    }
    return(mean_aff)
  }, cast_obj = cast_obj, sim_mat = sim_mat)
}

##Across cluster affinity
##An assymetric distance, where the distance from a to b
##is the average affinity of elements in a to elements in b
##returns a long form dataframe
aff_cluster_between <- function(cast_obj, sim_mat){
  pairs <- expand.grid(1:length(cast_obj), 1:length(cast_obj))
  affs <- apply(pairs, 1, function(p, cast_obj, sim_mat){
    elems_a <- cast_obj[[p[1]]]
    elems_b <- cast_obj[[p[2]]]
    if(sum(elems_a) > 1 && sum(elems_b) > 1){
      return(mean(rowSums(sim_mat[elems_a, elems_b])/sum(elems_b)))
    } else {
      return(mean(sim_mat[elems_a, elems_b]))
    }

  }, cast_obj = cast_obj, sim_mat = sim_mat)
  return(data.frame(pairs, affs))
}

##Assign site to cluster
##takes an affinity matrix, 
## rows are sites to assign to clusters,
## and cols are sites that make up the cast_obj
## so the dimensions of the matrix
## are n by length(cast_obj[[1]])
## the rows are not expected to exist in
## cast_obj

##find which clusters a site could belong, given aff_thres
##returns a data.frame with site, cluster, affinity
predict_clust <- function(cast_obj, new_sim_mat){
  do.call("rbind", lapply(1:length(cast_obj), function(clust, cast_obj, new_sim_mat){
    if(sum(cast_obj[[clust]]) > 1){
      return(data.frame(x_row = 1:nrow(new_sim_mat),
                        clust = clust,
                        aff = rowSums(new_sim_mat[, cast_obj[[clust]]])/sum(cast_obj[[clust]])))
    } else {
      return(data.frame(x_row = 1:nrow(new_sim_mat),
                        clust = clust,
                        aff = new_sim_mat[, cast_obj[[clust]]]))
    }
  }, cast_obj = cast_obj, new_sim_mat = new_sim_mat))
}
##heirarchical grouping of cast objects
##the aff_thres_list will be sorted
##this function is recursive
cast_h <- function(sim_mat, aff_thres_list, stabilize = FALSE, max_stab_iter = 20){
  assertthat::assert_that(nrow(sim_mat) > 1)
  aff_thres_list <- sort(aff_thres_list, decreasing = TRUE)
  next_level <- cast_alg(sim_mat = sim_mat, aff_thres = aff_thres_list[1])
  if(stabilize){
    next_level <- cast_stabilize(sim_mat = sim_mat, cast_obj = next_level,
                                 max_iter = max_stab_iter,
                                 aff_thres = aff_thres_list[1])
  }
  if(length(aff_thres_list) > 1){
    tmp <- aff_cluster_between(cast_obj = next_level, sim_mat = sim_mat)
     new_sim_mat <- matrix(tmp$affs, nrow = length(next_level), ncol = length(next_level))
    return(c(cast_h(sim_mat = new_sim_mat,
                    aff_thres_list = aff_thres_list[-1],
                    stabilize = stabilize,
                    max_stab_iter = max_stab_iter),
             list(list(clust = next_level, affs = new_sim_mat)) ))
  } else {
    tmp <- aff_cluster_between(cast_obj = next_level, sim_mat = sim_mat)
    new_sim_mat <- matrix(tmp$affs, nrow = length(next_level), ncol = length(next_level))
    return(list(list(clust = next_level, affs = new_sim_mat)))
  }
}

##unrolls the cast clusters to the desired depth
##also recursive
cast_h_reorder <- function(cast_h_obj, depth = length(cast_h_obj), reordering = NULL){
  assertthat::assert_that(depth > 0)
  assertthat::assert_that(depth <= length(cast_h_obj))

  if(is.null(reordering)){
    reordering <- 1:length(cast_h_obj[[1]]$clust)
  }
  layers <- seq(1, depth )
  reordered <- do.call("c", lapply(cast_h_obj[[1]]$clust[reordering], which))
  if(depth > 1){
    return(cast_h_reorder(cast_h_obj[-1], depth = depth - 1, reordering = reordered))
  } else {
    return(reordered)
  }
}
##Scratch area
no_eval <- function(){

library(gradientForest)
data(CoMLsimulation)
predictor.vars <- colnames(Xsimulation)
response.vars <- colnames(Ysimulation)
demobootstrap <- bootstrapGradientForest(x = data.frame(Ysimulation,Xsimulation),
                                        predictor.vars = predictor.vars,
                                        response.vars = response.vars,
                                        nbootstrap = 50, #small, should be 500 or 1000 in a real experiment
                                        compact = T, nbin = 200,
                                        transform = NULL,
                                        corr.threshold = 0.5,
                                        maxLevel = floor(log2(length(Ysimulation)*0.368/2)),
                                        trace = TRUE
                                        )
library(ggplot2)
plots <- gg_bootstrapGF(demobootstrap, n_curves = 10, debug = TRUE)
lapply(names(plots), function(pl, gg_boot) {
  ggsave(filename = paste0("demobootstrap_", pl, "_curves.png"), plot = gg_boot[[pl]],
  units = "cm",
  width = 16,
  height = 9,
  dpi = 300,
  scale = 2
  )
}, gg_boot = plots)

##figure out what happens when a varirable is unused
## nbootstrap <- 100

##   gf_bootstrap <- future.apply::future_lapply(X = 1:nbootstrap, function(i){

##     ##Fit GF with a single tree, but otherwise identical calls.

##     ##In bootstrapping, the input data are randomly resampled.
##     ##GF internally bootstraps each tree as part of the random forest algorithm
##     ##therefore, each loop will produce a bootstrapped GF model

##     valid_model <- FALSE
##     tries <- 1
##     while(!valid_model & tries < 10){
##       if(tries > 1) message(paste0("i: [", i, "], try: ", tries))
##       valid_model <- TRUE
##       gf_list <-  tryCatch({
##         gradientForest::gradientForest(data = data.frame(Ysimulation,Xsimulation),
##                                                 predictor.vars=predictor.vars,
##                                                 response.vars=response.vars,
##                                                 ntree = 1,
##                                                 mtry=NULL,
##                                                 transform=NULL,
##                                                 maxLevel=0,
##                                                 corr.threshold=0.5,
##                                                 compact=TRUE,
##                                                 nbin=200,
##                                                 trace=FALSE
##                                                 )
##       }, error = function(e){
##         message(e)
##         message("GF model failed to fit, i: [", i, "], tries: [", tries,"]")
##         tries <- tries + 1
##         valid_model <- FALSE
##       })
##     }
##     ## message(class(gf_list))
##     ## message(unique(gf_list$res$var))
##     ## message(length(unique(gf_list$res$var)))
##     if(length(unique(gf_list$res$var)) < length(predictor.vars)){
##       message("Some vars not found, i: [", i, "], vars: [", unique((gf_list$res$var)),"]")
##     }
##     return(gf_list)

##   })

## x_cumimp_i <- gradientForest::cumimp(gf_bootstrap[[30]], "D")$x

## test out predicted gf
##use demobootstrap
## predict with random data?
new_x <- data.frame(lapply(1:10,
                           function(x){
                             out <- list()
                             out[[LETTERS[x]]]<- runif(n = 100, min =  min(Xsimulation[, x])-0.5, max = max(Xsimulation[, x])+0.5)
                           return(out)
}))



predicted_data <- predict.bootstrapGradientForest(demobootstrap,
                                                  new_x,
                                                  type = c("mean", "variance", "points"),
                                                  extrap="compress",
                                                  extrap_pow = 1/16
                                                  )

library(ggplot2)
ggsave("demobootstrap_variance.png", ggplot2::ggplot(data = predicted_data[predicted_data$type == "variance",], mapping =
                                                                                                                  ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
  units = "cm",
  width = 16,
  height = 9,
  dpi = 300,
  scale = 2
  )


ggsave("demobootstrap_mean.png", ggplot2::ggplot(data = predicted_data[predicted_data$type == "mean",], mapping =
                                                       aes(x = x, y = y)) +
         ggplot2::geom_line() +
         ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
)

pred_mean_var <- merge(predicted_data[predicted_data$type == "mean", c("var",  "x", "y" )],
                       predicted_data[predicted_data$type == "variance",c("var",  "x", "y")],
                       by = c("var", "x"))
names(pred_mean_var)[3:4] <- c("mean", "variance")
pred_mean_var["upper"] <-pred_mean_var$mean + sqrt(pred_mean_var$variance)
pred_mean_var["lower"] <-pred_mean_var$mean - sqrt(pred_mean_var$variance)
ggsave("demobootstrap_mean_variance.png", ggplot2::ggplot(data =pred_mean_var[pred_mean_var$var == "A",], mapping =
                                                            aes(x = x, y = mean)) +
         ggplot2::geom_line() +
         ggplot2::geom_line(data = pred_mean_var[pred_mean_var$var == "A",c("x", "upper")], mapping = aes(x = x, y = upper),
                              colour = "red") +
         ggplot2::geom_line(data = pred_mean_var[pred_mean_var$var == "A",c("x", "lower")], mapping = aes(x = x, y = lower),
                            colour = "red"),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
)

ggsave("demobootstrap_points.png", ggplot2::ggplot(data = predicted_data[predicted_data$type == "points",], mapping =
                                                       aes(x = x, y = y, group = gf_model, color = as.factor(gf_model) )) +
         ggplot2::geom_line(size = 0.1) +
         ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
)
ggsave("demobootstrap_points_quantile.png", ggplot2::ggplot(data = predicted_data[predicted_data$type == "points",], mapping =
                                                     aes(x = x, y = y )) +
         ggplot2::geom_quantile(quantiles = c(0.05, 0.5, 0.95), method = "rqss", lambda = 0.1) +
         ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
)
ggsave("demobootstrap_points_smooth.png", ggplot2::ggplot(data = predicted_data[predicted_data$type == "points",], mapping =
                                                              aes(x = x, y = y )) +
         ggplot2::geom_smooth() +
         ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
)

# ggsave("demobootstrap_points_stat.png", ggplot2::ggplot(data = predicted_data[predicted_data$type == "points",], mapping =
#                                                            aes(x = x, y = y)) +
#           ggplot2::stat_summary( geom="line") +
#           ggplot2::facet_wrap("var", scales = "free_x")
# )


##now what?
##cluster?
## Hotellings! Probably use my block hotellings method, given that I already have mean and variance.
##Then cluster with WGCNA


##Hotellings matrix now.

## there are a few packages that implement hotellings T^2 test, but most of them assume a point
##cloud, which was why I wrote a method that applied pre-existing means and variance
## values.
library(rmethods)
##need to make a wide form matrix for hotellings_bulk
library(tidyr) #I am not using this within my package, so it is not a dependency.
predict_wide <- tidyr::pivot_wider(data = predicted_data[predicted_data$type %in% c("mean", "variance"),],
                                         names_from =  c("var"),
                                         values_from = c("y"),
                                         id_cols = c("x_row", "type")  )
pred_vars <- pred_names(demobootstrap)
predict_wide_means <- predict_wide[predict_wide$type == "mean",]
predict_wide_variance <- predict_wide[predict_wide$type == "variance",]
p_mat <- rmethods:::hotellings_bulk(
             means = predict_wide_means[, pred_vars],
             res_sq = predict_wide_variance[, pred_vars]
           )

p_mat_long <- tidyr::pivot_longer(data = data.frame(x = as.numeric(row.names(p_mat)), p_mat),
                                  cols = -x,
                                  names_to = c("y"),
                                  names_prefix = "X",
                    values_to = c("p"))
p_mat_long$y <- as.integer(p_mat_long$y)


ggplot2::ggsave("p_mat.png",
                ggplot2::ggplot(data = p_mat_long,
                                mapping = aes(x = x, y = y, fill = p)) +
                ggplot2::geom_tile(),
                units = "cm",
                width = 16,
                height = 9,
                dpi = 300,
                scale = 2
                )


##Bring in WGCNA!

library(WGCNA)

##following the WGCNA tutorials

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(p_mat, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
## this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
## Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##My results show that I should use 4
pow <- 4
net = blockwiseModules(p_mat,
                       power = pow,
                       TOMType = "unsigned",
                       minModuleSize = 5, #changed from default, I only have 100 samples
                       reassignThreshold = 0,
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "generated_data",
                       verbose = 3)

table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

nGenes = ncol(p_mat)
nSamples = nrow(p_mat)

dissTOM = 1-TOMsimilarityFromExpr(p_mat, power = pow)
plotTOM = dissTOM^7
diag(plotTOM) = NA
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

##Now to repeat on Pilbara
##recompute
## gf_model_hash <- list( dir = "/vmshare/phd/projects/aus_bioregions/experiments/2019-04-24-1424_bioregion_pilbara_more_datasets/archivist/", hash = "29833da39091a8846e5210e764297928") 
library(foreach)
library(ggplot2)
dataset_names <- c( "DEC_Rfs",
                   "DOF_EGsvy",
                   "DOF_Trap",
                   "DOF_Trawl",
                   "NWS_EoT",
                   "Sled",
                   "Soviet_Trawl",
                   "Trawl",
                   "Warehouse_Trl"
                   )

datasets <- foreach(dataset = dataset_names, .combine = c) %do% {
  ret <- list()
  ret[[dataset]] <- list()
  tmp <- load(paste0("/vmshare/phd/qris_sandbox/Q1215/phildyer_wip/Roland_Demersal_dataset/", dataset, "_X.Rdata"))
  ret[[dataset]][["X"]] <- get(tmp)
  tmp <- load(paste0("/vmshare/phd/qris_sandbox/Q1215/phildyer_wip/Roland_Demersal_dataset/", dataset, "_Y.Rdata"))
  ret[[dataset]][["Y"]] <- get(tmp)
  return(ret)
}

load("/vmshare/phd/qris_sandbox/Q1215/phildyer_wip/Roland_Demersal_dataset/Phys_grid.RData")
env_vars <-  c("MT_SST_AV","CRS_NO3_AV","CRS_T_AV","CRS_S_SR","CRS_T_SR")
## env_vars <-  c(
##   "GorgBATHY", "GorgSLOPE",  "RBN_BSTRESS", "dbS_CRBNT", 
##   "dbS_GRAVEL", "dbS_SAND", "dbS_MUD", "dbS_ROCK", "dbS_GRNSZ", 
##   "dbS_SORTG", "CRS_NO3_AV", "CRS_NO3_SR", "CRS_PO4_AV", "CRS_PO4_SR", 
##   "CRS_O2_AV", "CRS_O2_SR", "CRS_S_AV", "CRS_S_SR", "CRS_T_AV", 
##   "CRS_T_SR", "CRS_SI_AV", "CRS_SI_SR", "SW_K490_AV", "SW_K490_SR", 
##   "SW_CHLA_AV", "SW_CHLA_SR", "MT_SST_AV", "MT_SST_SR", "SW_BIR_AV", 
##   "SW_BIR_SR", "VGPM_AV", "VGPM_SR", "PAR_AV", "PAR_SR", 
##   "EPOC_AV", "EPOC_SR", "EFF_COVER")

                                        #These variables have only a single value across all sites
                                        #"TERAN_CHAN", "TERAN_PASS", "TERAN_PEAK", 
                                        #"TERAN_PIT", "TERAN_PLAN", "TERAN_RIDG""FT_AVG_EFF",

                                        #These variables weren't used by Roland
                                        #"PT_AVG_EFF", "GorgASPECT",

spatial_vars <- c("GRD_LON", "GRD_LAT")

n_vars <- length(env_vars)

gf_trees <- 50
gf_bins <- 201

n_vars <- length(env_vars)


extrapolate_gf <- "compress"
extrap_pow <- 1/32


dataset_names <- c( "Trawl"
                   )
gf_models <- foreach(dataset = dataset_names, .packages = c("gradientForest"), .combine = c) %do% {
  ret <- list()
  ret[[dataset]] <- bootstrapGradientForest(cbind(datasets[[dataset]][["X"]],
                                                  datasets[[dataset]][["Y"]]),
                                            predictor.vars = env_vars,
                                            response.vars = names(datasets[[dataset]][["Y"]]),
                                            nbootstrap=gf_trees,
                                            compact = F, nbin = gf_bins,
                                            transform = NULL,
                                            corr.threshold = 0.5,
                                            maxLevel = floor(log2(dim(datasets[[dataset]][["Y"]])[1]*0.368/2)),
                                            trace = TRUE
                                            )
  return(ret)
}

##don't have a method for combining bootstraps yet
## gf_combined <- do.call("combinedGradientForest", c(gf_models, list(standardize="after")))

## TODO: filter out variables that don't get used
## TODO: figure out why many variables are creating a very poor fit
## TODO: see if offsets need to use model weights
## DONE: restrict env to 1.1*max and 0.9*min
##   Works best using quantiles, I took 10% and 90%
## TODO: 14,000 points may crash my laptop
##   Restricting to quantiles gives only 5000 points.
##   I also reduced it to 2000
## Some variables seem to have a single split
##using dataset with many sites
gf_trawl <- gf_models$Trawl

gf_plots <- gg_bootstrapGF(gf_trawl, n_curves = 100, debug = TRUE)
lapply(names(gf_plots), function(pl, gg_boot) {
  ggsave(filename = paste0("pilbara_", pl, "_curves.png"), plot = gg_boot[[pl]],
         units = "cm",
         width = 16,
         height = 9,
         dpi = 300,
         scale = 2
         )
}, gg_boot = gf_plots)

imp_vars <- pred_names(gf_trawl)
Phys_grid_complete <- Phys_grid[complete.cases(Phys_grid[, env_vars]), ]
Phys_grid_use <- Phys_grid_complete

grid_range <- apply(Phys_grid_use[,env_vars], 2, range)
sample_range <- apply(datasets$Trawl$X[,env_vars], 2, range)
## limits <- rbind(sample_range[1,]*0.9, sample_range[2, ]*1.1 )

quantile(datasets$Trawl$X[,"CRS_NO3_AV"], probs = c(0.1, 0.9))

limits <- apply(datasets$Trawl$X[,env_vars], 2,
                function(x){
                  quantile(x, probs = c(0.1,0.9))
                })

phys_cropped <- do.call("cbind", lapply(env_vars,
                             function(env_v,
                                      Phys_grid_use,
                                      limits){
                               out <- Phys_grid_use[, env_v]
                               out_na <- out > limits[2, env_v] | out < limits[1, env_v]
                               print(sum(out_na))
                               out[out_na] <- NA
                               out_df <- data.frame(out)
                               names(out_df) <- env_v
                               return(out_df)
                             },
                             limits = limits,
                             Phys_grid_use = Phys_grid_use)
                      )

Phys_grid_x_row <- cbind(Phys_grid_complete, x_row = seq_along(Phys_grid_complete$GRID_ID))

phys_cropped_spa <- cbind(Phys_grid_x_row[,c("x_row", "GRD_LON", "GRD_LAT" )] , phys_cropped )

phys_cropped_spa_cut <- phys_cropped_spa[complete.cases(phys_cropped_spa), ]
subset_rows <- sample.int(n = nrow(phys_cropped_spa_cut), size = nrow(phys_cropped_spa_cut))

##issue with x_rows
##x_rows are for Phys_grid before any clipping
##predict.bootstrapGradientForest is unaware of them
## and just gives 1:nrow of input
##Map grid cells to predict cells.
## Need to create new x_row that aligns with rows of
## phys_crop_spa_cut
phys_predict <- predict.bootstrapGradientForest(gf_trawl,
                                ## Phys_grid_use[1:14000, imp_vars],
                                                phys_cropped_spa_cut[subset_rows, imp_vars],
                                type = c("mean", "variance", "points"),
                                extrap=extrapolate_gf,
                                extrap_pow = extrap_pow
                                )
phys_cropped_spa_cut_x_row <- phys_cropped_spa_cut[subset_rows, ]
phys_cropped_spa_cut_x_row$x_row <- seq_along(phys_cropped_spa_cut_x_row$x_row)
phys_pred_spat <- merge(phys_predict[!is.na(phys_predict$y), ], phys_cropped_spa_cut_x_row[,c("x_row", "GRD_LAT", "GRD_LON")])

## e_vars <- names(importance(gf_combined))[1:n_vars]
##plot the predictions
ggsave("pilbara_boot_variance.png", ggplot2::ggplot(data = phys_pred_spat[phys_pred_spat$type == "variance",], mapping =
                                                       aes(x = x, y = y)) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
  units = "cm",
  width = 16,
  height = 9,
  dpi = 300,
  scale = 2
  )

ggsave("pilbara_boot_data.png", ggplot2::ggplot(data = phys_pred_spat[phys_pred_spat$type == "variance",], mapping =
                                                                                                                 aes(x = x)) +
                                    ggplot2::geom_histogram() +
                                ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
       )

ggsave("pilbara_boot_mean.png", ggplot2::ggplot(data = phys_pred_spat[phys_pred_spat$type == "mean",], mapping =
                                                       aes(x = x, y = y)) +
         ggplot2::geom_line() +
         ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
)
pred_mean_var <- merge(phys_pred_spat[phys_pred_spat$type == "mean", c("var",  "x", "y" )],
                       phys_pred_spat[phys_pred_spat$type == "variance",c("var",  "x", "y")],
                       by = c("var", "x"))
names(pred_mean_var)[3:4] <- c("mean", "variance")
pred_mean_var["upper"] <-pred_mean_var$mean + sqrt(pred_mean_var$variance)
pred_mean_var["lower"] <-pred_mean_var$mean - sqrt(pred_mean_var$variance)
plot_var <-  "CRS_NO3_AV"
ggsave("pilbara_boot_mean_variance.png", ggplot2::ggplot(data =pred_mean_var[pred_mean_var$var ==plot_var,], mapping =
                                                            aes(x = x, y = mean)) +
         ggplot2::geom_line() +
         ggplot2::geom_line(data = pred_mean_var[pred_mean_var$var == plot_var,c("x", "upper")], mapping = aes(x = x, y = upper),
                              colour = "red") +
         ggplot2::geom_line(data = pred_mean_var[pred_mean_var$var == plot_var,c("x", "lower")], mapping = aes(x = x, y = lower),
                            colour = "red"),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
)

ggsave("pilbara_boot_points.png", ggplot2::ggplot(data = phys_pred_spat[phys_pred_spat$type == "points",], mapping =
                                                       aes(x = x, y = y, group = gf_model, color = as.factor(gf_model) )) +
         ggplot2::geom_line(size = 0.1) +
         ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
)
ggsave("pilbara_boot_points_quantile.png", ggplot2::ggplot(data = phys_pred_spat[phys_pred_spat$type == "points",], mapping =
                                                     aes(x = x, y = y )) +
         ggplot2::geom_quantile(quantiles = c(0.05, 0.5, 0.95), method = "rqss", lambda = 0.1) +
         ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
)
ggsave("pilbara_boot_points_smooth.png", ggplot2::ggplot(data = phys_pred_spat[phys_pred_spat$type == "points",], mapping =
                                                              aes(x = x, y = y )) +
         ggplot2::geom_smooth() +
         ggplot2::facet_wrap(ggplot2::vars(var), scales = "free_x"),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
)

##cluster

library(rmethods)
##need to make a wide form matrix for hotellings_bulk
library(tidyr) #I am not using this within my package, so it is not a dependency.
predict_wide <- tidyr::pivot_wider(data = phys_pred_spat[phys_pred_spat$type %in% c("mean", "variance"),],
                                         names_from =  c("var"),
                                         values_from = c("y"),
                                         id_cols = c("x_row", "type", "GRD_LAT", "GRD_LON")  )
pred_vars <- pred_names(gf_trawl)
predict_wide_means <- predict_wide[predict_wide$type == "mean",]
predict_wide_variance <- predict_wide[predict_wide$type == "variance",]
p_mat <- rmethods:::hotellings_bulk(
             means = predict_wide_means[, pred_vars],
             res_sq = (predict_wide_variance[, pred_vars])^1
           )

##p_mat is 1 for similar, 0 for different.
##p_mat <- 1-p_mat

p_mat_long <- tidyr::pivot_longer(data = data.frame(x = as.numeric(row.names(p_mat)), p_mat),
                                  cols = -x,
                                  names_to = c("y"),
                                  names_prefix = "X",
                    values_to = c("p"))
p_mat_long$y <- as.integer(p_mat_long$y)

##This is not wise, plotting a 2000x2000 grid
##it works though
## ggplot2::ggsave("pilbara_p_mat.png",
##                 ggplot2::ggplot(data = p_mat_long,
##                                 mapping = aes(x = x, y = y, fill = p)) +
##                 ggplot2::geom_raster(),
##                 units = "cm",
##                 width = 16,
##                 height = 9,
##                 dpi = 300,
##                 scale = 2
##                 )


##Bring in WGCNA!

library(WGCNA)

##following the WGCNA tutorials

powers = c(c(1:10), seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(p_mat, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
## this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
## Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##My results show a knee at 16
pow <- 16
net = blockwiseModules(p_mat,
                       power = pow,
                       TOMType = "unsigned",
                       minModuleSize = 10, 
                       reassignThreshold = 0,
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "pilbara_data",
                       verbose = 3)

table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

nGenes = ncol(p_mat)
nSamples = nrow(p_mat)

dissTOM = 1-TOMsimilarityFromExpr(p_mat, power = pow)
plotTOM = dissTOM^7
diag(plotTOM) = NA
png(filename = "pilbara_heatmap.png",
units = "cm",
width = 16,
height = 9,
res = 300,
)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

##Now plot in space
clust_df <- data.frame(x_row = seq_along(net$colors), clust = net$colors)
clust_spat <- unique(merge(x = phys_pred_spat[,c("x_row", "GRD_LON", "GRD_LAT")], y = clust_df, by = "x_row" ))

ggsave("pilbara_clust.png", ggplot2::ggplot(data = clust_spat,
                                            mapping =
                                              ggplot2::aes(x = GRD_LON, y = GRD_LAT, fill = as.factor(clust)) )+
                                  ggplot2::geom_raster(mapping = aes(fill = as.factor(clust))),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
       )

##pilbara with corrected WGCNA usage

##My block size is ~5000, I will do all at once

sft = pickSoftThreshold.fromSimilarity(p_mat, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
## this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
## Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##here pow is much lower, 3 gives a peak
pow <- 3
adj_mat <- adjacency.fromSimilarity(p_mat, type = "unsigned", power = pow)

##TOMDenom "mean" is experimental
tom <- TOMsimilarity(adj_mat, TOMType = "unsigned", TOMDenom = "mean")

diss_tom = 1 - tom

dendros = fastcluster::hclust(as.dist(diss_tom), method = "average")

##I took this code from the source of WGCNA, so I need to provide values 
blockLabels = try(cutreeDynamic(dendro = dendros, 
                                deepSplit = 2, #default
                                cutHeight = 0.995, #default
                                minClusterSize = min(20, ncol(p_mat)/2), #default
                                method = "hybrid", #hardcode
                                distM = diss_tom,
                                maxCoreScatter = NULL,#default 
                                minGap = NULL,#default
                                maxAbsCoreScatter = NULL,#default
                                minAbsGap = NULL,#default
                                minSplitHeight = NULL,#default
                                minAbsSplitHeight = NULL,#default
                                externalBranchSplitFnc = list(),#default
                                minExternalSplit = numeric(0),#default
                                externalSplitOptions = list(),#default
                                externalSplitFncNeedsDistance = logical(0), #default
                                assumeSimpleExternalSpecification = FALSE,#hardcode
#
                                pamStage = TRUE, #default
                                pamRespectsDendro = FALSE,
                                verbose = 3-3,
                                indent = 0 + 2),
                  silent = FALSE);


clust_df <- data.frame(x_row = seq_along(blockLabels), clust = blockLabels)
clust_spat <- unique(merge(x = phys_pred_spat[,c("x_row", "GRD_LON", "GRD_LAT")], y = clust_df, by = "x_row" ))

ggsave("pilbara_clust.png", ggplot2::ggplot(data = clust_spat,
                                            mapping =
                                              ggplot2::aes(x = GRD_LON, y = GRD_LAT, fill = as.factor(clust)) )+
                            ggplot2::geom_raster(mapping = aes(fill = as.factor(clust))),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
       )
##

aff_thres <- 0.95
cast_obj <- cast_alg(p_mat, aff_thres)

stable_cast <- cast_stabilize(cast_obj, aff_thres, p_mat, max_iter = 8)

##Row by cluster
clustering <- do.call("rbind",
                      lapply(seq_along(stable_cast), function(i, stable_cast) {
                        d <- which(stable_cast[[i]])
                        data.frame(x_row = d, clust = i)
                      }, stable_cast = stable_cast)
                      )

##plot clustering
clust_spat <- unique(merge(x = phys_pred_spat[,c("x_row", "GRD_LON", "GRD_LAT")], y = clustering, by = "x_row" ))

##reorder p_mat by clusters
clust_reorder <- do.call("c", lapply(stable_cast, which))
p_mat_cast <- p_mat[clust_reorder,clust_reorder]
colnames(p_mat_cast) <- 1:nrow(p_mat_cast)
p_mat_long_cast <- tidyr::pivot_longer(data = data.frame(x = 1:nrow(p_mat_cast), p_mat_cast),
                                  cols = -x,
                                  names_to = c("y"),
                                  names_prefix = "X",
                                  values_to = c("p"))
p_mat_long_cast$y <- as.integer(p_mat_long_cast$y)
ggplot2::ggsave("pilbara_p_mat_cast.png",
                ggplot2::ggplot(data = p_mat_long_cast,
                                mapping = aes(x = x, y = y, fill = p)) +
                ggplot2::geom_raster(),
                units = "cm",
                width = 16,
                height = 9,
                dpi = 300,
                scale = 2
                )

##before the stability stage
clust_reorder <- do.call("c", lapply(cast_obj, which))
p_mat_cast <- p_mat[clust_reorder,clust_reorder]
colnames(p_mat_cast) <- 1:nrow(p_mat_cast)
p_mat_long_cast <- tidyr::pivot_longer(data = data.frame(x = 1:nrow(p_mat_cast), p_mat_cast),
                                       cols = -x,
                                       names_to = c("y"),
                                       names_prefix = "X",
                                       values_to = c("p"))
p_mat_long_cast$y <- as.integer(p_mat_long_cast$y)
ggplot2::ggsave("pilbara_p_mat_cast_unstable.png",
                ggplot2::ggplot(data = p_mat_long_cast,
                                mapping = aes(x = x, y = y, fill = p)) +
                ggplot2::geom_raster(),
                units = "cm",
                width = 16,
                height = 9,
                dpi = 300,
                scale = 2
                )


ggsave("pilbara_clust_cast.png", ggplot2::ggplot(data = clust_spat,
                                            mapping =
                                              ggplot2::aes(x = GRD_LON, y = GRD_LAT, fill = as.factor(clust)) )+
                            ggplot2::geom_raster(mapping = aes(fill = as.factor(clust))),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
       )

corner_size <- 500
p_mat_long_corner <- tidyr::pivot_longer(data = data.frame(x = as.numeric(row.names(p_mat[1:corner_size,1:corner_size])), p_mat[1:corner_size,1:corner_size]),
                                  cols = -x,
                                  names_to = c("y"),
                                  names_prefix = "X",
                                  values_to = c("p"))
p_mat_long_corner$y <- as.integer(p_mat_long_corner$y)
names(p_mat_long_corner) <- c("Source", "Target", "edge")
write_delim(p_mat_long_corner, delim = ";", path = "p_mat_corner.csv")



##see which sites appear in more than 1 cluster
library(tidyverse)
assignments <- predict_clust(cast_obj, p_mat)
assignments <- assignments %>%
  filter(aff > aff_thres) %>%
  count(x_row)
##merge and plot by "ecotone"
ecotones <- merge(assignments, phys_pred_spat[,c("x_row", "GRD_LON", "GRD_LAT")], by =  "x_row")
ggsave("pilbara_ecotones.png", ggplot2::ggplot(data = ecotones,
                                            mapping =
                                              ggplot2::aes(x = GRD_LON, y = GRD_LAT, fill = as.factor(n)) )+
                            ggplot2::geom_raster(mapping = aes(fill = as.factor(n))),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
       )

## dendrogram of clusters, according to the between cluster affinity, to reorder
## the clusters for plotting
aff_bt <- aff_cluster_between(cast_obj, p_mat)
ggsave("pilbara_clust_between.png", ggplot2::ggplot(data = aff_bt,
                                                    mapping =
                                                      ggplot2::aes(x = Var1, y = Var2, fill = affs) )+
                                    ggplot2::geom_raster(mapping = aes(fill = affs)),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
       )
aff_bt_stable <- aff_cluster_between(stable_cast, p_mat)
ggsave("pilbara_clust_between_stable.png", ggplot2::ggplot(data = aff_bt_stable,
                                                 mapping =
                                                   ggplot2::aes(x = Var1, y = Var2, fill = affs) )+
                                 ggplot2::geom_raster(mapping = aes(fill = affs)),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
       )

aff_bt[aff_bt$Var1 == aff_bt$Var2, "affs"]
do.call("c", aff_wt)


aff_wt <- aff_clust_inner(cast_obj, p_mat)
aff_wt_stable <- aff_clust_inner(stable_cast, p_mat)

library(fastcluster)
aff_bt_square <- as.dist(as.matrix(tidyr::pivot_wider(aff_bt, id_cols = "Var2", values_from = "affs" ,names_from = "Var1" )[,-1]))
dendro <- fastcluster::hclust(aff_bt_square)

n_clust <- length(cast_obj)
##a bit complicated
## I want cluster 1 to become cluster 3
reorder_frame <- data.frame(old = 1:n_clust, new = dendro$order)
tmp <- merge(aff_bt, reorder_frame, by.x = "Var2", by.y = "old")
names(tmp)[names(tmp) == "new"] <- "Var2_new"
aff_bt_reorder <- merge(tmp, reorder_frame, by.x = "Var1", by.y = "old", )
names(aff_bt_reorder)[names(aff_bt_reorder) == "new"] <- "Var1_new"
ggsave("pilbara_clust_between_dendro.png", ggplot2::ggplot(data = aff_bt_reorder,
                                                    mapping =
                                                      ggplot2::aes(x = Var1_new, y = Var2_new, fill = affs) )+
                                    ggplot2::geom_raster(mapping = aes(fill = affs)),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
       )
##weak, trying a different method
##level 2 heirarchy
l2_mat <- matrix(aff_bt$affs, nrow = n_clust, ncol = n_clust)
l2_clust <- cast_alg(l2_mat, 0.80)

l2_clust_reorder <- do.call("c", lapply(l2_clust, which))
l2_reorder <- l2_mat[l2_clust_reorder,l2_clust_reorder]
colnames(l2_reorder) <- 1:nrow(l2_reorder)
l2_reorder_long <- tidyr::pivot_longer(data = data.frame(x = 1:nrow(l2_reorder), l2_reorder),
                                       cols = -x,
                                       names_to = c("y"),
                                       names_prefix = "X",
                                       values_to = c("p"))
l2_reorder_long$y <- as.integer(l2_reorder_long$y)
ggplot2::ggsave("pilbara_l2_reorder_clusters.png",
                ggplot2::ggplot(data = l2_reorder_long,
                                mapping = aes(x = x, y = y, fill = p)) +
                ggplot2::geom_raster(),
                units = "cm",
                width = 16,
                height = 9,
                dpi = 300,
                scale = 2
                )


heir_clust <- cast_h(sim_mat = p_mat, seq(0.5, 0.95, 0.05))
dep <- 8
p_reorder <- cast_h_reorder(heir_clust, depth = dep)
aff_l <- heir_clust[[dep+1]]$affs
aff_l_cast <- aff_l[p_reorder,p_reorder]
colnames(aff_l_cast) <- 1:nrow(aff_l_cast)
aff_l_long_cast <- tidyr::pivot_longer(data = data.frame(x = 1:nrow(aff_l_cast), aff_l_cast),
                                       cols = -x,
                                       names_to = c("y"),
                                       names_prefix = "X",
                                       values_to = c("p"))
aff_l_long_cast$y <- as.integer(aff_l_long_cast$y)
ggplot2::ggsave(sprintf("pilbara_aff_d_%i_cast_unstable_dendro.png", dep),
                ggplot2::ggplot(data = aff_l_long_cast,
                                mapping = aes(x = x, y = y, fill = p)) +
                ggplot2::geom_raster(),
                units = "cm",
                width = 16,
                height = 9,
                dpi = 300,
                scale = 2
                )

dep <- 10
p_reorder <- cast_h_reorder(heir_clust, depth = dep)
aff_l_cast <- p_mat[p_reorder,p_reorder]
colnames(aff_l_cast) <- 1:nrow(aff_l_cast)
aff_l_long_cast <- tidyr::pivot_longer(data = data.frame(x = 1:nrow(aff_l_cast), aff_l_cast),
                                       cols = -x,
                                       names_to = c("y"),
                                       names_prefix = "X",
                                       values_to = c("p"))
aff_l_long_cast$y <- as.integer(aff_l_long_cast$y)
ggplot2::ggsave("pilbara_aff_full_cast_unstable_dendro.png",
                ggplot2::ggplot(data = aff_l_long_cast,
                                mapping = aes(x = x, y = y, fill = p)) +
                ggplot2::geom_raster(),
                units = "cm",
                width = 16,
                height = 9,
                dpi = 300,
                scale = 2
                )

##Again, but with stabilization steps
heir_clust <- cast_h(sim_mat = p_mat, seq(0.5, 0.95, 0.05), stabilize = TRUE, max_stab_iter = 5)
dep <- 8
p_reorder <- cast_h_reorder(heir_clust, depth = dep)
aff_l <- heir_clust[[dep+1]]$affs
aff_l_cast <- aff_l[p_reorder,p_reorder]
colnames(aff_l_cast) <- 1:nrow(aff_l_cast)
aff_l_long_cast <- tidyr::pivot_longer(data = data.frame(x = 1:nrow(aff_l_cast), aff_l_cast),
                                       cols = -x,
                                       names_to = c("y"),
                                       names_prefix = "X",
                                       values_to = c("p"))
aff_l_long_cast$y <- as.integer(aff_l_long_cast$y)
ggplot2::ggsave(sprintf("pilbara_aff_d_%i_cast_stable_dendro.png", dep),
                ggplot2::ggplot(data = aff_l_long_cast,
                                mapping = aes(x = x, y = y, fill = p)) +
                ggplot2::geom_raster(),
                units = "cm",
                width = 16,
                height = 9,
                dpi = 300,
                scale = 2
                )

dep <- 10
p_reorder <- cast_h_reorder(heir_clust, depth = dep)
aff_l_cast <- p_mat[p_reorder,p_reorder]
colnames(aff_l_cast) <- 1:nrow(aff_l_cast)
aff_l_long_cast <- tidyr::pivot_longer(data = data.frame(x = 1:nrow(aff_l_cast), aff_l_cast),
                                       cols = -x,
                                       names_to = c("y"),
                                       names_prefix = "X",
                                       values_to = c("p"))
aff_l_long_cast$y <- as.integer(aff_l_long_cast$y)
ggplot2::ggsave("pilbara_aff_full_cast_stable_dendro.png",
                ggplot2::ggplot(data = aff_l_long_cast,
                                mapping = aes(x = x, y = y, fill = p)) +
                ggplot2::geom_raster(),
                units = "cm",
                width = 16,
                height = 9,
                dpi = 300,
                scale = 2
                )

##Scree plots

ggsave("pilbara_scree_stable.png",
       ggplot(data.frame(y = sort(sapply(stable_cast, sum), decreasing = TRUE), x = seq_along(cast_obj)), aes(x=x, y=y)) +
       geom_point(),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
       )
ggsave("pilbara_scree_unstable.png",
       ggplot(data.frame(y = sort(sapply(cast_obj, sum), decreasing = TRUE), x = seq_along(cast_obj)), aes(x=x, y=y)) +
  geom_point(),
  units = "cm",
  width = 16,
  height = 9,
  dpi = 300,
  scale = 2
  )

##before the stability stage
clust_reorder <- do.call("c", lapply(cast_obj[dendro$order], which))
p_mat_cast <- p_mat[clust_reorder,clust_reorder]
colnames(p_mat_cast) <- 1:nrow(p_mat_cast)
p_mat_long_cast <- tidyr::pivot_longer(data = data.frame(x = 1:nrow(p_mat_cast), p_mat_cast),
                                       cols = -x,
                                       names_to = c("y"),
                                       names_prefix = "X",
                                       values_to = c("p"))
p_mat_long_cast$y <- as.integer(p_mat_long_cast$y)
ggplot2::ggsave("pilbara_p_mat_cast_unstable_dendro.png",
                ggplot2::ggplot(data = p_mat_long_cast,
                                mapping = aes(x = x, y = y, fill = p)) +
                ggplot2::geom_raster(),
                units = "cm",
                width = 16,
                height = 9,
                dpi = 300,
                scale = 2
                )



##k-means, with threshold by p-values
##I will fit a k_medoid, then look at the largest p value in the group
k_range <- 2:30

library(cluster)
phys_pred_mean <- phys_pred_spat %>%
  filter(type == "mean") %>%
  select(x_row, var, y ) %>%
  tidyr::pivot_wider(names_from = c("var"), values_from = c("y"))

k_clusts <- lapply(k_range, function(k, phys_grid, sim_mat){
  clust <- clara(x = phys_pred_mean[ , env_vars], k = k, samples = 100, sampsize = 500, rngR = TRUE, pamLike = TRUE)
  p_vals <- sim_mat[clust$i.med, clust$i.med]
  p_limits <- range(p_vals[upper.tri(p_vals)])
  return(list(clust = clust, p_vals = p_vals, p_limits = p_limits))
}, phys_grid = phys_pred_mean, sim_mat = p_mat)

p_range_set <- do.call("rbind", lapply(k_clusts, function(k_c){
  k <- length(k_c$clust$i.med)
  return(data.frame(k = k, type = c("min", "max"), p = k_c$p_limits))
})
)

ggsave("pilbara_k_medoid_variance.png",
       ggplot(p_range_set, aes(x= k, y = p, group = type, colour = as.factor(type))) +
       geom_line(),
       units = "cm",
       width = 16,
       height = 9,
       dpi = 300,
       scale = 2
       )


}
