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

#' ggplot object of the similarity matrix
#'
#' Without a cast object, it will just plot them
#' similarity matrix as is, with no particular order.
#'
#' With a cast object, the similarity matrix will
#' be reordered to put sites in the same cluster
#' together.
#'
#' If secondary_sort is true, then the clusters will
#' be further reordered to put similar clusters together,
#' using a dendrogram sort.
#'
#' if highlight = TRUE, then rectangles will be drawn
#' around each cluster.
gg_sim_mat <- function(sim_mat,
                       cast_ob = NULL,
                       sort_between = TRUE,
                       sort_within = TRUE,
                       highlight = FALSE
                       ){

  if(!is.null(cast_ob)) {
    ##order by cast object

    if(sort_within){
      cast_ob <- lapply(cast_ob, function(clust, sim_mat){
        ##sort by strength of affinity to cluster
        order(rowMeans(sim_mat[clust, clust, drop = FALSE]))
        return(clust[order])
      },  sim_mat = sim_mat)
    }
    if(sort_between) {
      aff_btw <- aff_cluster_between(sim_mat = sim_mat, cast_obj = cast_ob)

      aff_btw_wide <- matrix(aff_btw$affs, sqrt(nrow(aff_btw)), sqrt(nrow(aff_btw)))

      aff_dist <- dist(aff_btw_wide)
      aff_sort <- hclust(aff_dist)

      clust_reorder <- aff_sort$order

    } else {
      clust_reorder <- seq(1, length(cast_ob))
    }
    site_reorder <- do.call("c", cast_ob[clust_reorder])

    ##add cluster rectangles
    if(highlight) {
      rects <- do.call("rbind", lapply(seq(1, length(cast_ob)), function(i, cast_ob_i){
        if(i > 1){
          start <- do.call(sum, lapply(1:(i-1), function(j, cast_ob_j){
            length(cast_ob_j[[j]])
          }, cast_ob_j = cast_ob_i)) + 1
        } else{
          start <- 1
        }
        end <- start + length(cast_ob_i[[i]])
        return(data.frame(xmin = start, xmax = end))
      }, cast_ob_i = cast_ob[clust_reorder]))
    }

  } else {
    site_reorder <-  seq.int(1, nrow(sim_mat))
  }
  sim_mat_reorder <- sim_mat[site_reorder, site_reorder]
  sim_mat_long <- data.frame(expand.grid(1:nrow(sim_mat), 1:ncol(sim_mat)), as.vector(as.matrix(sim_mat_reorder)))
  names(sim_mat_long) <- c("x", "y", "p")
  p <- ggplot2::ggplot(data = sim_mat_long,
                  mapping = aes(x = x, y = y, fill = p)) +
    ggplot2::geom_raster()
  if(highlight){
    p <- p + ggplot2::annotate(geom = "rect", xmin = rects$xmin, ymin = rects$xmin, xmax = rects$xmax, ymax = rects$xmax,  colour = "red", fill = NA)
  }
  return(p)
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
    gf_names <- paste("F",1:ngear,sep="")
  if (any(empty <- gf_names==""))
    gf_names[empty] <- paste("F",1:ngear,sep="")[empty]

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

  gf_combine_dist <- combine_gfbootstrap_dist(gf_combine, x_samples = x_samples)

  gf_combine$offsets <- gfbootstrap_offsets(gf_combine_dist)

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
#' 
pred_names_combined <- function(obj) {
  unique(do.call("c", future.apply::future_lapply(obj$gf_list,
                                                  function(x){
                                                    ##now each is a combined GF
                                                    names(x$CU)
                                                  })
                 ))
}
                                        #             (2) List it as "suggests" in your DESCRIPTION file and
                                        #precede each use with "if(require(desiredR_ForgePackage))".  If
                                        #"require" returns TRUE, you do what you want.  Else issue an error message.

#' Combined gf_dist
#'
#' Calculate the distance between bootstrapped
#' samples of combined gradientForest models.
#'
#'
combine_gfbootstrap_dist <- function(gf_combine,
                                     x_samples = 100){
  assertthat::assert_that(inherits(gf_combine, "combinedBootstrapGF"))

  pred_vars <- pred_names_combined(gf_combine) 
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
#' Returns a ggplot object, ready for plotting
gg_combined_bootstrapGF <- function(x,
                           vars = pred_names_combined(x),#names(importance(x$gf_list[[1]], type = "Weighted", sorted = TRUE)),
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
      
      if (is.element(pred, names(x$gf_list[[gf]]$CU)) ){
        curve_data <- gradientForest::cumimp(x$gf_list[[gf]], pred, weight = x$weight)
      } else {
        message(paste0("gg_combined_bootstrapgf: Tree [", gf,
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

#' Predict combinedBootstrapGF
#'
#' Generic S3 function for predicting combinedBootstrapGF
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
#' extrap can take on 2 values, c(TRUE, FALSE)
#'
#' extrap = TRUE extrapolates, but applies compression so predicted
#' values approach an asymptote as X moves very far from training values.
#' extrap_pow is only used with extrap = TRUE, and sets the
#' compression power. 1/4 gives the 4th root, 1 is equivalent to linear extrapolation
#' and 0 is equivalent to clipping to the maximum
#'
#'
#' extrap = FALSE return Na outside of the training values. extrap = FALSE only uses
#' models that observed the training data when calculating mean and variance
#'
#' ... arguments passed to cumimp
#' If `type` contains more than one string, a list of
#' data.frames will be returned, otherwise just a data.frame.
#' 
#' TODO: This code is rather spaghetti, eg. activity 4 depends on 2 and 3,
#' then activity 5 depends on 1 and 4
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
    newdata <- object$gf_list[[1]]$X[,pred_names_combined(object)]
  } else {
    newnames <- names(newdata)
    gf_preds <- pred_names_combined(object)
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

aff_func <- function(x, u, sim_mat){
  ##Takes the vector x of rows and the row u of sim_match
  ##and returns the affinity from u to x
  assertthat::assert_that(length(x) == nrow(sim_mat))
  assertthat::assert_that(is.logical(x))
  assertthat::assert_that(is.numeric(u))

  ##Simple, this is just elements of sim_mat
  return(sim_mat[u, ]*x)
}

aff_clust_all <- function(sim_mat, cast_ob){
  do.call(cbind, lapply(1:length(cast_ob), 
                        function(clust, sim_mat, cast_ob){

                          ##get the affinity to each cluster
                          if(length(cast_ob[[clust]]) > 0){
                            return(rowMeans(sim_mat[, cast_ob[[clust]], drop = FALSE]))
                          } else {
                            ##empty clusters have 0 affinity
                            return(rep(0, nrow(sim_mat)))
                          }
                        }, cast_ob = cast_ob, sim_mat = sim_mat))
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
    valid_seeds <- spares

    aff_local <- rep(0, n)

    ##Do the next steps until stability is reached
    repeat{
      is_change <- FALSE



      ##find initial cluster point
      if(all(!new_clust)){
        maxima <- which(spares & valid_seeds)
        new_elem <- maxima[sample.int(length(maxima), 1)] ##intention: take one of the max at random
        next_seed <- new_elem
        new_clust[new_elem] <- TRUE
        spares[new_elem] <- FALSE
        is_change <- TRUE
        debug_spares_expected <- debug_spares_expected - 1
        ## print(list(new_elem, sum(spares), debug_spares_expected, "first"))
        assertthat::assert_that(debug_spares_expected == sum(spares))
        ##update aff_local
        ## aff_local[new_clust | spares] <- aff_local[new_clust | spares] + aff_func(new_clust | spares, new_elem, sim_mat)
        aff_local <- aff_local + aff_func(new_clust | spares, new_elem, sim_mat)
      }

      ##addition stage
      ##affinity seems to be used to bias group size
      ##Strong affinity is needed to join a large group

      ##Short-circuit evaluation, won't try to find max of an empty set
      ##trialling a slightly different approach.
      ##Find any high affinity elements.
      ##high affinity, sum of affinities into cluster (which is aff_local)
      ##exceeds thres * size of cluster
      ##So now it is a rolling mean. At each iteration, you add to the sum,
      ## then only divide by n when you need to test.
      while(any(spares) && (max(aff_local[spares]/sum(new_clust)) >= aff_thres)){
        maxima <- which(aff_local == max(aff_local[spares]) & spares)
        new_elem <- maxima[sample.int(length(maxima), 1)] ##intention: take one of the max at random
        new_clust[new_elem] <- TRUE
        spares[new_elem] <- FALSE
        debug_spares_expected <- debug_spares_expected - 1
        ## print(list(new_elem, sum(spares), debug_spares_expected, "add"))
        assertthat::assert_that(debug_spares_expected == sum(spares))
        is_change <- TRUE
        ##update aff_local
        ## aff_local[new_clust | spares] <- aff_local[new_clust | spares] + aff_func(new_clust | spares, new_elem, sim_mat)
        aff_local <- aff_local + aff_func(new_clust | spares, new_elem, sim_mat)
      }

      ##Removal stage
      while(any(new_clust) && (min(aff_local[new_clust]/sum(new_clust)) < aff_thres)){
        minima <- which(aff_local == min(aff_local[new_clust]) & new_clust)
        new_elem <- minima[sample.int(length(minima), 1)] ##intention: take one of the max at random
        new_clust[new_elem] <- FALSE
        spares[new_elem] <- TRUE
        debug_spares_expected <- debug_spares_expected + 1
        ## print(list(new_elem, sum(spares), debug_spares_expected, "remove"))
        assertthat::assert_that(debug_spares_expected == sum(spares))
        is_change <- TRUE
        ##update aff_local
        aff_local <- aff_local - aff_func(new_clust | spares, new_elem, sim_mat)
      }

      if(all(!new_clust)){
        ##cluster is empty
        valid_seeds[next_seed] <- FALSE
        message("seeds left: [", sum(valid_seeds), "]")
        if(all(!valid_seeds)){
          ##no more valid seeds exist, all have been tried, create a leftovers cluster
          new_clust <- spares
          spares <- rep(FALSE, n)
          break
        }
      }

      if(!is_change){
        break
      }
    }
    message("cluster assigned of size [", sum(new_clust), "]")
    message("[", sum(spares), "] sites left to assign.")
    assertthat::assert_that(debug_spares_expected == (debug_spares_loop - sum(new_clust)))
    clust[[clust_id]] <- which(new_clust)
    clust_id <- clust_id + 1


  }

  return(clust)
}



#' Stabilize cluster membership
#'
#' cast_stabilize() works iteratively, updating one site at a time
#' and recalculating the affinities each time.
#' cast_stabilize() will only make an update if it either maintains
#' all the average within cluster affinities above the threshold, or
#' increases the within cluster affinities if one is already below the
#' threshold.
#'
#' cast_stabilize_batch() does a bulk update, finding all sites that
#' could be assigned to another cluster before recalulating the affinities.
#' cast_stabilize_batch() always proceeds without checking the impact
#' the update has on within-cluster affinities.
cast_stabilize_batch <- function(cast_obj, aff_thres, sim_mat, max_iter = 20){
  iter <- 1

  while(iter <= max_iter){
    ##For each vertex, find affinity to other clusters
    ##Cast_obj is not huge, but constantly recreating it when I am just
    ##flipping bools seems wasteful
    ##This approach tests each site, and updates clustering at end
    updates <- lapply(seq_along(cast_obj[[1]]), function(u, cast_obj, sim_mat){
      clust_id <- which(sapply(cast_obj, function(clust, u){
        u %in% clust
      }, u = u))
      assertthat::assert_that(length(clust_id) == 1)
      ##u belongs to clust_id
      clust_aff <- lapply(cast_obj, function(clust, u, sim_mat){
          ##get the affinity to each cluster
          if(any(clust)){
            return(mean(sim_mat[u, clust]))
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
        cast_obj[[ updates[upd,"old"] ]] <- cast_obj[[ updates[upd,"old"] ]][cast_obj[[ updates[upd,"old"] ]] != updates[upd,"u"]  ] 
        cast_obj[[updates[upd,"new"]]] <- c(cast_obj[[updates[upd,"new"]]], updates[upd,"u"])
      }
    } else {
      break
    }
    iter <- iter + 1
  }
  return(cast_obj)
}
cast_stabilize <- function(cast_obj, aff_thres, sim_mat, max_iter = nrow(sim_mat)*2){
  iter <- 1

  #matrix, mean affinity of each site to each cluster
  ##two cols must be updated each iteration
  clust_aff <- aff_clust_all(sim_mat, cast_obj)
  ##maps element to cluster
  ##can be updated with a single value change each loop
  clust_ind <- lapply(seq_along(cast_obj), function(clust, cast_obj){
    data.frame(elem = cast_obj[[clust]], clust = clust)
  }, cast_obj = cast_obj)
  clust_ind <- do.call(rbind, clust_ind)
  clust_ind <- clust_ind[order(clust_ind$elem), ]

  ##sites that should not be moved, to preserve aff_thres
  locked_sites <- c()

  while(iter <= max_iter){
    ##single step updates require some conserved data.

    ##find maximal out of cluster affinity

    ##only have to update clusters that change

    ##candidates are sites from other clusters that exceed threshold

    ##data.frame, element with highest affinity gain by changing clusters
    ##must be recalculated every loop
    ind <- clust_ind$elem + (clust_ind$clust - 1) * nrow(clust_ind)
    clust_current <- clust_aff[ind]
    gain <- clust_aff - clust_current
    if(length(locked_sites ) > 0) {
      gain[locked_sites, ] <- 0
    }
    max_ind <- which.max(gain)
    max_ind_zero <- (max_ind -1)
    i <- 1 + (max_ind_zero %% nrow(clust_ind)) #modulo is zero based indexing, R is 1 based indexing
    from <- clust_ind$clust[i]
    to <- (max_ind - i)/nrow(clust_ind) + 1
    upd <- data.frame(i = i, from = from, to = to, gain = gain[max_ind])

    ##Apply updates
    if(upd$gain > 0){
      ## message("applying update [", iter,"]")
      ## message(paste(names(upd), collapse =" - "))
      ## message(paste(upd, collapse = " - "))


      to_clust <- c(cast_obj[[upd$to]], upd$i)
      clust_aff_to <- mean(sim_mat[to_clust, to_clust])
      from_clust <- cast_obj[[ upd$from ]][cast_obj[[ upd$from ]] != upd$i  ] 
      clust_aff_from <- mean(sim_mat[from_clust, from_clust])

      if(clust_aff_to >= aff_thres && clust_aff_from >= aff_thres){
        cast_obj[[upd$to]] <- to_clust
        cast_obj[[upd$from]] <- from_clust

        clust_ind[upd$i, 2] <- upd$to
        clust_aff[, upd$to ] <- rowMeans(sim_mat[, cast_obj[[upd$to]], drop = FALSE] )
        clust_aff[, upd$from ] <-rowMeans(sim_mat[, cast_obj[[upd$from]], drop = FALSE])

        locked_sites <- c()

      } else {
        if(length(locked_sites ) > 0) {
          locked_sites <- c(locked_sites, upd$i)
        } else {
          locked_sites <- upd$i
        }
        message("clust_aff_to [", clust_aff_to ,"] may have dropped below threshold [", aff_thres, "]")
        message("clust_aff_from [", clust_aff_from ,"] may have dropped below threshold [", aff_thres, "]")
      }

    } else {
      break
    }
    iter <- iter + 1
  }
  return(cast_obj)
}


#' Two functions that extend cast to reduce the
#' number of clusters
#'
#' cast_alg_cautious stabilizes the matrix after every cluster is added,
#' so we don't create a lot of small clusters on the edges.
#' I expect to use all the sites faster.
#'
#' cast_compact takes a stabilized cast object, and deletes
#' one cluster, assignes the sites to nearest clusters,
#' stabilizes again, and
#' checks whether the average cluster affinity
#' remains above the threshold. If the threshold is still valid,
#' keep the changed cast algorithm, otherwise revert and try another cluster.
cast_alg_cautious <- function(sim_mat, aff_thres, max_iter = 20){
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
    valid_seeds <- spares

    aff_local <- rep(0, n)

    ##Do the next steps until stability is reached
    repeat{
      is_change <- FALSE



      ##find initial cluster point
      if(all(!new_clust)){
        maxima <- which(spares & valid_seeds)
        new_elem <- maxima[sample.int(length(maxima), 1)] ##intention: take one of the max at random
        next_seed <- new_elem
        new_clust[new_elem] <- TRUE
        spares[new_elem] <- FALSE
        is_change <- TRUE
        debug_spares_expected <- debug_spares_expected - 1
        ## print(list(new_elem, sum(spares), debug_spares_expected, "first"))
        assertthat::assert_that(debug_spares_expected == sum(spares))
        ##update aff_local
        ## aff_local[new_clust | spares] <- aff_local[new_clust | spares] + aff_func(new_clust | spares, new_elem, sim_mat)
        aff_local <- aff_local + aff_func(new_clust | spares, new_elem, sim_mat)
      }

      ##addition stage
      ##affinity seems to be used to bias group size
      ##Strong affinity is needed to join a large group

      ##Short-circuit evaluation, won't try to find max of an empty set
      while(any(spares) && (max(aff_local[spares]/sum(new_clust)) >= aff_thres)){
        maxima <- which(aff_local == max(aff_local[spares]) & spares)
        new_elem <- maxima[sample.int(length(maxima), 1)] ##intention: take one of the max at random
        new_clust[new_elem] <- TRUE
        spares[new_elem] <- FALSE
        debug_spares_expected <- debug_spares_expected - 1
        ## print(list(new_elem, sum(spares), debug_spares_expected, "add"))
        assertthat::assert_that(debug_spares_expected == sum(spares))
        is_change <- TRUE
        ##update aff_local
        ## aff_local[new_clust | spares] <- aff_local[new_clust | spares] + aff_func(new_clust | spares, new_elem, sim_mat)
        aff_local <- aff_local + aff_func(new_clust | spares, new_elem, sim_mat)
      }

      ##Removal stage
      while(any(new_clust) && (min(aff_local[new_clust]/sum(new_clust)) < aff_thres)){
        minima <- which(aff_local == min(aff_local[new_clust]) & new_clust)
        new_elem <- minima[sample.int(length(minima), 1)] ##intention: take one of the max at random
        new_clust[new_elem] <- FALSE
        spares[new_elem] <- TRUE
        debug_spares_expected <- debug_spares_expected + 1
        ## print(list(new_elem, sum(spares), debug_spares_expected, "remove"))
        assertthat::assert_that(debug_spares_expected == sum(spares))
        is_change <- TRUE
        ##update aff_local
        aff_local <- aff_local - aff_func(new_clust | spares, new_elem, sim_mat)
      }

      if(all(!new_clust)){
        ##cluster is empty
        valid_seeds[next_seed] <- FALSE
        message("seeds left: [", sum(valid_seeds), "]")
        if(all(!valid_seeds)){
          ##no more valid seeds exist, all have been tried, create a leftovers cluster
          new_clust <- spares
          spares <- rep(FALSE, n)
          break
        }
      }

      if(!is_change){
        break
      }
    }



    message("cluster assigned of size [", sum(new_clust), "]")
    message("[", sum(spares), "] sites left to assign.")
    assertthat::assert_that(debug_spares_expected == (debug_spares_loop - sum(new_clust)))

    ##now, stabilize the clusters.
    clust[[clust_id]] <- which(new_clust)
    if(any(spares)){
      clust[[clust_id + 1]] <- which(spares)
    }
    clust <- cast_stabilize(cast_obj = clust, aff_thres = aff_thres, sim_mat = sim_mat, max_iter = max_iter)
    if(any(spares)){
      spares <- rep(FALSE, n)
      spares[clust[[clust_id +1]] ] <- TRUE
    }
    message("[", sum(spares), "] sites left to assign after stabilize step.")
    clust_id <- clust_id + 1
  }

  return(clust)
}
cast_compact <- function(cast_ob, sim_mat, aff_thres, max_iter = nrow(sim_mat)*2){
  ##sort by size
  cluster_size <- data.frame(clust = 1:length(cast_ob), size = sapply(cast_ob, length))
  cluster_size <- cluster_size[order(cluster_size$size),]
  i <- 1
  repeat{
    if(length(cast_ob) <= 1){
      ##Stop at one cluster
      break
    }
    cast_test <- cast_ob
    ##assign affinities if cluster i to other clusters
    clust_id <- cluster_size$clust[i]
    updates <- lapply(seq_along(cast_ob[[clust_id]]), function(u, cast_ob, sim_mat){
      clust_aff <- sapply(cast_ob, function(clust, site, sim_mat){
          ##get the affinity to each cluster
          if(length(clust) > 0){
            return(mean(sim_mat[site, clust]))
          } else {
            ##empty clusters have 0 affinity
            return(0)
          }
      }, site = cast_ob[[clust_id]][u], sim_mat = sim_mat)
      new_clust <- which(clust_aff == max(clust_aff[-clust_id]))[1]
      return(data.frame(site = cast_ob[[clust_id]][u], old = clust_id, new = new_clust))
      }, cast_ob = cast_ob, sim_mat = sim_mat)
    updates <- do.call("rbind", updates)
    ##Apply updates
      message("reassigned [", nrow(updates),"] samples from cluster [", clust_id, "]")
      for(upd in 1:nrow(updates) ){
        cast_test[[ updates[upd,"old"] ]] <- cast_test[[ updates[upd,"old"] ]][cast_test[[ updates[upd,"old"] ]] != updates[upd,"site"]  ] 
        cast_test[[updates[upd,"new"]]] <- c(cast_test[[updates[upd,"new"]]], updates[upd,"site"])
      }

    cast_test[[clust_id]] <- NULL

    ##test new cluster affinities
    new_aff <- do.call("c", aff_clust_inner(cast_obj = cast_test, sim_mat = sim_mat))
    message("min_new_aff: [", min(new_aff), "] and threshold of [", aff_thres, "]. attempt: [", i, "]")
    cast_test_stab <- cast_stabilize(cast_obj = cast_test, aff_thres = aff_thres, sim_mat = sim_mat, max_iter = max_iter )
    new_aff_stab <- do.call("c", aff_clust_inner(cast_obj = cast_test_stab, sim_mat = sim_mat))
    message("min_new_aff_stab: [", min(new_aff_stab), "] and threshold of [", aff_thres, "]")
    if(min(new_aff_stab) >= aff_thres && min(new_aff) >= aff_thres){
      cast_ob <- cast_test_stab
      cluster_size <- data.frame(clust = 1:length(cast_ob), size = sapply(cast_ob, length))
      cluster_size <- cluster_size[order(cluster_size$size),]
      i <- 0
    }

    i <- i + 1
    if(i > nrow(cluster_size)){
      break
    }
  }
  return(cast_ob)
}

##within cluster affinity
##returns a list
aff_clust_inner <- function(cast_obj, sim_mat){
  lapply(seq_along(cast_obj), function(clust, cast_obj, sim_mat){
    elem_cor <- sim_mat[cast_obj[[clust]], cast_obj[[clust]] ]
    mean_aff <- mean(elem_cor)
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
    if(length(elems_a) > 1 && length(elems_b) > 1){
      return(mean(rowSums(sim_mat[elems_a, elems_b])/length(elems_b)))
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
    if(length(cast_obj[[clust]]) > 1){
      return(data.frame(x_row = 1:nrow(new_sim_mat),
                        clust = clust,
                        aff = rowSums(new_sim_mat[, cast_obj[[clust]]])/length(cast_obj[[clust]])))
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
  reordered <- do.call("c", cast_h_obj[[1]]$clust[reordering])
  if(depth > 1){
    return(cast_h_reorder(cast_h_obj[-1], depth = depth - 1, reordering = reordered))
  } else {
    return(reordered)
  }
}

#' Calculate membership matix from cast object
#'
#' Helper function for hubert_gamma(),
#'
#' Produces a matrix where element (i,j)
#' is 1 if site i and site j are in the same
#' cluster, otherwise 0.
#'
#' is_long toggles between long form and wide form
#' matrix for return values.
#' For very large areas, a full nxn matrix may
#' not fit it memory. However, the matrix becomes
#' sparse, as k increases, so long form should always work.
membership_mat <- function(cast_ob, is_long = FALSE){
  n_sites <- do.call(sum, lapply(cast_ob, length))

  mat_long  <-do.call(rbind, lapply(cast_ob, function(clust){
    mat_long_clust <- expand.grid(clust, clust)
    return(mat_long_clust)
  })
  )
  names(mat_long) <- c("i", "j")

  if(is_long){
    return(data.frame(mat_long, clust = 1))
  } else {
    ##I need side effects to avoid creating k large matricies
    mat_wide <- matrix(0, nrow = n_sites, ncol = n_sites)
    for (clust in cast_ob) {
      mat_wide[clust, clust] <- 1
    }
    return(mat_wide)
  }
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
#' @param sim_mat must be a square matrix equivalent,
#' @param 
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

  h_gamma <- mean(sim_mat * member_mat) #not exact, Tseng and Kao 2003 calculate this over just the upper triangular matrix
  return(h_gamma)
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
#' @param sim_mat can be either a square matrix equivalent,
#' or sparse, with three columns. First column is i, then j, then similarity
#' @param 
hubert_gamma_out_of_scope <- function(sim_mat, member_mat, is_long = FALSE, norm_z = TRUE){
  if(dim(sim_mat)[1] != dim(sim_mat)[2]){
    if (dim(sim_mat)[2] == 3) {
      is_long_sim <- TRUE
      k_est_s <- max(sim_mat[,1:2])
    } else {
      stop(paste0("hubert_gamma: sim_mat must be either square or data.frame(i, j, sim). Dims: ", dim(sim_mat), ".\n"))
    }
  } else {
    is_long_sim <- FALSE
    k_est_s <- max(dim(sim_mat))
  }

  if(dim(member_mat)[1] != dim(member_mat)[2]){
    if (dim(member_mat)[2] == 2) {
      is_long_member <- TRUE
      k_est_m <- max(sim_mat[,1:2])
    } else {
      stop(paste0("hubert_gamma: member_mat must be either square or data.frame(i, j). Dims: ", dim(member_mat), ".\n"))
      }
    } else {
      is_long_member <- FALSE
      k_est_m <- max(dim(member_mat))
    }

  assertthat::assert_that(k_est_m == k_est_s)
  k_est <- k_est_m
  if(norm_z){
    if(is_long_sim){
      mean_s <- sum(sim_mat[,3]) / (k_est^2)

      non_zero <- nrow(sim_mat[,3])
      n_zero <- k_est^2 - non_zero
      var_s <- (sum((sim_mat[,3]-mean_s)^2) + mean_s^2 * n_zero) / k_est^2
      sd_s <- sqrt(var_s)
      sim_mat[,3] <- (sim_mat[,3] - mean_s) / sd_s
    } else {
      mean_s <- mean(sim_mat)
      sd_s <- sd(sim_mat)
      sim_mat <- (sim_mat - mean_s) / sd_s
    }
    if(is_long_member){
      non_zero <- nrow(member_mat[,3])
      n_zero <- k_est^2 - non_zero
      mean_m <- non_zero / (k_est^2)

      var_m <-  ((1 - mean_m)^2 * non_zero + mean_m^2 * n_zero) / k_est^2
      sd_m <- sqrt(var_s)
      clust_score <- list(same = ((1-mean_m) / sd_m), other = ((0-mean_m) / sd_m))
    } else {
      mean_m <- mean(member_mat)
      sd_m <- sd(member_mat)
      member_mat <- (member_mat - mean_m) / sd_m
    }
  } else {
    mean_m <- 0
    sd_m <- 1
    mean_s <- 0
    sd_s <- 1
  }
 
  if (is_long_sim & is_long_member){
    ## hard case, need to align indicies
    return()
  }

  if (!is_long_sim & !is_long_member){
    h_gamma <- mean(sim_mat * member_mat) #not exact, Tseng and Kao 2003 calculate this over just the upper triangular matrix
    return(h_gamma)
  }

  if (is_long_sim & !is_long_member){
    ##unroll the matrix
    ##R matrix column major, so i increments first
    sim_mat_ind <- sim_mat[,1] + (sim_mat[,2]-1)*k_est
    
    member_mat
    h_gamma <- 0
    return(h_gamma)
  }

}
