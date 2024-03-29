% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gfbootstrap.R
\name{pred_names.combinedBootstrapGF}
\alias{pred_names.combinedBootstrapGF}
\title{Get predictor names from a combined gfbootstrap object}
\usage{
\method{pred_names}{combinedBootstrapGF}(obj)
}
\arguments{
\item{obj}{combinedBootstrapGF object}
}
\value{
vector of character strings, naming the predictors
}
\description{
Gradient Forest will drop variables that are
never used in any splits. During bootstrapping,
it is likely that at least one model drops a predictor.
}
\details{
pred_names() searches all models and finds the
full set of predictors, even if some models only
use a subset.
}
