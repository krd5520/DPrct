% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ITTtable.R
\name{ITTtable}
\alias{ITTtable}
\title{ITT table for multiple models}
\usage{
ITTtable(
  data,
  reg.models = NULL,
  response.vars = NULL,
  families = "gaussian",
  treat.vars,
  covariate.vars = NULL,
  control.var = NULL,
  bonferroni.npvals = NULL,
  mult.test.correct = c("treatments", "responses"),
  pivot.xtreat = TRUE,
  only.control.mean = TRUE,
  add.pval.stars = TRUE,
  model.names = NULL,
  stderr.func = NULL,
  include.df = FALSE,
  ...
)
}
\arguments{
\item{data}{is the data to use to fit the model}

\item{reg.models}{vector of strings or formulas for the regression models. If
\code{reg.models==NULL} (default), then \code{treat.vars} and \code{covariate.vars},
cannot be \code{NULL} and are used to make a linear regression model for each
\code{response.vars}. If \code{reg.models} is only one formula or string, then
the predictors of that formula are used with each \code{response.vars} to make
a vector of models.}

\item{response.vars}{a vector of names, indicators, or indices of the response
variables within \code{data}.}

\item{families}{a vector of the family for each \code{reg.models} \link{\code{glm}}. The
\code{length(families)} must be 1 or equal to \code{length(reg.models)} (or \code{response.vars}).
If the length is 1, then the same distribution family is used for all models. Default is
"gaussian".}

\item{treat.vars}{a vector of variable names, indices, or logicals for the
columns of the data that are treatment variables}

\item{covariate.vars}{a vector of names, indicators, or indices of the covariate
variables within \code{data} to create the regression model formulas if \code{reg.models==NULL}.}

\item{control.var}{the string for the name of the control variable in the data.
If \code{control.var==NULL} (default), then the control mean will not be included in the table}

\item{bonferroni.npvals}{is a numeric value for the number of p-values to use
a Bonferroni correction to adjust for multiple testing. The default is
the number of treatment variables and level coefficients in the model}

\item{mult.test.correct}{a vector with elements "treatments" or "responses". If
\code{bonferroni.npvals==NULL}, then \code{mult.test.correct} is used to determine
the bonferroni correction. If \code{mult.test.correct=="treatments"}, then
\code{bonferroni.npvals} is the number of treatment variables/levels. If
\code{mult.test.correct=="responses"}, then \code{bonferroni.npvals} is the
length of \code{reg.models} or \code{response.vars}. If \code{mult.test.correct}
is the vector of "treatments" and "responses", then \code{bonferroni.npvals} is
the number of treatment variables/levels plus the number of treatment variables/levels.}

\item{pivot.xtreat}{if \code{TRUE}, then the resulting data.frame is pivoted wider
(see \link{\code{pivot_wider}}) to have "Estimate", "StdErr", "AdjPvalue" columns for
each treatment variable/level (control.var if applicable). If \code{control.var!=NULL}
and \code{only.control.mean==TRUE}, then the control mean (not the StdErr or AdjPvalue)
is also a column. If \code{pivot.xtreat==FALSE},then data frame is left in long form
with a "Treatment" column.}

\item{only.control.mean}{is indicator if the wider pivot data.frame only includes
the StdErr and AdjPvalue columns.}

\item{add.pval.stars}{if \code{TRUE}, then standard errors will be formatted with
3 decimal places in square brackets followed by one "*" if unadjusted p-value is <0.1,
two if <0.05, and three if <0.01. For example, \code{"[0.110]**"}. If
\code{add.pval.stars==FALSE}, then standard error will remain as numeric elements.}

\item{model.names}{(optional) is a vector of names for each regression model.}

\item{stderr.func}{(optional) is a function that takes the \link{glm} model output
and returns the standard error of the estimated treatment coefficients. If NULL
then the standard errors from the summary table will be used.}

\item{...}{additional inputs will be passed to \code{\link{glm}}}
}
\value{
a data.frame with the coefficient estimates, standard errors, bonferroni
adjusted p-values for each treatment variable/level and each regression
model (or response variables).
}
\description{
ITT table for multiple models
}
