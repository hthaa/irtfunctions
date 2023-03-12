#' Root Mean Squared Error of Approximation (RMSEA).
#' @description Calculates the Root Mean Squared Error of Approximation (RMSEA) by way of the Chi-squared statistic, its Degrees of Freedom, and the sample size.
#' @param chi2 Chi-squared.
#' @param df The Chi-squared statistics Degrees of Freedom.
#' @param ss The sample size.
#' @return A single value: The RMSEA statistic.
RMSEA <- function(chi2, df, ss) {
  return(sqrt((chi2 - df)/(df * (ss - 1))))
}

CBA <- function(x) {
  #Calculates Cronbachs alpha for a data-frame or a matrix. Must be ordered so
  #that rows are respondents, and collumns are items.
  #
  # Args:
  #   x: a data frame of item responses.
  #
  # Returns:
  #   Cronbachs alpha for a data-frame or a matrix.

  n <- ncol(x)
  Vi <- apply(x, 2, var)
  Vtest <- var(rowSums(x))
  return((n / (n - 1)) * (1 - (sum(Vi) / Vtest)))
}

MAE <- function(x, y) {
  #Calculates the "Mean Absolute Error" of objects "x" and "y".
  #
  # Args:
  #   x: value or vector of values.
  #   y: value or vector of values.
  #
  # Returns:
  #   The mean absolute error between values or vectors "x" and "y".

  MAE.output <- mean(abs(x - y))
  return(MAE.output)
}

RMSE <- function(x, y) {
  #Calculates the "Root Mean Square Error" of objects "x" and "y". The RMSE is
  #conceptually related to MAE, but is more sensitive to outliers.
  #
  # Args:
  #   x: value or vector of values.
  #   y: value or vector of values.
  #
  # Returns:
  #   The root mean square error between values or vectors "x" and "y".

  RMSE.output <- sqrt(mean(((x - y)^2)))
  return(RMSE.output)
}

PCR <- function(b = 0, a = 1, c = 0, d = 1, theta = 0) {
  #Calculates the Probability of providing a response of "1" to a particular
  #item with particular properties under 1, 2 and 3 parameter logistic IRT
  #models (difficulty, discrimination and guessing), given the respondents
  #value of Theta.
  #
  # Args:
  #   difficulty: The difficulty parameter of the particular item.
  #   discrimination: The discrimination parameter of the particular item.
  #   guessing: The guessing parameter of the particular item.
  #   theta: The ability parameter of the particular respondent.
  #
  # Returns:
  #   A value specifying the probability of a particular individual giving a
  #   response of "1" to a particular item. As with all probabilities, the
  #   range of values is between 0 and 1.

  p <- c + (d - c) * (1 / (1 + (exp(-a * (theta - b)))))
  return(p)
}


ThetaVectorGenerator <- function(nresp        = 100,
                                 m.theta      = 0,
                                 sd.theta     = 1,
                                 min.theta    = -3,
                                 max.theta    = 3,
                                 distribution = "normal",
                                 force.range  = FALSE,
                                 f.r.method   = "censor") {
  #Takes input from a caller and generates a vector of Theta values, based on
  #the function-callers specifications.
  #
  # Args:
  #   nresp: The sample size, or length of the vector that is to be produced.
  #   m.theta: The mean of Theta in the population from which a sample of Theta
  #    values are drawn.
  #   sd.theta: The standard deviation of Theta in the population from which
  #    the sample is drawn.
  #   min.theta: The minimum value Theta is allowed to take.
  #   max.theta: The maximum value Theta is allowed to take.
  #    distribution: String. The form the population distribution. If "normal",
  #    sample from a population with a normally distributed Theta. If "uniform,
  #    sample from a population where theta is uniformly distributed. Default
  #    is "normal", option is "uniform".
  #   force.range: Logical. Force range of Theta if distribution = "normal".
  #    Default set to FALSE.
  #   f.r.method: The method by which specified range of Theta is enforced. If
  #    "censor", values above maximum are censored to the value of max.theta,
  #    and if below minimum, values are censored to the value of min.theta. If
  #    "resample", out-of-range entries are discarded and replaced with a new
  #    sampled value from the population, repeated until entry is no longer out
  #    of range. Default is "censor", option is "resample".
  #
  # Returns:
  #   A vector of values drawn from a population conforming to specifications.

  theta.vec <- NULL

  #Generate vector with normally distributed theta.
  if (distribution == "normal") {
    theta.vec <- rnorm(nresp, m.theta, sd.theta)
  }

  #Generate vector with normally distributed theta within restricted range.
  if (distribution == "normal" & force.range == TRUE) {
    if ((min.theta >= max.theta)) {
      stop("Sequence aborted: invalid range specified.")
    } else {
      #Check if theta vector values are within specified range, and modify
      # according to specifications if it is not.
      range.check <- TRUE
      while (range.check == TRUE) {
        if (f.r.method == "censor") {
          theta.vec <- ifelse(theta.vec > max.theta, max.theta, theta.vec)
          theta.vec <- ifelse(theta.vec < min.theta, min.theta, theta.vec)
          range.check <- FALSE
        }
        if (f.r.method == "resample") {
          theta.vec <- ifelse(theta.vec > max.theta | theta.vec < min.theta,
                              rnorm(1, m.theta, sd.theta), theta.vec)
          if (min(theta.vec) >= min.theta & max(theta.vec) <= max.theta) {
            range.check = FALSE
          }
        } else {
          stop("Sequence aborted: invalid f.r.method specified.")
        }
      }
    }
  }

  #Generate vector with uniformly distributed theta within restricted range.
  if (distribution == "uniform") {
    if (min.theta >= max.theta) {
      stop("Sequence aborted: invalid Theta range specified.")
    } else {
      theta.vec <- runif(ss, min.theta, max.theta)
    }
  }

  #Return assembled vector of Theta values.
  return(theta.vec)
}

ItemPropertiesGenerator <- function(nitem = 10,
                                    mindif = -3, maxdif = 3, rdifdt = FALSE,
                                    mindis = 1,  maxdis = 1, rdisdt = TRUE,
                                    mingss = 0,  maxgss = 0, rgssdt = TRUE,
                                    minmax = 1,  maxmax = 1, rmaxdt = TRUE) {
  #Generates matrix of item IRT-parameter (discrimination = b, difficulty = a,
  #guessing = c) values, based on specifications.
  #
  # Args:
  #   nitem: The number of items (or rows) for which parameters are to be
  #    specified.
  #   mindif: The minimum value the difficulty (b) parameter will take.
  #   maxdif: The maximum value the difficulty (b) parameter will take.
  #   rdifdt: Whether the difficulty parameter for each item will take on a
  #    random value between the minimum and maximum. Default set to FALSE. If
  #    TRUE, the difficulty parameter will vary randomly between the minimum
  #    and maximum. Default set to FALSE. If FALSE, the difficulty parameter
  #    will increase by a fixed amount from one item to the next.
  #   mindis: The minimum value the discrimination parameter (a) will take.
  #   maxdis: The maximum value the discrimination parameter (a) will take.
  #   rdisdt: Whether the discrimination parameter for each item will take on a
  #    random value between the minimum and maximum. Default set to TRUE. If
  #    FALSE, the discrimination parameter will increase by a fixed amount from
  #    one item to the next.
  #   mingss: The minimum value the guessing (c) parameter will take.
  #   maxgss: The maximum value the guessing (c) parameter will take.
  #   rgssdt: Whether the guessing parameter for each item will take on a
  #    random value between the minimum and maximum. Default set to TRUE. If
  #    TRUE, the guessing parameter will vary randomly between the minimum and
  #    maximum. If FALSE, the guessing parameter will increase by a fixed
  #    amount from one item to the next.
  #
  # Returns:
  #   A matrix of item parameters of the form:
  #         [, b] [, a] [, c] [, d]
  #   [1, ]
  #   [2, ]
  #    ...
  #   [J, ]

  item.props <- cbind(
    cbind(
      b = if (mindif - maxdif != 0) {
        if (rdifdt == FALSE) {
          seq(mindif, maxdif, (maxdif - mindif) / (nitem - 1))
        } else {
          runif (nitem, mindif, maxdif)
        }
      } else {
        rep(mindif, nitem)
      }
      ,
      a = if (mindis - maxdis != 0) {
        if (rdisdt == FALSE) {
          seq(mindis, maxdis, (maxdis - mindis) / (nitem - 1))
        } else {
          runif(nitem, mindis, maxdis)
        }
      } else {
        rep(mindis, nitem)
      }
    )
    ,
    c = if (mingss - maxgss != 0) {
      if (rgssdt == FALSE) {
        seq(mingss, maxgss, (maxgss - mingss) / (nitem - 1))
      } else {
        runif(nitem, mingss, maxgss)
      }
    } else {
      rep(mingss, nitem)
    }
  )
  item.props <- cbind(item.props, d = if (minmax - maxmax != 0) {
    if (rmaxdt == FALSE) {
      seq(minmax, maxmax, (maxmax - minmax) / (nitem - 1))
    } else {
      runif(nitem, minmax, maxmax)
    }
  } else {
    rep(minmax, nitem)
  }
  )
  return(item.props)
}

DichotomousDataGenerator <- function(iptable = itemprops, thetavector = theta) {
  #Generates a data-frame where collumns = number of items and rows = number of
  #respondents. Responses take the form of dichotomous values (0/1), where the
  #values of each individual observation is probabilistically generated based on
  #the respondents probability of responding "1" to the specific item with the
  #properties "b", "a", and "c".
  #
  # Args:
  #   iptable: table of item properties with columns [, b], [, a], [, c], and
  #    rows [1, ], [2, ], ..., [J, ].
  #   thetavector: vector of values representing each individuals' value of
  #    "theta".
  #
  # Dependencies: The PCR function.
  #
  # Returns:
  #   A data-frame of values 0 and 1, with number of columns = total number of
  #    items, and number of rows = rotal number of respondents.

  nitem <- nrow(iptable)
  nresp <- length(thetavector)
  dichotomous.data <- matrix(nrow = nresp, ncol = nitem)
  for (itemcount in 1:nitem) {
    pitemr <- PCR(iptable[itemcount, 1], iptable[itemcount, 2],
                  iptable[itemcount, 3], iptable[itemcount, 4], thetavector)
    itemresponse <- ifelse(runif(nresp, 0, 1) > pitemr, 0, 1)
    dichotomous.data[, itemcount] <- itemresponse
  }
  return(dichotomous.data)
}

ProbItemCorrect <- function(iptable = itemprops, thetavector = theta) {
  nitem <- nrow(iptable)
  nresp <- length(thetavector)
  dichotomous.data <- matrix(nrow = nresp, ncol = nitem)
  for (itemcount in 1:nitem) {
    pitemr <- PCR(iptable[itemcount, 1], iptable[itemcount, 2],
                  iptable[itemcount, 3], iptable[itemcount, 4], thetavector)
    dichotomous.data[, itemcount] <- pitemr
  }
  return(dichotomous.data)
}

ECA.lee <- function(a = l.cut, b = iptable, c = theta, d = 1) {
  #Calculates the marginal "Expected Classification Accuracy" under the Lee
  #approach using the cacIRT package.
  #
  # Args:
  #   a: The cutoff value of the theta-estimate used for classification.
  #   b: Table of item properties with collumns [, b], [, a], and [, c].
  #   c: Vector of theta values.
  #   d: Scaling parameter. Default set to 1 (no scaling). Setting the scaling
  #      parameter to 1.702, for example, would yield an approximated normal
  #      ogive scaled model.
  #
  # Dependencies: The cacIRT package.
  #
  # Returns:
  #   The marginal "Expected Classification Accuracy" under the Lee approach.

  if (!require("cacIRT")) install.packages("cacIRT")
  class.Lee(a, matrix(c(b[, 2], b[, 1], b[, 3]), ncol = 3), c,
            D = d)$Marginal[1, 1]
}

ECA.rud <- function(a = r.cut, b = iptable, c = theta, d = 1) {
  #Calculates the marginal "Expected Classification Accuracy" under the Rudner
  #approach using the cacIRT package.
  #
  # Args:
  #   a: The cutoff value of the theta-estimate used for classification.
  #   b: Table of item properties with collumns [, b], [, a], and [, c].
  #   c: Vector of theta values.
  #   d: Scaling parameter. Default set to 1 (no scaling). Setting the scaling
  #      parameter to 1.702, for example, would yield an approximated normal
  #      ogive scaled model.
  #
  # Dependencies: The cacIRT package.
  #
  # Returns:
  #   The marginal "Expected Classification Accuracy" under the Rudner approach.

  if (!require("cacIRT")) install.packages("cacIRT")
  class.Rud(a, matrix(c(b[, 2], b[, 1], b[, 3]), ncol = 3), c,
            D = d)$Marginal[1, 1]
}

OnePL <- function(data = irtdata,
                  ability = theta,
                  one.pl.f = FALSE,
                  one.pl.c = oneplc,
                  summary.output = TRUE,
                  item.props = itemprops,
                  discrete.ability = TRUE,
                  rmse.output = TRUE,
                  mae.output = TRUE,
                  theta.se.output = TRUE,
                  theta.pcr.corr = TRUE,
                  classification.output = TRUE,
                  nresp = length(theta)) {
  #Fits a One-Parameter Logistic Model to a dichotomus data-set using the LTM
  #package. Returns either a row of information regarding the overall fit of the
  #model estimates with the true input parameters, or a data-frame containing
  #information regarding the precision of the model in reproducing each of the
  #respondents true parameters.
  #
  # Args:
  #   data: The input data-set to which the model is to be fit.
  #   one.pl.f: Estimate a joint parameter value for each items' discrimination
  #    parameter (i.e., equality constraint on discrimination parameter
  #    estimate across items.
  #   one.pl.c: Specify a joint discrimination parameter value for all of the
  #    items. For a Rasch model, a parameter value of 1 would be specified.
  #   summary.output: Logical. If true, return vector of model summary fit
  #    statistics. If false, return table of individual respondent parameter
  #    estimates.
  #   item.props: Table of item properties of form [j, b], [j, a], and [j, c]
  #    Required for calculating accuracy of item difficulty estimates.
  #   discrete.ability: Vector of discrete theta parameters (e.g., -2, 0, 1).
  #    Required for calculating classification accuracy.
  #   rmse.output: Ask for the models root mean square error of its estimates.
  #   mae.output: Ask for the models mean absolute error of its estimates.
  #   theta.se.output: Ask for the standard error of respondent theta
  #    estimates. If summary.output = TRUE, returns the mean standard error of
  #    the theta estimates. If summary.output = FALSE, returns the vector of
  #    individual theta estimate standard errors.
  #   classification.output: Ask for output calculating the fit of the models
  #    accuracy of classifications (e.g. rate of false positives and negatives)
  #
  # Dependencies:
  #   The LTM package.
  #   The RMSE function.
  #   The MAE function.
  #
  # Returns:
  #   Either a row of information summarizing the performance of the 1PL model
  #   in reproducing known parameter values, or a data-frame containing each of
  #   the respondent parameter estimates.

  if (!require("ltm")) install.packages("ltm")
  #If equality constraint is not to be imposed, fit model with specified item
  # discrimination parameters.
  if (one.pl.f == FALSE) {
    one.param <- rasch(data, IRT.param = TRUE,
                       constraint = cbind(ncol(data) + 1, one.pl.c))
  } else {
  #If equality constraint is to be imposed, then do the following instead:
    one.param <- rasch(data, IRT.param = TRUE)
  }
  #Retrieve and store model-estimated factor scores.
  one.p.fac.score <- factor.scores(one.param, data)$score.dat$z1
  #Make a vector of discrete ability estimates.
  discrete.1pl.ab.est <- floor(one.p.fac.score)
  #Make two vectors registering false-negatives and false-positives.
  false.positive.1pl <- ifelse(discrete.ability >= discrete.1pl.ab.est, 0, 1)
  false.negative.1pl <- ifelse(discrete.ability <= discrete.1pl.ab.est, 0, 1)
  #If request summary output, assemble and return vector with summarizing
  # info regarding model precision.
  if (summary.output == TRUE) {
    summary.data <- NULL
     #If user asks for Root Mean Square Error (RMSE) output, assemble:
     if (rmse.output == TRUE) {
       rmse.data <- c("rmse.b.1pl" = RMSE(item.props[, 1], coef(one.param)[, 1]),
                      "rmse.a.1pl" = RMSE(item.props[, 2], coef(one.param)[, 2]),
                      "rmse.c.1pl" = RMSE(item.props[, 3], 0),
                      "rmse.theta.1pl" = RMSE(ability, one.p.fac.score)
                      )
       summary.data <- c(summary.data, rmse.data)
     }
    #If user asks for Mean Absolute Error (MAE) output, assemble:
    if (mae.output == TRUE) {
      mae.data <- c("mae.b.1pl" = MAE(item.props[, 1], coef(one.param)[, 1]),
                    "mae.a.1pl" = MAE(item.props[, 2], coef(one.param)[, 2]),
                    "mae.c.1pl" = MAE(item.props[, 3], 0),
                    "mae.theta.1pl" = MAE(ability, one.p.fac.score)
                    )
      summary.data <- c(summary.data, mae.data)
    }
    #If user ask for mean standard errors for theta estimates, assemble:
    if (theta.se.output == TRUE) {
      summary.data <- c(summary.data,
                        "mean.theta.SE.1pl" = mean(one.p.fac.score)
                        )
    }
    #If user asks for correlation between theta estimate and number of correct
    # responses, assemble:
    if (theta.pcr.corr == TRUE) {
      summary.data <- c(summary.data,
                        "theta.pcr.corr.1pl" = cor(rowMeans(data),
                                                   one.p.fac.score)
                        )
    }
    #If user asks for classification accuracy statistics, assemble:
    if (classification.output == TRUE) {
      ca.data <- c("FP.1pl.sum" = sum(false.positive.1pl),
                   "FP.1pl.mean" = mean(false.positive.1pl),
                   "FN.1pl.sum" = sum(false.negative.1pl),
                   "FN.1pl.mean" = mean(false.negative.1pl),
                   "miss.class.1pl.SUM" = sum(false.positive.1pl) +
                     sum(false.negative.1pl),
                   "miss.class.1pl.MEAN" = ((sum(false.positive.1pl) +
                                              sum(false.negative.1pl)) /
                                              nrow(data))
                   )
      summary.data <- c(summary.data, ca.data)
    }
    return(summary.data)
    #If user does not request summary info, return individual estimates.
  } else {
    individual.estimates <- c(one.p.fac.score, discrete.1pl.ab.est)
    individual.estimates <- c(individual.estimates, false.positive.1pl)
    individual.estimates <- c(individual.estimates, false.negative.1pl)
    return(individual.estimates)
  }
}

TwoPL <- function(data = irtdata,
                  ability = theta,
                  summary.output = TRUE,
                  item.props = itemprops,
                  discrete.ability = discreteability,
                  rmse.output = TRUE,
                  mae.output = TRUE,
                  theta.se.output = TRUE,
                  theta.pcr.corr = TRUE,
                  classification.output = TRUE,
                  nresp = length(theta)) {
  #Fits a Two-Parameter Logistic Model to a dichotomus data-set using the LTM
  #package. Returns either a row of information regarding the overall fit of the
  #model estimates with the true input parameters, or a data-frame containing
  #information regarding the precision of the model in reproducing each of the
  #respondents true parameters.
  #
  # Args:
  #   data: The input data-set to which the model is to be fit.
  #   summary.output: Logical. If true, return vector of model summary fit
  #    statistics. If false, return table of individual respondent parameter
  #    estimates.
  #   item.props: Table of item properties of the form [, b], [, a], and [, c]
  #    Required for calculating accuracy of item difficulty estimates.
  #   discrete.ability: Vector of discrete theta parameters (e.g., -2, 0, 1).
  #    Required for calculating classification accuracy.
  #   rmse.output: Ask for the models root mean square error of its estimates.
  #   mae.output: Ask for the models mean absolute error of its estimates.
  #   theta.se.output: Ask for the standard error of respondent theta
  #    estimates. If summary.output = TRUE, returns the mean standard error of
  #    the theta estimates. If summary.output = FALSE, returns the vector of
  #    individual theta estimate standard errors.
  #   classification.output: Ask for output calculating the fit of the models
  #    accuracy of classifications (e.g. rate of false positives and negatives)
  #
  # Dependencies:
  #   The LTM package.
  #   The RMSE function.
  #   The MAE function.
  #
  # Returns:
  #   Either a row of information summarizing the performance of the 2PL model
  #   in reproducing known parameter values, or a data-frame containing each of
  #   the respondent parameter estimates.

  if (!require("ltm")) install.packages("ltm")
  #Fit 2pl model to data.
  two.param <- ltm(data~z1, IRT.param = TRUE)
  #Retrieve and store model-estimated factor scores.
  two.p.fac.score <- factor.scores(two.param, data)$score.dat$z1
  #Make a vector of discrete ability estimates.
  discrete.2pl.ab.est <- floor(two.p.fac.score)
  #Make two vectors registering false-negatives and false-positives.
  false.positive.2pl <- ifelse(discrete.ability >= discrete.2pl.ab.est, 0, 1)
  false.negative.2pl <- ifelse(discrete.ability <= discrete.2pl.ab.est, 0, 1)
  #If request summary output, assemble and return vector with summarizing
  # info regarding model precision.
  if (summary.output == TRUE) {
    summary.data <- NULL
    #If user asks for Root Mean Square Error (RMSE) output, assemble:
    if (rmse.output == TRUE) {
      rmse.data <-
          c("rmse.b.2pl" = RMSE(item.props[, 1], coef(two.param)[, 1]),
            "rmse.a.2pl" = RMSE(item.props[, 2], coef(two.param)[, 2]),
            "rmse.c.2pl" = RMSE(item.props[, 3], 0),
            "rmse.theta.2pl" = RMSE(ability, two.p.fac.score)
          )
      summary.data <- c(summary.data, rmse.data)
    }
    #If user asks for Mean Absolute Error (MAE) output, assemble:
    if (mae.output == TRUE) {
      mae.data <-
          c("mae.b.2pl" = MAE(item.props[, 1], coef(two.param)[, 1]),
            "mae.a.2pl" = MAE(item.props[, 2], coef(two.param)[, 2]),
            "mae.c.2pl" = MAE(item.props[, 3], 0),
            "mae.theta.2pl" = RMSE(ability, two.p.fac.score)
          )
      summary.data <- c(summary.data, mae.data)
    }
    #If user ask for mean standard errors for theta estimates, assemble:
    if (theta.se.output == TRUE) {
      summary.data <- c(summary.data, "mean.theta.SE.2pl" =
                              mean(two.p.fac.score))
    }
    #If user asks for correlation between theta estimate and number of correct
    # responses, assemble:
    if (theta.pcr.corr == TRUE) {
      summary.data <- c(summary.data,
                            "theta.pcr.corr.2pl" = cor(rowMeans(data),
                                                       two.p.fac.score))
    }
    #If user asks for classification accuracy statistics, assemble:
    if (classification.output == TRUE) {
      ca.data <-
          c("FP.2pl.sum" = sum(false.positive.2pl),
            "FP.2pl.mean" = mean(false.positive.2pl),
            "FN.2pl.sum" = sum(false.negative.2pl),
            "FN.2pl.mean" = mean(false.negative.2pl),
            "miss.class.2pl.SUM" = sum(false.positive.2pl) +
              sum(false.negative.2pl),
            "miss.class.2pl.MEAN" = (sum(false.positive.2pl) +
              sum(false.negative.2pl)) / nrow(data)
          )
      summary.data <- c(summary.data, ca.data)
    }
    return(summary.data)
    #If user does not request summary info, return individual estimates.
  } else {
    individual.estimates <- c(two.p.fac.score, discrete.2pl.ab.est)
    individual.estimates <- c(individual.estimates, false.positive.2pl)
    individual.estimates <- c(individual.estimates, false.negative.2pl)
    return(individual.estimates)
  }
}

ThreePL <- function(data = irtdata,
                    ability = theta,
                    summary.output = TRUE,
                    item.props = itemprops,
                    discrete.ability = discreteability,
                    rmse.output = TRUE,
                    mae.output = TRUE,
                    theta.se.output = TRUE,
                    theta.pcr.corr = TRUE,
                    classification.output = TRUE,
                    nresp = length(theta)
                    ) {
  #Fits a Three-Parameter Logistic Model to a dichotomus data-set using the LTM
  #package. Returns either a row of information regarding the overall fit of the
  #model estimates with the true input parameters, or a data-frame containing
  #information regarding the precision of the model in reproducing each of the
  #respondents true parameters.
  #
  # Args:
  #   data: The input data-set to which the model is to be fit.
  #   summary.output: Logical. If true, return vector of model summary fit
  #    statistics. If false, return table of individual respondent parameter
  #    estimates.
  #   item.props: Table of item properties of the form [, b], [, a], and [, c]
  #    Required for calculating accuracy of item difficulty estimates.
  #   discrete.ability: Vector of discrete theta parameters (e.g., -2, 0, 1).
  #    Required for calculating classification accuracy.
  #   rmse.output: Ask for the models root mean square error of its estimates.
  #   mae.output: Ask for the models mean absolute error of its estimates.
  #   theta.se.output: Ask for the standard error of respondent theta
  #    estimates. If summary.output = TRUE, returns the mean standard error of
  #    the theta estimates. If summary.output = FALSE, returns the vector of
  #    individual theta estimate standard errors.
  #   classification.output: Ask for output calculating the fit of the models
  #    accuracy of classifications (e.g. rate of false positives and negatives)
  #
  # Dependencies:
  #   The LTM package.
  #   The RMSE function.
  #   The MAE function.
  #
  # Returns:
  #   Either a row of information summarizing the performance of the 2PL model
  #   in reproducing known parameter values, or a data-frame containing each of
  #   the respondent parameter estimates.

  if (!require("ltm")) install.packages("ltm")
  #Fit 3pl model to data.
  three.param <- tpm(data, type = "latent.trait", IRT.param = TRUE)
  #Retrieve and store model-estimated factor scores.
  three.p.fac.score <- factor.scores(three.param, data)$score.dat$z1
  #Make a vector of discrete ability estimates.
  discrete.3pl.ab.est <- floor(three.p.fac.score)
  #Make two vectors registering false-negatives and false-positives.
  false.positive.3pl <- ifelse(discrete.ability >= discrete.3pl.ab.est, 0, 1)
  false.negative.3pl <- ifelse(discrete.ability <= discrete.3pl.ab.est, 0, 1)
  #If request summary output, assemble and return vector with summarizing
  # info regarding model precision.
  if (summary.output == TRUE) {
    summary.data <- NULL
    #If user asks for Root Mean Square Error (RMSE) output, assemble:
    if (rmse.output == TRUE) {
      rmse.data <-
        c("rmse.b.3pl" = RMSE(item.props[, 1], coef(three.param)[, 1]),
          "rmse.a.3pl" = RMSE(item.props[, 2], coef(three.param)[, 2]),
          "rmse.c.3pl" = RMSE(item.props[, 3], coef(three.param)[, 3]),
          "rmse.theta.3pl" = RMSE(ability, three.p.fac.score)
          )
      summary.data <- c(summary.data, rmse.data)
    }
    #If user asks for Mean Absolute Error (MAE) output, assemble:
    if (mae.output == TRUE) {
      mae.data <-
        c("mae.b.3pl" = MAE(item.props[, 1], coef(three.param)[, 1]),
          "mae.a.3pl" = MAE(item.props[, 2], coef(three.param)[, 2]),
          "mae.c.3pl" = MAE(item.props[, 3], coef(three.param)[, 3]),
          "mae.theta.3pl" = MAE(ability, three.p.fac.score)
          )
      summary.data <- c(summary.data, mae.data)
    }
    #If user ask for mean standard errors for theta estimates, assemble:
    if (theta.se.output == TRUE) {
      summary.data <- c(summary.data, "mean.theta.SE.3pl" =
                              mean(three.p.fac.score))
    }
    #If user asks for correlation between theta estimate and number of correct
    # responses, assemble:
    if (theta.pcr.corr == TRUE) {
      summary.data <- c(summary.data,
                            "theta.pcr.corr.3pl" = cor(rowMeans(data),
                                                       three.p.fac.score))
    }
    #If user asks for classification accuracy statistics, assemble:
    if (classification.output == TRUE) {
      ca.data <-
        c("FP.3pl.sum" = sum(false.positive.3pl),
          "FP.3pl.mean" = mean(false.positive.3pl),
          "FN.3pl.sum" = sum(false.negative.3pl),
          "FN.3pl.mean" = mean(false.negative.3pl),
          "miss.class.3pl.SUM" = sum(false.positive.3pl) +
            sum(false.negative.3pl),
          "miss.class.3pl.MEAN" = (sum(false.positive.3pl) +
                                     sum(false.negative.3pl)) / nrow(data)
          )
      summary.data <- c(summary.data, ca.data)
      }
    return(summary.data)
    #If user does not request summary info, return individual estimates.
  } else {
    individual.estimates <- c(three.p.fac.score, discrete.3pl.ab.est)
    individual.estimates <- c(individual.estimates, false.positive.3pl)
    individual.estimates <- c(individual.estimates, false.negative.3pl)
    return(individual.estimates)
  }
}

IRTSim <- function(
  #Specify types of models to fit.
  Fit.OnePL = TRUE,    #Fit 1PL models to data.
  Fit.TwoPL = TRUE,    #Fit 2PL models to data.
  Fit.ThreePL = FALSE, #Fit 3PL models to data (3PL is error prone).

  #Specifications relating to respondents (rows).
  Ab.Feed = NULL,      #Feed ability vector. Renders resp. options below moot.
  From.I.Resp,         #The minimum number of respondents (rows).
  To.I.Resp,           #The maximum number of respondents (rows).
  I.Resp.Reps,         #The number of repetitions (determines increments).

    From.Theta.Mean = 0,  #The minimum value for population theta mean.
    To.Theta.Mean = 0,    #The maximum value for population theta mean.
    Theta.Mean.Reps = 1,  #The number of repetitions (determines increments).

    From.Theta.SD = 1,    #The minimum value for population theta SD.
    To.Theta.SD = 1,      #The maximum value for population theta SD.
    Theta.SD.Reps = 1,    #The number of repetitions (determines increments).

  #Specifications relating to items (collumns).
  It.Feed = NULL,     #Feed item-parameters. Renders item options below moot.
  From.J.Item,        #The minimum number of items (collumns).
  To.J.Item,          #The maximum number of items (collumns).
  J.Item.Reps,        #The number of repetitions (determines increments).

    From.Dis = 1,     #The minimum value for the discrimination parameter.
    To.Dis = 1,       #The maximum value for the discrimination parameter.
    Dis.Reps = 1,     #The number of repetitions (determines increments).

    From.Dif = -3,    #The minimum value for the difficulty parameter.
    To.Dif = 3,       #The maximum value for the difficulty parameter.
    Dif.Reps = 1,     #The number of repetitions (determines increments).

    From.Gss = 0,     #The minimum value for the guessing parameter.
    To.Gss = 0,       #The maximum value for the guessing parameter.
    Gss.Reps = 1,     #The number of repetitions (determines increments).

  #Specifications relating to output.
  ECA = FALSE,        #Include measures of estimated classification accuracy.
    Rudner = FALSE,   #Include Rudner's estimate of classification accuracy.
      r.cut = 0,      #Specify theta-cutoff for Rudner ECA.
    Lee = FALSE,      #Include Lee's estimate of classification accuracy.
      l.cut = .5      #Specify sumscore cutoff for the Lee ECA.

  ) {

  if (ECA == TRUE) {
    if (!require("cacIRT")) install.packages("cacIRT")
  }

  #Here, the number of ability alterations are performed during the simulations
  # is determined. Ab.Feed trumps all other specifications, meaning that if the
  # user specifies a pre-assembled vector of ability values, only a single
  # ability-condition is specified.
   ##Possible future addition; add several different pre-assembled ability
   ## vectors, and specify number of conditions based on how many different
   ## pre-assembled vectors are fed to the function.
  if (is.null(Ab.Feed) == FALSE) {
    resp.reps <- 1
  } else {
    resp.reps <- tonresp - fromnresp + 1
  }
  #Here, the number of item alterations are performed during the simulations is
  # determined. It.Feed trumps all other specifications, meaning that if the
  # user specifies a pre-assembled vector of item parameter values, only a
  # single item-condition is specified.
   ##Possible future addition; add several different pre-assembled ability
   ## tables, and specify number of conditions based on how many different
   ## pre-assembled tables are fed to the function.
  if (is.null(It.Feed) == FALSE) {
    item.reps <- 1
  } else {
    item.reps <- tonitem - fromnitem + 1
  }
  #Here, the number of discrimination-parameter value alterations are performed
  # during the simulations are determined. If min-max = 0, only the "From.Dis"
  # specification is operative.
  disreps  <- todis - fromdis + 1
  if (fromdis == 1 & fromdis - todis != 0) {
    disfeedmin <- mindis
    disfeedmax <- maxdis
  } else {
    disfeedmin <- fromdis
    disfeedmax <- fromdis
  }

  for (resp.rep in 1:resp.reps) {
    #Generate the theta-vector that will underlie the observations, and which the
    # estimated theta parameter values will be compared against.
    if (is.null(Ab.Feed) == TRUE) {
      theta <- ThetaVectorGenerator()
    } else {
      theta <- Ab.Feed
    }
    #TODO: Loop
    resp.feed <- resp.feed + inc.resp.by
  }

  for (item.rep in 1:item.reps) {
    if (is.null(It.Feed) == TRUE) {
      #Generate the table of item properties that will underlie the observations,
      # and which the estimated item parameter values will be compared against.
      item.props <- ItemPropertiesGenerator()
    }

    if (is.null(It.feed) == FALSE | is.matrix(It.feed) == FALSE) {
      stop("Item parameter feed must be of class 'matrix', with rows 1:J, and
           3 collumns (difficulty [b], discrimination [a], and guessing [c]).
           You can use the ItemPropertiesGenerator() function to generate
           a matrix of appropriate form.")
    }

    if (ncol(It.feed) != 3) {
      stop("Input matrix of item parameters must take the form
           [j, b] [j, a] [j, c]")
    }
    if (is.null(It.Feed) == FALSE) {
      item.props <- It.feed
    }
    #TODO: Loop
    item.feed <- item.feed + inc.item.by
  }

  #This is the main part of the assembler, constituting the inner-most loop. It
  # is the loop that generates data and fits models to that data, extracting
  # different summary-measures of precision pertaining to a whole host of
  # parameters.
  for (rep in 1:Reps) {
    if (is.null(Ab.Feed) == TRUE) {
      theta <- ThetaVectorGenerator()
    } else {
      theta <- Ab.Feed
    }
    if (is.null(It.Feed) == TRUE) {
      iptable <- ItemPropertiesGenerator()
    } else {
      iptable <- It.Feed
    }
    #First couple of collumns of output reserved for summary-information
    # regarding the conditions underlying the observations.
    output <-
        c(
          "nitem"= nitem, "nresp" = nresp, "M.theta" = mean(theta),
          "SD.theta" = sd(theta), "Mean.b" = mean(iptable[, 1]),
          "SD.b" = sd(iptable[, 1]), "Mean.a" = mean(iptable[, 2]),
          "SD.a" = sd(iptable[, 2]), "Mean.c" = mean(iptable[, 3]),
          "SD.c" = sd(iptable[, 3]), "c.b.alpha" = cronbach.alpha(data)$alpha,
          "theta.pcr.corr" = cor.test(rowMeans(data), ability)
          )
    if (ECA == TRUE) {
      output <- c(output, ECA.Lee())
    }
    if (Do.OnePL == TRUE) {
      output <- c(output, OnePL())
    }
    if (Do.TwoPL == TRUE) {
      output <- c(output, TwoPL())
    }
    if (Do.ThreePL == TRUE) {
      output <- c(output, ThreePL())
    }
  }
}
