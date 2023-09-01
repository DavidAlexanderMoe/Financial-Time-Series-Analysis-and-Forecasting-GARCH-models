################################################################################
##
## File:    TSA-Predict-Student-Functions.R
## 
## Purpose: Prediction functions.
##
## Created: 2018.03.03
##
## Version: 2019.11.08
## 
################################################################################

################################################################################
## Functions to manage predictions from Arima()
################################################################################

.model <- 
function(object)
{    
  ##############################################################################
  ##
  ## Arguments
  ##  object:  
  ##
  ## Value
  ##
  ##############################################################################

  ## FUNCTION:
     
  #### Extract
  coef <- names(object$coef)
  xreg <- colnames(object$xreg)
  #### Make
  coef <- coef[ !( coef %in% xreg ) ] 
  ind.drift <- "drift" %in% xreg
  xreg <- xreg[ xreg != "drift" ]
  #### Answer
  list( drift = ind.drift, xreg = xreg, coef = coef )
}
# ------------------------------------------------------------------------------


.predict.t1 <- 
function(nobs, J, n.ahead)
{    
  ##############################################################################
  ## Purpose
  ##  Compute the argument "t" to get ex-post forecasts (obtained from 
  ##   .predict(..., fixed.n.ahead = TRUE)) as function of:
  ##   - nobs = number of observations in the time series
  ##   - J = how many ex-post forecasts to compute  
  ##   - n.ahead = horizon = how many steps ahead to go 
  ##  The function is implemented because 'hand' calculation of "t" is a bit 
  ##   tricky for students.
  ## Arguments
  ##  nobs:    (numeric[1]) number of observations.
  ##  J:       (numeric[1]) Number of ex-post forecasts to compute. 
  ##  n.ahead: (numeric[1]) Number of steps ahead (horizon).
  ## Value
  ##  The argument x extended some steps ahead.
  ##############################################################################

  ## FUNCTION:
     
  #### Answer
  nobs - J - n.ahead + 1
}
# ------------------------------------------------------------------------------


.predict.naive <- 
function(fit, J, n.ahead, 
  g = c("id", "sqrt", "log", "log10"))
{
  ##############################################################################
  ## Purpose
  ##  Extract ex-post 'naive' forecasts from the time series.
  ##  The function is implemented because 'hand' extraction of 'naive' forecasts
  ##   form data is a bit tricky for students.
  ## Arguments
  ##  fit:     (Arima) A model estimated from Arima().
  ##  J:       (numeric[1]) Number of ex-post forecasts to compute. 
  ##  n.ahead: (numeric[1]) Number of steps ahead (horizon).
  ##  g: (character[1]) type of transformation. Only "id", "sqrt", "log", 
  ##   "log10" are currently implemented.
  ## Value
  ##  (numeric[J]) 'naive' forecasts.
  ##############################################################################

  ## FUNCTION:
     
  #### Set naive depending on the estimated model
  d <- fit$arma[6]
  D <- fit$arma[7]
  S <- fit$arma[5]
  naive <- if (D > 0) { S } else if (d > 0) { 1 } else { 0 }
  
  #### Compute based on naive
  y <- fit$x
  if (naive == 0)
  {
    ans <- rep.int(mean(y), J)
  }
  else (naive > 0)
  {
    naive <- naive * (n.ahead %/% naive + (n.ahead %% naive > 0))
    ans <- y[(NROW(y) - J - naive + 1) : (NROW(y) - naive) ]
  }
  #### Transform
  gDef <- c("id", "sqrt", "log", "log10")
  g <- g[1]
  ####
  if ( !( g %in% gDef ) )
  {
    x1 <- paste( paste0("'", gDef, "'"), collapse = ", ")
    stop("Argument 'g' have to be one among ", x1)
  }
  else if ( g %in% c("sqrt", "log", "log10") )
  {
    #### Set the auxiliary function
    f.inv <- if ( g == "sqrt" ) { .sqrt.inv }
    else if ( g == "log" ) { .log.inv }
    else if ( g == "log10" ) { .log10.inv }
    
    #### To the original scale
    ans <- f.inv( ans )
  }    
  
  #### Answer
  ans
}
# ------------------------------------------------------------------------------


.predict.changing.h <- 
function(object, n.ahead, t, y, xreg)
{  
  ##############################################################################
  ## Arguments:
  ##  object:        An object created by Arima().
  ##  n.ahead:       Number of steps ahead (horizon).
  ##  t:             Forecast origin. Forecasts start from Time from t + 1. 
  ##                  It must be <= length(y).
  ##  y:             Data used to compute predictions. Use this argument if you 
  ##                  aim at getting forecasts beyond the estimation period (in 
  ##                  this case NROW(y) must be > NROW(object$x)).  
  ##  xreg:          External data used to compute predictions. Use this 
  ##                  argument if you aim at getting forecasts beyond the 
  ##                  estimation period (see argument 'y'). 
  ##############################################################################
  
  #### Settings
  nobs <- NROW(y)
  ind.xreg <- NROW(xreg) > 0
  ind.in  <- 1 : t
  ind.out <- (t + 1) : (t + n.ahead)
  
  #### Extract
  y1 <- y[ind.in]
  if ( ind.xreg )
  {
    xreg1 <- xreg[ ind.in, , drop = FALSE]
    xreg2 <- xreg[ind.out, , drop = FALSE]
  }
  else
  {
    xreg1 <- xreg2 <- NULL
  }

  #### Fit  
  ## Prepare to avoid issues with drift
  object1 <- object
  object1$x <- object1$x[ind.in]     
  object1$xreg <- object1$xreg[ind.in, , drop = FALSE] 
  ## Cases
  ind <- NROW( .model(object)$xreg ) > 0
  fit1 <- if ( ind )
  {
    Arima(y = y1, model = object1, xreg = xreg1)
  }
  else
  {
    Arima(y = y1, model = object1)
  }
    
  #### Predict
  x1 <- forecast::forecast(object = fit1, h = n.ahead, level = 95, xreg = xreg2)
  se <- (x1$upper - x1$mean) / qnorm(0.975)
  
  #### Store
  pred1 <- cbind(t = ind.out, mean = as.numeric(x1$mean), se = se)

  #### Answer
  as.data.frame(pred1)
}
# ------------------------------------------------------------------------------


.predict.fixed.h <- 
function(object, n.ahead, t, y, xreg)
{  
  #############################################################################
  ## Arguments:
  ##  object:        An object created by Arima()
  ##  n.ahead:       Number of steps ahead (horizon).
  ##  t:             Time from which to start forecasts (it must be <= length(y))
  ##  y:             Data used to compute predictions. Use this argument if you 
  ##                 aim at getting forecasts beyond the estimation period (in 
  ##                 this case NROW(y) must be > NROW(object$x)).  
  ##  xreg:          External data used to compute predictions. Use this argument 
  ##                 if you aim at getting forecasts beyond the estimation 
  ##                 period (see argument 'y')).
  #############################################################################
  
  #### Settings
  nobs <- NROW(y)
  ind.xreg  <- NROW(xreg) > 0
  mod <- .model(object)
  ind.drift <- mod$drift
  
  #### Full drift
  if (ind.drift) 
  {
    drift <- object$xreg[, "drift"]
    x.drift <- seq(from = drift[1], by = 1, length.out = nobs)
  }
  
  #### External regressors
  xreg.full <- if ( ind.xreg ) 
  {
    xreg[, mod$xreg, drop = FALSE]
  }  
  else
  {
    object$xreg[, mod$xreg, drop = FALSE]
  }  
  if ( ind.drift )
  {
    xreg.full <- if (NCOL(xreg.full) == 0)
    {
      cbind(drift = x.drift)
    }
    else
    {
      cbind(drift = x.drift, xreg.full)
    }
  }
  
  #### Initialize
  ##
  # n1 <- nobs - t + 1             ## Fixed 2019.10.13
  # ind <- t : nobs
  # n1 <- nobs - t                 ## Fixed 2020.06.05
  # ind <- t : (t + n1 - 1)
  n1 <- (nobs - n.ahead) - (t - 1)
  ind <- t : (nobs - n.ahead)
  if ( n1 <= 0 ) { return(NULL) }
  ## Pred
  pred1 <- matrix(NA, n1, 3)
  colnames(pred1) <- c("t", "mean", "se")

  #### Cycle
  pos <- 1
  for( t in ind )
  { 
    #### Extract
    ind.in  <- 1 : t
    ind.out <- (t + 1) : (t + n.ahead)
    y1 <- y[ind.in]
    if ( ind.xreg )
    {
      #### Old: fixed 2020.12.01
      # xreg1 <- xreg[ ind.in, , drop = FALSE]
      # xreg2 <- xreg[ind.out, , drop = FALSE]
      #### xreg1 and xreg2 must not include a drift column
      ind <- colnames(xreg.full) != "drift"
      xreg1 <- xreg.full[ ind.in, ind, drop = FALSE]
      xreg2 <- xreg.full[ind.out, ind, drop = FALSE]
    }
    else
    {
      xreg1 <- xreg2 <- NULL
    }

    #### Fit  
    ## Prepare to avoid issues with drift
    object1 <- object
    ## Cases
    object1$x <- object1$x[ind.in]
    ## object1$xreg <- object1$xreg[ind.in, , drop = FALSE]
    object1$xreg <- xreg.full[ind.in, , drop = FALSE]
    ind <- NROW( .model(object1)$xreg ) > 0
    fit1 <- if ( ind )
    {
      Arima(y = y1, model = object1, xreg = xreg1)
    }
    else
    {
      Arima(y = y1, model = object1)
    }
    
    #### Predict
    x1 <- forecast::forecast(object = fit1, h = n.ahead, level = 95, xreg = xreg2)
    se <- (x1$upper - x1$mean) / qnorm(0.975)
    #### Store
    pred1[pos, ] <- c(t + n.ahead, x1$mean[n.ahead], se[n.ahead])
    pos <- pos + 1
  }
  
  #### Answer
  as.data.frame(pred1)
}
# ------------------------------------------------------------------------------


.predict <- 
function(object, n.ahead, t, 
  y = NULL, xreg = NULL, fixed.n.ahead = TRUE)
{  
  ##############################################################################
  ## Arguments:
  ##  object:        An object created by Arima().
  ##  n.ahead:       Number of steps ahead (horizon).
  ##  t:             Time from which to start forecasts (it must be <= 
  ##                  length(y))
  ##  y:             Data used to compute predictions. Use this argument if we 
  ##                  aim at getting forecasts beyond the estimation period (in 
  ##                  this case NROW(y) must be > NROW(object$x)).  
  ##  xreg:          External data used to compute predictions. Use this 
  ##                  argument if we aim at getting forecasts beyond the 
  ##                  estimation period (see argument 'y').
  ##  fixed.n.ahead: Whether the number of steps ahead is retained fixed
  ##############################################################################

  ## FUNCTION:

  #### Settings
  fixed.n.ahead <- as.logical(fixed.n.ahead[1])

  #### Data
  if ( NROW(y) == 0 )
  {
    y <- object$x
  }
  if ( NROW(xreg) == 0 )
  {
    xreg <- .model(object)$xreg
    if ( NROW(xreg) > 0 )
    {
      xreg <- object$xreg[, xreg, drop = FALSE]
    }
  }

  #### Settings
  nobs <- NROW(y)

  #### Check
  if ( n.ahead <= 0 )
  {
    stop("Argument 'n.ahead' must be a positive integer")
  } 
  n.ahead <- round(n.ahead[1])
  if ( t > nobs )
  {
    stop("Argument 't' must be <= NROW(y)")
  } 

  #### n.ahead not fixed
  pred1 <- if (!fixed.n.ahead)
  {
    .predict.changing.h(object = object, n.ahead = n.ahead, t = t, 
      y = y, xreg = xreg)
  }
  else
  {
    .predict.fixed.h(object = object, n.ahead = n.ahead, t = t, 
      y = y, xreg = xreg)
  }
  
  #### Answer
  list( t = t, n.ahead = n.ahead, fixed.n.ahead = fixed.n.ahead, 
    pred = as.data.frame(pred1) )
}
# ------------------------------------------------------------------------------


.oeff.4.predict <- 
function(object, n.ahead, 
  delta = 0.7, type = c(".", ""))
{
  ##############################################################################
  ## Arguments:
  ##  object:        An object created by tso()
  ##  n.ahead:       Number of steps ahead.
  ##  delta:         Delta setting used in the tso() call.
  ##  type:          "" for predict()
  ##                 "." for .predict()
  ##############################################################################

  #### Settings
  type <- type[1]
  n.ahead <- n.ahead[1]
  
  #### Immediate return if no outliers are included
  if (  NROW(object$outliers) == 0 | !(type %in% c("", ".")) )
  {
    NULL
  }
    
  #### Select
  nobs <- NROW(object$y)
  freq <- frequency(object$y)
  pars <- coefs2poly(object$fit)
  mo   <- object$outliers
  
  #### Effects
  xreg <- outliers.effects(mo = mo, n = nobs + n.ahead, weights = FALSE, 
    delta = delta, pars = pars, freq = freq)
  if ( type == "" )
  {
    xreg <- xreg[(nobs + 1) : (nobs + n.ahead), , drop = FALSE]
  }
  
  #### Answer
  xreg
}
# ------------------------------------------------------------------------------


################################################################################
## Functions for prediction checking
################################################################################

.ErrorMeasures.old <- 
function(y, fit, 
  naive = "mean")
{
  ##############################################################################
  ## Arguments:
  ##  y:      (numeric[n]) time series of data.
  ##  fit:    (numeric[m]) forecasts.
  ##  naive:  forecasts to compute scaled measures.
  ## Remarks: 
  ##  1) Let m = NROW(fit) and n = NROW(y). In case m < n, only the last m 
  ##     elements in y are used to compute error measures
  ##############################################################################

  ## FUNCTION:
  
  #### Settings
  naive <- naive[1]
  
  #### Errors
  nf <- NROW(fit)
  y1 <- y[ (NROW(y) - nf + 1) : NROW(y) ]
  u  <- y1 - fit

  #### Error measures
  ME   <- mean( u )
  MAE  <- mean( abs(u) )
  RMSE <- sqrt( mean( u^2 ) )
  #### Percentage error measures
  if ( all(y1 > 0) )
  {
    ur  <- u / y1
    MPE    <- mean( ur )
    MAPE   <- mean( abs( ur ) )
    RMSPE  <- sqrt( mean( ur^2 ) )
    # urs <- y1 / ( 0.5 * (y1 + fit) )
    # SyMAPE <- mean( abs(urs) )
    # SyRMSPE <- sqrt( mean( urs^2 ) )
  }
  else
  {
    MPE    <- NULL
    MAPE   <- NULL
    RMSPE  <- NULL
    SMAPE  <- NULL
    RSMSPE <- NULL
  }

  #### Scaled error measures
  if (naive == "mean")
  {
    u1 <- y1 - mean(y1)
  }
  else if (naive == "past")
  {
    u1 <- diff(x = y, lag = 1)
    u1 <- u1[ (NROW(u1) - nf + 1) : NROW(u1) ]
  }
  else if ( is.numeric(naive) )
  {
    u1 <- diff(x = y, lag = naive)
    u1 <- u1[ (NROW(u1) - nf + 1) : NROW(u1) ]    
  }

  ####
  ScMAE  <- MAE / mean( abs(u1) )
  ScRMSE <- RMSE / sqrt( mean( u1^2 ) )
  
  ####
  c(ME = ME, MAE = MAE, RMSE = RMSE, 
    MPE = MPE, MAPE = MAPE, RMSPE = RMSPE, 
    ScMAE = ScMAE, ScRMSE = ScRMSE)
}
# ------------------------------------------------------------------------------


.ErrorMeasures <- 
function(y, fit, naive)
{ 
  #### Errors
  nf <- NROW(fit)
  y1 <- y[ (NROW(y) - nf + 1) : NROW(y) ]
  u  <- y1 - fit

  #### Error measures
  ME   <- mean( u )
  MAE  <- mean( abs(u) )
  RMSE <- sqrt( mean( u^2 ) )
  #### Percentage error measures
  if ( all(y1 > 0) )
  {
    ur  <- u / y1
    MPE    <- mean( ur )
    MAPE   <- mean( abs( ur ) )
    RMSPE  <- sqrt( mean( ur^2 ) )
  }
  else
  {
    MPE    <- NULL
    MAPE   <- NULL
    RMSPE  <- NULL
  }

  #### Scaled error measures
  u1 <- y1 - naive
  ScMAE  <- MAE / mean( abs(u1) )
  ScRMSE <- RMSE / sqrt( mean( u1^2 ) )
  
  ####
  c(ME = ME, MAE = MAE, RMSE = RMSE, 
    MPE = MPE, MAPE = MAPE, RMSPE = RMSPE, 
    ScMAE = ScMAE, ScRMSE = ScRMSE)
}
# ------------------------------------------------------------------------------


.MincerZarnowitz <- 
function(y, fit, 
  msg = "")
{
  #### Estimate
  lm1 <- lm( y ~ fit ) 
  #### vcov
  vcov    <- vcovHC(x = lm1, type = "const")
  vcovHC  <- vcovHC(x = lm1)
  vcovHAC <- vcovHAC(x = lm1)
  #### Coef
  cnames <- c("estimate", "s.e.", "tstat", "pvalue")
  coef    <- coeftest(x = lm1, vcov. = vcov);    colnames(coef) <- cnames
  coefHC  <- coeftest(x = lm1, vcov. = vcovHC);  colnames(coefHC) <- cnames
  coefHAC <- coeftest(x = lm1, vcov. = vcovHAC); colnames(coefHAC) <- cnames
  coef <- data.frame(coef[, 1, drop = FALSE], coef[, 2:4], HC = coefHC[, 2:4], 
    HAC = coefHAC[, 2:4], check.names = FALSE)
  #### Test
  test <- linearHypothesis(model = lm1, 
    hypothesis.matrix = c("(Intercept) = 0", "fit = 1"), vcov. = vcov)
  testHC <- linearHypothesis(model = lm1, 
    hypothesis.matrix = c("(Intercept) = 0", "fit = 1"), vcov. = vcovHC)
  testHAC <- linearHypothesis(model = lm1, 
    hypothesis.matrix = c("(Intercept) = 0", "fit = 1"), vcov. = vcovHAC)
  #### Print
  if ( msg != "" )
  {
    cat(msg)
    print( coef )
    cat(
      "       F stat:", test$"F"[2], ",", 
      "df: (", test$Df[2], ",", test$Res.Df[2], "),",
      "p-value:", test$"Pr(>F)"[2], "\n", 
      "(HC)  F stat:", testHC$"F"[2], ",", 
      "df: (", testHC$Df[2], ",", testHC$Res.Df[2], "),",
      "p-value:", testHC$"Pr(>F)"[2], "\n", 
      "(HAC) F stat:", testHAC$"F"[2], ",", 
      "df: (", testHAC$Df[2], ",", testHAC$Res.Df[2], "),",
      "p-value:", testHAC$"Pr(>F)"[2], "\n")
  }
  #### Answer
  list(x = test, xHC = testHC, xHAC = testHAC) 
}
# ------------------------------------------------------------------------------


.DieboldMariano <- 
function(e1, e2, h, power, msg = "")
{
  #### Compute
  x1 <- dm.test(e1 = e1, e2 = e2, h = h, power = power)
  #### Print
  if ( msg != "" )
  {
    cat(msg,
      "Horiz:", x1$parameter["Forecast horizon"], 
      ", Loss fct pow:", x1$parameter["Loss function power"], 
      ", Stat (L1-L2):", x1$statistic, "\n")
  }
  #### Answer
  x1
}
# ------------------------------------------------------------------------------


################################################################################
## Functions to manage variable transformations
################################################################################

.loglik <- 
function(fit, 
  g = c("id", "sqrt", "log", "log10"))
{
  ##############################################################################
  ## DESCRIPTION
  ##  Return the log-likelihood and the AIC for the original scale (y) of the 
  ##   variable when the model is estimated on a different scale (w).
  ##
  ## ARGUMENTS
  ##  fit: () a model fitted with stats::arima or forecast::Arima.
  ##  g: (character[1]) type of transformation.
  ##
  ## VALUE
  ##  (list) with components
  ##   $g (character[1]) 
  ##   $loglik (numeric[1]) 
  ##   $aic (numeric[1]) 
  ##############################################################################

  ## FUNCTION:

  #### Settings
  gDef <- c("id", "sqrt", "log", "log10")
  g <- g[1]
  
  #### Extract
  loglik <- fit$loglik
  aic <- fit$aic
  bic <- fit$bic
  
  ####
  if ( !( g %in% gDef ) )
  {
    x1 <- paste( paste0("'", gDef, "'"), collapse = ", ")
    stop("Argument 'g' have to be one among ", x1)
  }
  else if ( g %in% c("sqrt", "log", "log10") )
  {
    #### Time series (in the modeling scale)
    w <- fitted(fit) + residuals(fit)

    #### Set two auxiliary function
    if ( g == "sqrt" )
    {
      f.inv <- .sqrt.inv
      f.der <- .sqrt.der
    }
    else if ( g == "log" )
    {
      f.inv <- .log.inv
      f.der <- .log.der
    }
    else if ( g == "log10" )
    {
      f.inv <- .log10.inv
      f.der <- .log10.der
    }
    
    #### y (original scale) 
    y <- f.inv( w )
     
    #### Adjust
    x1 <- sum( log( f.der(y) ) )
    loglik1 <- loglik + x1
    aic1 <- aic - 2 * (loglik1 - loglik)
    bic1 <- bic - 2 * (loglik1 - loglik)
    loglik <- loglik1 
    aic <- aic1 
    bic <- bic1 
  }    
  
  #### Answer
  list(g = g, loglik = loglik, aic = aic, bic = bic)
}
# ------------------------------------------------------------------------------


.mean <- 
function(mean, sd,  
  g = c("id", "sqrt", "log", "log10"))
{
  ##############################################################################
  ## DESCRIPTION
  ##  Return the conditional mean of the original scale (y) of the variable.
  ##
  ## ARGUMENTS
  ##  mean: (numeric) mean of the variable in the modeling scale (w).
  ##  sd: (numeric[1]) standard deviation of the error component.
  ##  g: (character[1]) type of transformation.
  ##
  ## VALUE
  ##  (list) with components
  ##   $g (character[1]) 
  ##   $loglik (numeric[1]) 
  ##   $aic (numeric[1]) 
  ##############################################################################

  ## FUNCTION:

  #### Settings
  gDef <- c("id", "sqrt", "log", "log10")
  g <- g[1]
  
  ####
  if ( !( g %in% gDef ) )
  {
    x1 <- paste( paste0("'", gDef, "'"), collapse = ", ")
    stop("Argument 'g' have to be one among ", x1)
  }
  else if ( g %in% c("sqrt", "log", "log10") )
  {
    #### Set two auxiliary function
    mean <- if ( g == "sqrt" )
    {
      mean^2 + sd^2 
    }
    else if ( g == "log" )
    {
      exp( mean + 0.5 * sd^2 )
    }
    else if ( g == "log10" )
    {
      l <- log(10)
      mean <- mean * l
      sd <- sd * l
      exp( mean + 0.5 * sd^2 )      
    }    
  }
  
  #### Answer
  list(g = g, mean = mean)
}
# ------------------------------------------------------------------------------


.pred.bands <- 
function(pred, 
  alpha = 0.05,
  g = c("id", "sqrt", "log", "log10") )
{
  ##############################################################################
  ## DESCRIPTION
  ##  Return the mean, median and confidence bands of the original scale (y) 
  ##  of the variable.
  ##
  ## ARGUMENTS
  ##  pred: () and object produced by .predict().
  ##  alpha: (numeric[1]) sum of probabilities in the tails.
  ##  g: (character[1]) type of transformation. Only "id", "sqrt", "log", 
  ##   "log10" are currently implemented.
  ##
  ## VALUE
  ##  (list) with components
  ##   $g (character[1]) 
  ##   $mean (numeric)
  ##   $median (numeric)
  ##   $lower (numeric)
  ##   $upper (numeric)
  ##   $alpha (numeric)
  ##############################################################################

  ## FUNCTION:
  
  #### Settings
  g <- g[1]
  
  #### Mean/Median
  median <- pred$pred$mean
  mean   <- .mean(mean = pred$pred$mean, sd = pred$pred$se, g = g)$mean

  #### Bands
  z <- qnorm(1 - alpha/2)
  lower <- median - z * pred$pred$se
  upper <- median + z * pred$pred$se

  #### Set the inverse
  if ( g != "id" )
  {
    f.inv <- if ( g == "sqrt" ) { .sqrt.inv }
      else if ( g == "log" ) { .log.inv }
      else if ( g == "log10" ) { .log10.inv }
      else ( stop("Only g in 'id', 'sqrt', 'log', 'log10' are currently implemented") )
  
    #### Adjust median and bands
    median <- f.inv(median)
    lower  <- f.inv(lower)
    upper  <- f.inv(upper)
  }

  #### Answer
  list(g = g, t = pred$pred$t, mean = mean, median = median, 
    lower = lower, upper = upper, se = pred$pred$se, alpha = alpha) 
}
# ------------------------------------------------------------------------------


################################################################################
## Auxiliary functions
################################################################################

.sqrt.inv <- 
function(x)
{
  x^2  
}
# ------------------------------------------------------------------------------


.sqrt.der <- 
function(x)
{
  0.5 / sqrt(x)  
}
# ------------------------------------------------------------------------------


.log.inv <- 
function(x)
{
  exp(x)  
}
# ------------------------------------------------------------------------------


.log.der <- 
function(x)
{
  1 / x  
}
# ------------------------------------------------------------------------------


.log10.inv <- 
function(x)
{
  10^x  
}
# ------------------------------------------------------------------------------


.log10.der <- 
function(x)
{
  1 / (log(10) * x)  
}
# ------------------------------------------------------------------------------

