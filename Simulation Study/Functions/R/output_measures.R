output_measures <- function(data, # output from check_convergence
                            xyvar, # list with names that the x and y variable start with for each model
                            filtered = TRUE,
                            bias_measure = "mean",
                            rel_bias = TRUE) {
  # save results (and add to later)
  results <- list(
    parameters = data$parameters, # parameters
    models = data$models, # models
    convergence = data$convergence, # convergence
    valid_vals = data$valid_vals
  )
  phi <- data$parameters$phi # matrix of true autoregressive and cross-lagged effects
  t <- data$parameters$timepoints # number of timepoints
  gammax <- phi[1, 2] # true CL on x
  gammay <- phi[2, 1] # true CL on y
  nmod <- length(data$models) # number of models
  modnames <- names(data$models) # names of models
  nvalid <- rep(NA, nmod)
  names(nvalid) <- modnames
  for(i in 1:nmod){
    if(!filtered){
      whichmod = 1:length(data$datasets)
    } else {
      whichmod = data[["valid_mods_index"]][[i]]
    }
    ndat = length(whichmod)
    nvalid[i] = ndat
  }
  results$nvalid = nvalid
  for (i in 1:nmod) { # loop over models
    name_mod <- modnames[i] # name of model
    ndat = nvalid[i]
    if(!filtered){
      whichmod = 1:length(data$datasets)
    } else {
      whichmod = data[["valid_mods_index"]][[i]]
    }
    if (name_mod != "gps") { # if not gps model (i.e. lavaan models)
      xvar <- xyvar[[name_mod]]$xvar # x variable
      yvar <- xyvar[[name_mod]]$yvar # y variable
      ## get parameter estimates (to obtain indices of parameters)
      parest.index <- data[[name_mod]][[1]] %>%
        parameterestimates()
      if (name_mod == "clpm_lag2") { # if CLPM with lag 2
        xind <- c(2, 6, 14, 22) # indices of x parameters
        yind <- c(3, 9, 17, 25) # indices of y parameters
      } else if (name_mod == "clpm_lag2_covar"){ # if not CLPM with lag 2
        xind <- c(82, 166, 254, 342)
        yind <- c(123, 209, 297, 385)
      } else {
        xind <- # indices of effects on x
          which( #which parameters follow conditions
            grepl(paste0("^", xvar), parest.index$lhs) & #effect to x
              grepl(paste0("^", yvar), parest.index$rhs) & #effect from y
              parest.index$op == "~" #is a regression
          )
        yind <- # indices of effects on y
          which( #which parameters follow conditions
            grepl(paste0("^", yvar), parest.index$lhs) & #effect to y
              grepl(paste0("^", xvar), parest.index$rhs) & #effect from x
              parest.index$op == "~" #is a regression
          )
      }
      estimatesx <- # get estimates of x parameters
        map(xind, function(x) # loop over indices
          map(data[[name_mod]][whichmod], function(y) # loop over datasets
            (parameterestimates(y)[x, "est"])) %>% unlist()) # get estimates
      stderr_x <- # get estimates of x parameters
        map(xind, function(x) # loop over indices
          map(data[[name_mod]][whichmod], function(y) # loop over datasets
            (parameterestimates(y)[x, "se"])) %>% unlist()) # get estimates
      estimatesy <- # get estimates of y parameters
        map(yind, function(x) # loop over indices
          map(data[[name_mod]][whichmod], function(y) # loop over datasets
            (parameterestimates(y)[x, "est"])) %>% unlist()) # get estimates
      stderr_y <- # get estimates of y parameters
        map(yind, function(x) # loop over indices
          map(data[[name_mod]][whichmod], function(y) # loop over datasets
            (parameterestimates(y)[x, "se"])) %>% unlist()) # get estimates
      ci.varx <- # confidence interval of x estimate
        map(xind, function(j) # loop over indices
          tibble( # create tibble
            ci.low = (map(data[[name_mod]][whichmod], function(i) # lower CI
              parameterestimates(i)[j, "ci.lower"]) %>% unlist()), # get lower CI
            ci.up = (map(data[[name_mod]][whichmod], function(i) # upper CI
              parameterestimates(i)[j, "ci.upper"]) %>% unlist()), # get upper CI
            ci.include = case_when(ci.low > gammax | # does CI include true value?
                                     ci.up < gammax ~ FALSE,
                                   .default = TRUE)
          ))
      ci.vary <- # confidence interval of y estimate
        map(yind, function(j) # loop over indices
          tibble( # create tibble
            ci.low = (map(data[[name_mod]][whichmod], function(i) { # lower CI
              parameterestimates(i)[j, "ci.lower"]}) %>% unlist()), # get lower CI
            ci.up = (map(data[[name_mod]][whichmod], function(i) { # upper CI
              parameterestimates(i)[j, "ci.upper"]}) %>% unlist()), # get upper CI
            ci.include = case_when(ci.low > gammay | # does CI include true value?
                                     ci.up < gammay ~ FALSE,
                                   .default = TRUE)
          ))
      power_x <- ## power of model
        map(xind, function(j){
          tibble(
            pval = (map(data[[name_mod]][whichmod], function(i) { # lower CI
              parameterestimates(i)[j, "pvalue"]}) %>% unlist()), # get lower CI
            power_05 = case_when(pval < 0.05 ~ TRUE,
                                 .default = FALSE),
            power_10 = case_when(pval < 0.1 ~ TRUE,
                                 .default = FALSE)
          )
        })
      power_y <- ## power of model
        map(yind, function(j){
          tibble(
            pval = (map(data[[name_mod]][whichmod], function(i) { # lower CI
              parameterestimates(i)[j, "pvalue"]}) %>% unlist()), # get lower CI
            power_05 = case_when(pval < 0.05 ~ TRUE,
                                 .default = FALSE),
            power_10 = case_when(pval < 0.1 ~ TRUE,
                                 .default = FALSE)
          )
        })
    } else { # if gps model
      ndat <- data$parameters$ndat # number of datasets
      xind <- 1:(t - 1) # indices of x parameters
      yind <- t:(2 * t - 2) # indices of y parameters
      estimatesx <- map(xind, function(a) { # loop over x indices
        map(data$gps, function(b) { # loop over models for datasets
          coefficients(b[[a]])[2] # get estimate
        }) %>% unlist() # unlist
      })
      stderr_x <- map(xind, function(a) { # loop over x indices
        map(data$gps, function(b) { # loop over models for datasets
          sqrt(vcov(b[[a]])[2, 2]) # get SE
        }) %>% unlist() # unlist
      })
      estimatesy <- map(yind, function(a) { # loop over y indices
        map(data$gps, function(b) { # loop over models for datasets
          coefficients(b[[a]])[2] # get estimate
        }) %>% unlist() # unlist
      })
      stderr_y <- map(yind, function(a) { # loop over y indices
        map(data$gps, function(b) { # loop over models for datasets
          sqrt(vcov(b[[a]])[2, 2]) # get SE
        }) %>% unlist() # unlist
      })
      ci.varx <- map(xind, function(a) { # loop over x indices
        tibble( # create tibble
          ci.low = (map(data$gps, function(i){ # lower CI
            confint(i[[a]])[2, 1]
          }) %>% unlist()),
          ci.up = (map(data$gps, function(i){ # upper CI
            confint(i[[a]])[2, 2]
          }) %>% unlist()),
          ci.include = case_when(ci.low > gammax | # does CI include true value?
                                   ci.up < gammax ~ FALSE,
                                 .default = TRUE)
        )
      })
      ci.vary <- map(yind, function(a) { # loop over x indices
        tibble( # create tibble
          ci.low = (map(data$gps, function(i){ # lower CI
            confint(i[[a]])[2, 1]
          }) %>% unlist()),
          ci.up = (map(data$gps, function(i){ # upper CI
            confint(i[[a]])[2, 2]
          }) %>% unlist()),
          ci.include = case_when(ci.low > gammay | # does CI include true value?
                                   ci.up < gammay ~ FALSE,
                                 .default = TRUE)
        )
      })
    }
    if(bias_measure == "mean"){
      estx <- map(estimatesx, mean) %>% unlist() # mean of x estimates
      esty <- map(estimatesy, mean) %>% unlist() # mean of y estimates
    } else if(bias_measure == "median"){
      estx <- map(estimatesx, median) %>% unlist() # median of x estimates
      esty <- map(estimatesy, median) %>% unlist() # median of y estimates
    }
    se_x <- map(stderr_x, mean) %>% unlist() # mean of x SE
    se_y <- map(stderr_y, mean) %>% unlist() # mean of y SE
    empsex <- map(estimatesx, sd) %>% unlist() # sd of x estimates (empirical SE)
    empsey <- map(estimatesy, sd) %>% unlist() # sd of y estimates (empirical SE)
    biasx <- map(estimatesx, ~ mean(. - gammax)) %>% unlist() # bias of x estimates
    biasy <- map(estimatesy, ~ mean(. - gammay)) %>% unlist() # bias of y estimates
    mcse.x <- empsex / sqrt(ndat) # Monte Carlo SE of x bias
    mcse.y <- empsey / sqrt(ndat) # Monte Carlo SE of y bias
    msex <- map(estimatesx, ~ mean((. - gammax) ^ 2)) %>% unlist() # MSE of x estimates
    msey <- map(estimatesy, ~ mean((. - gammay) ^ 2)) %>% unlist() # MSE of y estimates
    mcsemsex <- map(estimatesx, ~ sd((.-gammax)^2)/sqrt(ndat-1)) %>% unlist() # Monte Carlo SE of x MSE
    mcsemsey <- map(estimatesy, ~ sd((.-gammay)^2)/sqrt(ndat-1)) %>% unlist() # Monte Carlo SE of y MSE
    rmsex <- sqrt(msex) # RMSE
    rmsey <- sqrt(msey)
    mcse_rmse_x <- mcsemsex*(1/2)*1/(rmsex)
    mcse_rmse_y <- mcsemsey*(1/2)*1/(rmsey)
    cr.x <- unlist(map(ci.varx, ~ mean(.$ci.include))) # coverage rate of x estimates
    cr.y <- unlist(map(ci.vary, ~ mean(.$ci.include))) # coverage rate of y estimates
    mcse_cr_x <- sqrt(cr.x*(1-cr.x)/ndat)
    mcse_cr_y <- sqrt(cr.y*(1-cr.y)/ndat)
    power_x05 <- unlist(map(power_x, ~mean(.$power_05)))
    power_x10 <- unlist(map(power_x, ~mean(.$power_10)))
    power_y05 <- unlist(map(power_y, ~mean(.$power_05)))
    power_y10 <- unlist(map(power_y, ~mean(.$power_10)))
    mcse_power_x05 <- sqrt(power_x05*(1-power_x05)/ndat)
    mcse_power_x10 <- sqrt(power_x10*(1-power_x10)/ndat)
    mcse_power_y05 <- sqrt(power_y05*(1-power_y05)/ndat)
    mcse_power_y10 <- sqrt(power_y10*(1-power_y10)/ndat)
    res_all <- tibble( # create tibble
        variable = c(rep("x", times = (t - 1)), rep("y", times = (t - 1))), # variable (i.e. x or y)
        time = rep(2:t, times = 2), # timepoint
        true = c(rep(gammax, times = (t - 1)), rep(gammay, times = (t -
                                                                      1))), # true value
        estimate = c(estx, esty), # mean estimate
        bias_estimate = c(biasx, biasy), # bias
        bias_mcse = c(mcse.x, mcse.y), # Monte Carlo SE of bias
        se_emp = c(empsex, empsey), # empirical SE
        se_est = c(se_x, se_y), # mean SE
        mse_estimate = c(msex, msey), # MSE
        mse_mcse = c(mcsemsex, mcsemsey), # Monte Carlo SE of MSE
        rmse_estimate = c(rmsex, rmsey), #RMSE
        rmse_mcse = c(mcse_rmse_x, mcse_rmse_y), #MCSE of RMSE
        coverage_estimate = c(cr.x, cr.y), # coverage rate
        coverage_mcse = c(mcse_cr_x, mcse_cr_y), # mcse of coverage rate
        power05_estimate = c(power_x05, power_y05),
        power10_estimate = c(power_x10, power_y10),
        power05_mcse = c(mcse_power_x05, mcse_power_y05),
        power10_mcse = c(mcse_power_x10, mcse_power_y10)
      )
    if(rel_bias){
      res_all <- res_all %>%
        dplyr::mutate(
          rel_bias = bias_estimate / true, # relative bias
          rel_bias_mcse = bias_mcse/true
        )
    }
    results$output[[name_mod]] <- res_all # save results
  }
  return(results) # return results
}
