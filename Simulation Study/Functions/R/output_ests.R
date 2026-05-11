output_ests <- function(data, # output from check_convergence
                            xyvar, # list with names that the x and y variable start with for each model
                            filtered = TRUE) {
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
            (parameterestimates(y)[x, "est"])) %>% unlist()) %>% # get estimates %>%
      lapply(as.data.frame) %>%
        bind_rows(.id = "timepoint") %>%
        set_names(c("timepoint", "estimate")) %>%
        mutate(timepoint = as.numeric(timepoint)+1)
      estimatesy <- # get estimates of y parameters
        map(yind, function(x) # loop over indices
          map(data[[name_mod]][whichmod], function(y) # loop over datasets
            (parameterestimates(y)[x, "est"])) %>% unlist()) %>% # get estimates
        lapply(as.data.frame) %>%
        bind_rows(.id = "timepoint") %>%
        set_names(c("timepoint", "estimate")) %>%
        mutate(timepoint = as.numeric(timepoint)+1)
    }
    results$output[[name_mod]][["estimatesx"]] <- estimatesx # save results
    results$output[[name_mod]][["estimatesy"]] <- estimatesy # save results
  }
  return(results) # return results
}
