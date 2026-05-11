analyze_data <- function(data, # list of simulated datasets, parameters, (and optionally formulas for gps)
                         models,
                         progress = FALSE) { # list of models to fit
    nmod <- length(models) # number of models
    N <- data$parameters$N
    ndat <- data$parameters$ndat
    mods <- vector("list", ndat)
    for (i in 1:nmod) { # loop over models
      mod <- models[[i]] # model to fit
      if(progress) cli::cli_progress_bar(paste0("Fitting ", names(models)[[i]]), total = ndat)
      if (names(models)[i] != "gps") { # if model is not gps (i.e. if it is a lavaan model)
        for(j in 1:ndat){
          mods[[j]] <- lavaan::lavaan( # fit model
            model = mod, # model to fit
            sample.cov = data$covmats[[j]], # dataset
            sample.mean = data$means[[j]],
            sample.nobs = N,
            meanstructure = T, # mean
            int.ov.free = T, # intercepts of observed variables not fixed to 0
          )
          if(progress) cli::cli_progress_update()
        }
      } else { # if model is gps
        mods <- map(data$datasets, function(x) { # loop over datasets
          p1() # update progress bar
          map(data$forms_gps, function(form) { # loop over formulas with gps
            lm(form, data = x) # fit linear model
          })
        })
      }
      data[[names(models)[i]]] <- mods # save models
    }
    data$models <- models # save models
    return(data) # return data
}
