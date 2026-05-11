# function to check if gradient is low enough for convergence
check.grad <- function(grads) {
  return(!any(grads >= 0.001))
}
# function to make a matrix symmetric
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}
check_convergence <- function(data, # list with models to evaluate
                              latent = T, # logical vector indicating if model uses a latent variable
                              filter = T) { # logical indicating if we want to filter the output
  ndat <- data$parameters$ndat
  modnames <- names(data$models)[names(data$models)!="gps"] # get model names (except for gps model)
  nmod <- length(modnames) # number of models
  for (i in 1:nmod) { # iterate over models
    name_mod <- modnames[i] # get model name
    if (latent[i]) { # if model has a latent variable
      data$convergence[[name_mod]] <- # create a list with convergence information
        tibble( # create a tibble
          converged = map(data[[name_mod]], ~ lavInspect(object = .x, "converged")),
          grads.check = map(data[[name_mod]], ~ check.grad(lavInspect(
            .x, "optim.gradient"
          ))), # check if gradient is low enough
          cov.lat = map(data[[name_mod]], ~ lavInspect(.x, what = "cov.lv")), # get latent variable covariance matrix
          cov.obs = map(data[[name_mod]], ~ lavInspect(.x, what = "theta")) # get observed variable covariance matrix
        ) %>%
        mutate( # add new columns
          cov.lat = map(cov.lat, makeSymm), # make covariance matrix symmetric
          posdef.lat = map(cov.lat, is.positive.semi.definite), # check if covariance matrix of latent variables is positive semidefinite
          cov.obs = map(cov.obs, makeSymm), # make covariance matrix symmetric
          posdef.obs = map(cov.obs, is.positive.semi.definite), # check if covariance matrix of observed variables is positive semidefinite
          dataset_id = row_number() # add dataset id
        ) %>%
        unnest(cols = c("converged", "grads.check", "posdef.lat", "posdef.obs")) %>% # unnest columns
        dplyr::select(converged,
                      grads.check,
                      posdef.lat,
                      posdef.obs,
                      dataset_id) # select columns
    } else { # if model does not have a latent variable
      data$convergence[[name_mod]] <- # create a list with convergence information
        tibble( # create a tibble
          converged = map(data[[name_mod]], ~ lavInspect(object = .x, "converged")), # check if model converged
          grads.check = map(data[[name_mod]], ~ check.grad(lavInspect(
            .x, "optim.gradient"
          ))), # check if gradient is low enough
          cov.obs = map(data[[name_mod]], ~ lavInspect(.x, what = "theta")) # get observed variable covariance matrix
        ) %>%
        mutate( # add new columns
          cov.obs = map(cov.obs, makeSymm), # make covariance matrix symmetric
          posdef.obs = map(cov.obs, is.positive.semi.definite), # check if covariance matrix of observed variables is positive semidefinite
          dataset_id = row_number() # add dataset id
        ) %>%
        unnest(cols = c("converged", "grads.check", "posdef.obs")) %>% # unnest columns
        dplyr::select(converged, grads.check, posdef.obs, dataset_id) # select columns
    }
    if (filter) { # if we want to filter the output
      filtname <- paste0(name_mod, "_filt") # create a new name for the filtered data
      if (latent[i]) { # if model has a latent variable
        # get indices of converged models with positive definite covariance matrices
        converged.index <- data$convergence[[name_mod]] %>%
          dplyr::filter(converged & grads.check & posdef.lat & posdef.obs) %>% # check requirements
          pull(dataset_id) # get dataset ids
        data[["valid_mods_index"]][[name_mod]] <- converged.index
      } else{ # if model does not have a latent variable
        # get indices of converged models with positive definite covariance matrices
        converged.index <- data$convergence[[name_mod]] %>%
          dplyr::filter(converged & grads.check & posdef.obs) %>% # check requirements
          pull(dataset_id) # get dataset ids
        data[["valid_mods_index"]][[name_mod]] <- converged.index
      }
    }
  }
  return(data) # return data
}
