remove_invalid <- function(data){
  ndat <- data$parameters$ndat

  #### FOR DPM_FREE
  valid_mods <- data$valid_mods_index$dpm_free
    nvalid <- length(valid_mods)
    # save ar values
    ar_x <- matrix(NA, nrow = nvalid, ncol = 4)
    ar_y <- matrix(NA, nrow = nvalid, ncol = 4)

    # save cross-regressive
    cr_x <- matrix(NA, nrow = nvalid, ncol = 4)
    cr_y <- matrix(NA, nrow = nvalid, ncol = 4)
    for(i in 1:nvalid){
      pars <- parameterEstimates(data[["dpm_free"]][[valid_mods[i]]])
      ar_x[i,] <- pars[c(13, 17, 21, 25), "est"]
      ar_y[i,] <- pars[c(16, 20, 24, 28), "est"]

      cr_x[i,] <- pars[c(14, 18, 22, 26), "est"]
      cr_y[i,] <- pars[c(15, 19, 23, 27), "est"]
    }
    ar_valid <- tibble(x_under = ar_x < 1,
                       x_over = ar_x > -1,
                       y_under = ar_y < 1,
                       y_over = ar_y > -1) %>%
      apply(1, any)
    cr_all <- cbind(cr_x, cr_y)
    cr_mean <- apply(cr_all, 2, mean)
    cr_sd <- apply(cr_all, 2, sd)
    cr_lower_bound <- cr_mean - 10*cr_sd
    cr_upper_bound <- cr_mean + 10*cr_sd
    cr_valid <- rep(NA, nvalid)
      for(i in 1:nvalid){
        cr_i <- cr_all[i,]
        cr_valid[i] <- all((cr_i < cr_upper_bound) & (cr_i > cr_lower_bound))
      }
    data[["valid_vals"]][["dpm_free"]][["ar"]] <- ar_valid
    data[["valid_vals"]][["dpm_free"]][["cr"]] <- cr_valid
    valid_vals <- (ar_valid & cr_valid)
    data$valid_mods_index$dpm_free <- valid_mods[valid_vals]

    ### FOR RICLPM_FREE
    valid_mods <- data$valid_mods_index$riclpm_free
    nvalid <- length(valid_mods)
    # save ar values
    ar_x <- matrix(NA, nrow = nvalid, ncol = 4)
    ar_y <- matrix(NA, nrow = nvalid, ncol = 4)

    # save cross-regressive
    cr_x <- matrix(NA, nrow = nvalid, ncol = 4)
    cr_y <- matrix(NA, nrow = nvalid, ncol = 4)
    for(i in 1:nvalid){
      pars <- parameterEstimates(data[["riclpm_free"]][[valid_mods[i]]])
      ar_x[i,] <- pars[c(53, 57, 61, 65), "est"]
      ar_y[i,] <- pars[c(56, 60, 64, 68), "est"]

      cr_x[i,] <- pars[c(54, 58, 62, 66), "est"]
      cr_y[i,] <- pars[c(55, 59, 63, 67), "est"]
    }
    ar_valid <- tibble(x_under = ar_x < 1,
                       x_over = ar_x > -1,
                       y_under = ar_y < 1,
                       y_over = ar_y > -1) %>%
      apply(1, any)
    cr_all <- cbind(cr_x, cr_y)
    cr_mean <- apply(cr_all, 2, mean)
    cr_sd <- apply(cr_all, 2, sd)
    cr_lower_bound <- cr_mean - 10*cr_sd
    cr_upper_bound <- cr_mean + 10*cr_sd
    cr_valid <- rep(NA, nvalid)
    for(i in 1:nvalid){
      cr_i <- cr_all[i,]
      cr_valid[i] <- all((cr_i < cr_upper_bound) & (cr_i > cr_lower_bound))
    }
    data[["valid_vals"]][["riclpm_free"]][["ar"]] <- ar_valid
    data[["valid_vals"]][["riclpm_free"]][["cr"]] <- cr_valid
    valid_vals <- which((ar_valid & cr_valid))
    data$valid_mods_index$riclpm_free <- valid_mods[valid_vals]
    return(data)
}

