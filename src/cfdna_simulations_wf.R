mysandbox <- function(){
  patient <- "009"
  ca125 <- filter(cfdna$ca125, str_detect(patient_id, patient)) %>% 
    mutate(x = first(CA125),
           CA125norm = CA125 / x) %>% 
    mutate(CA125norm = ifelse(CA125norm > 1, 1, CA125norm)) %>% 
    mutate(population = 10^9 * CA125norm)
  
  N <- seq(min(ca125$days), max(ca125$days)-10, by = 4)
  
  pop <- interpolate_population(ca125, N)
  
  ca125 %>% 
    select(days, population)
  
  print(ca125, n = 100)
  ggplot(ca125, aes(x = days, y = CA125norm)) +
    geom_point()
  
  
  
  # Parameters
  P <- 3   # Number of populations
  G <- nrow(pop)  # Number of generations
  p_vector <- c(0.5, 0.3, 0.7)  # Initial allele frequencies for each population
  
  # Create a matrix of population sizes for each generation for each population
  # Here, for simplicity, we're assuming the population sizes change linearly over generations
  N_matrix <- matrix(nrow = G , ncol = P)
  for (population in 1:P) {
    N_matrix[, population] <- pop$population # Example: population size changes from 100 to 200
  }
}


interpolate_population <- function(data, new_days) {
  # Check if new_days is within the range of existing days
  if (any(new_days < min(data$days)) || any(new_days > max(data$days))) {
    stop("new_days should be within the range of existing days.")
  }
  
  # Perform interpolation
  interpolated_values <- approx(data$days, data$population, xout = new_days)
  
  # Create a dataframe with interpolated values
  interpolated_df <- data.frame(days = interpolated_values$x, population = as.integer(round(interpolated_values$y)))
  
  return(interpolated_df)
}

wright_fisher_dynamic <- function(P, G, p_vector, N_vector, days, clones = NULL) {
  
  if (is.null(clones)){
    clones <- LETTERS[1:P]
  }
  
  # Initialize a matrix to store allele frequencies for each population over generations
  allele_freq_matrix <- matrix(0, nrow = G, ncol = P)
  
  # Set initial frequencies for each population
  allele_freq_matrix[1, ] <- p_vector
  
  # Simulate each population separately
  #for (population in 1:P) {
  for (generation in 2:(G)) {
    # Population size for the current generation
    N <- N_vector[generation]
    
    # Number of copies of the allele in the current generation
    #num_alleles <- rbinom(1, as.integer(N), allele_freq_matrix[generation - 1, population])
    num_alleles <- rmultinom(1, as.integer(N_vector[generation]), allele_freq_matrix[generation - 1, ])
    
    # Update the allele frequency
    allele_freq_matrix[generation, ] <- num_alleles / N
  }
  #}
  
  allele_freq_df <- as.data.frame(allele_freq_matrix)
  names(allele_freq_df) <- clones
  
  allele_freq_df$days <- days
  
  allele_freq_df <- allele_freq_df %>% 
    pivot_longer(cols = -days, names_to = "clone", values_to = "frequency")
  
  return(allele_freq_df)
}

calculate_p_value <- function(observed_value, mean, sd, observed_sd = 0.05, alternative = "two.sided") {
  
  combined_mean <- observed_value - mean
  combined_sd <- sqrt(sd^2 + observed_sd^2)
  
  # Calculate the z-score
  z_score <- combined_mean / combined_sd
  
  # Calculate the Z-score
  #z_score <- (observed_value - mean) / sd
  
  # Calculate the p-value based on the alternative hypothesis
  if (alternative == "two.sided") {
    p_value <- 2 * pnorm(-abs(z_score))
  } else if (alternative == "less") {
    p_value <- pnorm(z_score)
  } else if (alternative == "greater") {
    p_value <- pnorm(-z_score)
  } else {
    stop("Invalid alternative hypothesis. Use 'two.sided', 'less', or 'greater'.")
  }
  
  return(p_value)
}

calculate_N <- function(df, Nstart, Nmrd) {
  # Identify the first and minimum values of CA125
  first_CA125 <- df$CA125[1]
  min_CA125 <- min(df$CA125)
  
  # Define N values for interpolation
  N_first <- Nstart
  N_min <- Nmrd
  
  # Calculate interpolation parameters
  slope <- (N_first - N_min) / (first_CA125 - min_CA125)
  intercept <- N_first - slope * first_CA125
  
  x <- c(first_CA125, min_CA125)
  y <- c(N_first, N_min)
  
  # Create a spline function
  spline_function <- splinefun(x, y, method = "natural")
  
  # Apply interpolation to calculate N
  df <- df %>%
    mutate(population_old = as.integer(slope * CA125 + intercept)) %>% 
    mutate(population = spline_function(CA125))
  
  return(df)
}

calculate_N_exponential <- function(df, Nstart, Nmrd) {
  # Identify the first and minimum values of CA125
  first_CA125 <- df$CA125[1]
  #first_CA125 <- max(df$CA125)
  min_CA125 <- min(df$CA125)
  
  # Define N values for interpolation
  N_first <- Nstart
  N_min <- Nmrd
  
  # Fit an exponential model
  # N = a * exp(b * CA125)
  # log(N) = log(a) + b * CA125
  log_N_first <- log(N_first)
  log_N_min <- log(N_min)
  
  # Calculate coefficients for the exponential model
  b <- (log_N_min - log_N_first) / (min_CA125 - first_CA125)
  a <- exp(log_N_first - b * first_CA125)
  
  # Apply the exponential model to calculate N
  df <- df %>%
    mutate(population = a * exp(b * CA125))
  
  return(df)
}

library(stats)

# Function to combine p-values using Fisher's method
combine_p_values_fisher <- function(p_values) {
  # Calculate the test statistic
  test_statistic <- -2 * sum(log(p_values))
  
  # Degrees of freedom
  df <- 2 * length(p_values)
  
  # Calculate the combined p-value
  combined_p_value <- pchisq(test_statistic, df, lower.tail = FALSE)
  
  return(combined_p_value)
}


sim_pop <- function(mypatient, Nstart = 1e9, Nmrd = 1e4, Nsims = 100){
  
  print(mypatient)
  
  svdfclone <- svdf_clone %>% 
    filter(str_detect(patient, mypatient)) %>% 
    filter(nclones == 1)
  
  clonedf_start <- svdfclone %>% 
    filter(timepoint == 1) %>% 
    mutate(cloneVAF = ifelse(cloneVAF == 0, purity/100, cloneVAF)) %>% 
    mutate(cloneVAF = cloneVAF / sum(cloneVAF))
  
  clonedf_end <- svdfclone %>% 
    filter(str_detect(patient, mypatient)) %>% 
    filter(nclones == 1) %>% 
    filter(purity > 0.0) %>% 
    arrange(days) %>%
    group_by(timepoint) %>% 
    mutate(maxf = max(cloneVAF)) %>% 
    ungroup() %>% 
    filter(maxf > 0) %>% 
    filter(timepoint == max(timepoint)) %>% 
    mutate(cloneVAF = cloneVAF / sum(cloneVAF)) %>% 
    select(clone_id, cloneVAF, days) %>% 
    rename(clone = clone_id)
  
  ca125_patient <- filter(ca125, str_detect(patient_id, mypatient))
  
  if (mypatient == "045"){
    ca125_patient <- bind_rows(data.frame(CA125 = max(ca125_patient$CA125), days = 1, patient_id = "SPECTRUM-OV-045"), ca125_patient)
  }
  
  ca125_patient <- calculate_N_exponential(ca125_patient, Nstart, Nmrd) %>% 
    mutate(population = ifelse(population > Nstart, Nstart, population)) %>% 
    mutate(population = ifelse(population < Nmrd, Nmrd, population))
  
  # fact <- (ca125$CA125[1] / min(ca125$CA125))
  # 
  # ca125 <- ca125 %>% 
  #   mutate(x = first(CA125),
  #          CA125norm = CA125^2 / x) %>% 
  #   mutate(CA125norm = CA125norm / first(CA125norm)) %>% 
  #   mutate(CA125norm = ifelse(CA125norm > 1, 1, CA125norm)) %>% 
  #   mutate(population = Nstart * CA125norm)
  
  N <- seq(min(ca125$days), max(svdfclone$days), by = 4)
  N[length(N) + 1] <- max(svdfclone$days)
  N <- N[N < max(ca125_patient$days)]
  N <- N[N > min(ca125_patient$days)]
  #ca125 <- filter(ca125, days <  max(N))
  
  pop <- interpolate_population(ca125_patient, N) %>% 
    filter(days <= max(clonedf_end$days))
  
  # Parameters
  G <- nrow(pop)  # Number of generations
  p_vector <- clonedf_start$cloneVAF  # Initial allele frequencies for each population
  P <- length(p_vector)   # Number of populations
  
  # Create a matrix of population sizes for each generation for each population
  # N_matrix <- matrix(nrow = G , ncol = P)
  # for (population in 1:P) {
  #   N_matrix[, population] <- pop$population # Example: population size changes from 100 to 200
  # }
  
  #create a vector of population sizes
  N_vector <- pop$population
  
  sim_df <- data.frame()
  for (x in 1:Nsims){
    allele_freq_df <- wright_fisher_dynamic(P, G, p_vector, N_vector, pop$days, clones = clonedf_start$clone_id)
    sim_df <- bind_rows(sim_df, allele_freq_df %>% mutate(Nsim = x))
  }
  
  sim_df_summary <- sim_df %>% 
    group_by(days, clone) %>% 
    summarise(f_mean = mean(frequency),
              f_sd = sd(frequency),
              f_low = quantile(frequency, 0.025),
              f_high = quantile(frequency, 0.975)) %>% 
    ungroup()
  
  cp <- config$clonecols
  names(cp) <- LETTERS[1:length(cp)]
  
  pvals <- inner_join(sim_df_summary %>% filter(days == max(days)), clonedf_end, by = "clone") %>% 
    mutate(pval = calculate_p_value(cloneVAF, f_mean, f_sd)) %>% 
    mutate(patient = mypatient)
  
  wfplot <- sim_df_summary %>% 
    ggplot(aes(x = days, fill = clone, col = clone)) +
    geom_line(aes(y = f_mean), linewidth = 0.2) +
    geom_ribbon(aes(ymin = f_low, ymax = f_high), alpha = 0.2, linewidth = 0.2) +
    geom_point(data = clonedf_end, aes(x = days, y = cloneVAF), shape = 4, size = 0.9) +
    #geom_linerange(data = clonedf_end, aes(x = days, ymin = cloneVAF - 2*0.05, ymax = cloneVAF + 2*0.05)) +
    scale_fill_manual(values = cp) +
    scale_color_manual(values = cp) +
    ylab("Clone frequency (simulations)") +
    xlab("Days since first surgery") +
    scale_y_continuous(limits = c(0, 1))
  
  ca125_patient <- filter(ca125_patient, days <= max(clonedf_end$days))
  
  ca125plot <- ca125_patient %>% 
    ggplot(aes(x = days, y = CA125)) +
    geom_line() +
    geom_point() +
    xlab("Days since first surgery")
  
  popplot <- ca125_patient %>% 
    ggplot(aes(x = days, y = population)) +
    geom_line() +
    geom_point() +
    xlab("Days since first surgery") +
    scale_y_log10() +
    ylab("Population size")
  
  pvalplot <- pvals %>% 
    ggplot(aes(x = clone, y = pval, fill = clone, col = clone)) +
    geom_col() +
    scale_fill_manual(values = cp) +
    scale_color_manual(values = cp) +
    logyscale() +
    ylab("p-value (consistent with WF model)") +
    geom_hline(yintercept = 0.05, lty = 2)
  
  myplot <- plot_grid(ca125plot, popplot,
                      wfplot, pvalplot,align = "v", axis = "tb",
                      ncol = 4, rel_widths = c(1,1,1.5,0.5)
                      )
  
  myplotlist <- list(ca125plot, popplot,
                     wfplot, pvalplot)
  
  print(myplot)
  
  Sys.sleep(3)
  
  return(list(plot = myplot, plotlist = myplotlist, pvals = pvals, ca125 = ca125, pval_combined = combine_p_values_fisher(pvals$pval)))
}
