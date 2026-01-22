# Load Libraries ----------------------------------------------------------

pacman:: p_load(
  tidyverse,
  deSolve,
  plotly,
  gridExtra,
  reshape2,
  dplyr,
  FME
)



# define equations --------------------------------------------------------

#D:Define the equations
#I:Initialize the compartments
#P: Define parameters
#T: Define time interval
#O : solve the equations 


# seasonality -------------------------------------------------------------


# Function to calculate seasonality based on rainfall and phase angle
calculate_seasonality <- function(times, rainfall_data, amp, phi) {
  # Create a time series of rainfall data
  rainfall_ts <- ts(rainfall_data$rainfall, frequency = 12)
  
  # Extend the time series to cover the entire simulation period
  years_needed <- ceiling(max(times) / 365)+ 1
  extended_rainfall <- rep(as.vector(rainfall_ts), years_needed)
  
  # Smooth the rainfall data using a moving average
  smoothed_rainfall <- zoo::rollmean(extended_rainfall, k = 3, fill = "extend")
  
  # Interpolate to get daily values
  daily_rainfall_func <- approxfun(seq_along(smoothed_rainfall) * (365/12), smoothed_rainfall)
  daily_rainfall <- daily_rainfall_func(times)
  
  # Calculate seasonality factor
  max_rainfall <- max(rainfall_data$rainfall)
  min_rainfall <- min(rainfall_data$rainfall)
  rainfall_range <- max_rainfall - min_rainfall
  
  # Incorporate phase angle into the seasonality calculation
  seasonality <- 1 + amp * (daily_rainfall - min_rainfall) / rainfall_range *
    cos(2 * pi * (times / 365.25 - phi))
  
  return(seasonality)
} 
# Get seasonality factor for current time point
seas_vec <- seasonality[t + 1]  # t + 1 because R indices start at 1



# model equations ---------------------------------------------------------


chikungunya_model <- function(
    t,y,parms){
  with(as.list(c(parms,y)),{
    
    idx <- pmin(length(seasonality_vec), pmax(1, floor(t) + 1))
    seas <- seasonality_vec[idx]
    #Age-specific compartments
    S1 <- y[1]; E1 <- y[2]  #Age Group 1 (<60 years): Susceptible, Exposed
    S2 <- y[3]; E2 <- y[4]  #Age group 2 (60+ years): Susceptible, Exposed
    
    #Age Independent compartments
    Im <- y[5] #Infectious (mild)
    Is <- y[6] #Infectious(severe)
    C <-  y[7] #Chronic
    R <-  y[8] #Recovered
    #V <- y[9] #Vaccinated
    
    #Mosquito compartments
    Sm <- y[9]; Em <- y[10]; Im_mosq <- y[11]
    
    #Total human populations by age group 
    N1 = S1 + E1 + prop_age1_Im * Im + prop_age1_Is * Is + prop_age1_C* C + prop_age1_R *R 
    N2 = S2 + E2 + (1- prop_age1_Im)*Im + (1- prop_age1_Is)*Is + (1- prop_age1_C)* C + (1- prop_age1_R) *R
    N_total <- N1 + N2
    
    #Total infectious humans 
    I_total <- Im +Is
    
    # Force of infection from mosquitoes to humans
    lambda_h <- b * beta_mh * Im_mosq / N_total
    
    # Force of infection from humans to mosquitoes
    lambda_m <- b * beta_hm * I_total/ N_total
    
    # Age-specific forces of infection
    lambda_h1 <- lambda_h * contact_rate_1 *seas
    lambda_h2 <- lambda_h * contact_rate_2 *seas
    
    #Age group 1 (<60 years)
    dS1 = mu *N1 - mu * S1 - lambda_h1 *S1 #mu; birth/death rate
    dE1 = lambda_h1* S1 - (sigma+mu)* E1
    
    #Age group 2 (60+ years)
    dS2 = mu * N2 - mu* S2 - lambda_h2 *S2
    dE2 = lambda_h2 *S2 - (sigma + mu)*E2
    
    #Infectious mild~ receives from both age groups based on severity probabilities
    dIm = (sigma * (1- severe_prob_1)* E1) + (sigma * (1 - severe_prob_2)* E2) - (mu + gamma_m + mu_ms)*Im
    
    #Infectious severe~ receives from both age groups based on severity probabilities
    dIs = (sigma * severe_prob_1 * E1)+ (sigma * severe_prob_2 * E2)+ (mu_ms * Im)- (mu + gamma_s + chronic_prob_s* gamma_s)* Is
    
    #Chronic 
    dC = (chronic_prob_s * gamma_s * Is)- (mu + gamma_c)* C
    
    #Recovered - receives from mild, severe, and chronic compartments
    dR = (gamma_m *Im)+ ((1- chronic_prob_s)* gamma_s *Is)+ (gamma_c * C)- (mu * R)
    
    #Mosquito equations
    dSm = Lambda_m  - mu_m * Sm - lambda_m * Sm
    dEm = lambda_m * Sm - (sigma_m + mu_m)* Em
    dIm_mosq = sigma_m * Em - mu_m * Im_mosq
    
    #Add dCumInc too show incidences
    #dCumInc = (sigma* E1) + (sigma * E2)
    dCum1 = sigma * E1
    dCum2 = sigma * E2
    
    # Return 11 derivatives for the ode solver
    list(c(dS1, dE1,dS2,dE2,dIm, dIs, dC, dR, dSm,dEm, dIm_mosq, dCum1, dCum2))
    
  })
}


# function to run simulation ----------------------------------------------

run_chikungunya_simulation <- function (params, initial_conditions, time_points){
  
  #Solve ODEs
  out<- ode(y= initial_conditions,
            times= time_points,
            func = chikungunya_model,
            parms = params,
            method = "lsoda")
  
  # Convert to dataframe
  out_df <- as.data.frame(out)
  #browser()
  
  #Add column names
  colnames(out_df) <- c ("time", "S1", "E1", "S2", "E2", "Im", "Is", "C", "R", "Sm", "Em", "Im_mosq", "Cum1", "Cum2")
  
  return(out_df)
}


# Function to plot results ------------------------------------------------

plot_chikungunya_results <- function(results){
  
  #Calculate derived quantities
  results$Total_Infectious <- results$Im + results$Is
  results$Total_Susceptible <- results$S1 + results$S2
  results$Total_Exposed <- results$E1 + results$E2
  results$Total_Mosquito <- results$Sm + results$Em + results$Im_mosq
  
  # Plot 1: Age-specific Susceptible and exposed
  age_specific_data <- results |>
    select(time, S1, E1,S2, E2) |>
    melt(id.vars = "time", variable.name = "Compartment", value.name = "Population")
  
  age_specific_data$Age_Group <- ifelse(grepl("1", age_specific_data$Compartment),
                                        "Age group 1 (<60)", "Age group 2 (60+)")
  age_specific_data$State <- ifelse(grepl("S", age_specific_data$Compartment),
                                    "Susceptible", "Exposed")
  
  p1 <- ggplot(age_specific_data, aes(x= time, y= Population, color = State))+
    geom_line(linewidth = 1.2)+
    facet_wrap(~Age_Group, scales ="free_y")+
    labs(title = "Age-Specific Susceptible and Exposed Populations",
         x= "Time(days)", y = "Population")+
    theme_minimal()+
    theme(legend.position = "bottom")
  
  
  #Plot 2: Age-independent disease states
  disease_data <- results |>
    select(time, Im, Is, C, R)|>
    melt(id.vars ="time",variable.name = "Compartment", value.name = "Population")
  
  disease_data$Compartment <- factor(disease_data$Compartment,
                                     levels = c("Im", "Is", "C","R"),
                                     labels = c("Infectious(mild)", "Infectious(severe)",
                                                "Chronic", "Recovered"))
  
  p2 <- ggplot(disease_data, aes(x= time, y = Population, color = Compartment))+
    geom_line(linewidth= 1.2)+
    labs(title = "Disease States (Age-Independent)",
         x = "Time (days)", y = "Population")+
    theme_minimal()+
    theme(legend.position = "bottom")
  
  #Plot 3: Mosquito population
  mosquito_data <- results |>
    select(time, Sm, Em, Im_mosq)|>
    melt(id.vars = "time", variable.name = "Compartment", value.name = "Population")
  
  mosquito_data$Compartment <- factor(mosquito_data$Compartment,
                                      levels = c("Sm","Em", "Im_mosq"),
                                      labels = c("Susceptible", "Exposed", "Infectious"))
  
  p3 <- ggplot(mosquito_data, aes(x= time, y = Population, color = Compartment)) +
    geom_line(linewidth = 1.2)+
    labs(title = "Mosquito Population",
         x= "Time(days)", y = "Population")+
    theme_minimal()+
    theme(legend.position = "bottom")
  
  #Plot 4: Key Epidemiological Indicators
  epi_data <- results |>
    select(time, Total_Infectious, C, Total_Susceptible)|>
    melt(id.vars ='time', variable.name = "Indicator", value.name = "Population")
  
  
  epi_data$Indicator <- factor(epi_data$Indicator,
                               levels = c("Total_Infectious", "C", "Total_Susceptible"),
                               labels = c("Total Infectious", "Chronic", "Total Susceptible"))
  
  p4 <- ggplot(epi_data, aes(x= time, y = Population, color = Indicator))+
    geom_line(linewidth= 1.2)+
    labs(title = "Key Epidemiological Indicators",
         x= "Time(days)", y= "Population")+
    theme_minimal()+
    theme(legend.position = "bottom")
  
  #Arrange plots 
  grid.arrange(p1,p2,p3,p4,ncol=2)
  
} 

# Function to calculate RO ------------------------------------------------

calculate_RO <- function(params){
  with(params,{
    # Average infectious period considerinng chronic progression
    avg_infectious_period_mild <- 1 /(gamma_m + mu)
    avg_infectious_period_severe <- 1 /(gamma_s + mu)
    
    # Weighted average based on severity probabilities
    avg_severe_prob <- (severe_prob_1 + severe_prob_2)/2
    avg_infectious_period <- (1 - avg_severe_prob)* avg_infectious_period_mild + 
      avg_severe_prob * avg_infectious_period_severe
    
    avg_infectious_period_mosquito <- 1 / mu_m
    prob_transmission <- (sigma_m / ( sigma_m + mu_m))
    
    RO <- (b^2 * beta_hm * beta_mh * prob_transmission * avg_infectious_period * avg_infectious_period_mosquito)/mu_m
    
    return(RO)
  })
}

# Function to calculate cumulative incidence by age group -----------------

calculate_cumulative_incidence <- function(results, params) {
  
  # Calculate new infections over time (approximation)
  results$new_infections_1 <- c(0, diff(results$E1)) * params$sigma
  results$new_infections_2 <- c(0, diff(results$E2)) * params$sigma
  
  # Cumulative incidence
  results$cum_incidence_1 <- cumsum(pmax(0, results$new_infections_1))
  results$cum_incidence_2 <- cumsum(pmax(0, results$new_infections_2))
  
  return(results)
}


# Main execution function -------------------------------------------------

main_simulation <- function() {
  
  # 1. Define time points and rainfall data FIRST
  time_points <- seq(0, 365*5, by = 1) 
  
  rainfall_data <- data.frame(
    month = 1:12,
    rainfall = c(0.38, 1.52, 6.49, 29.32, 64.98, 109.02, 176.81, 240.86, 148.25, 49.59, 3.58, 0.25)
  )
  
  # 2. Calculate seasonality values SECOND
  seasonality_values <- calculate_seasonality(
    times = time_points, 
    rainfall_data = rainfall_data, 
    amp = 0.7, 
    phi = 0.33
  )
  
  # 3. Define parameters list THIRD (including the seasonality vector we just made)
  params <- list(
    # Human parameters
    mu = 1/(62*365),
    sigma = 1/5,
    gamma_m = 1/7,
    gamma_s = 1/14,
    gamma_c = 1/60,
    severe_prob_1 = 0.1,
    severe_prob_2 = 0.3,
    chronic_prob_s = 0.25,
    mu_ms = 1/10,
    contact_rate_1 = 1.0,
    contact_rate_2 = 0.8,
    
    # Seasonality vector is now safely defined
    seasonality_vec = seasonality_values,
    
    # Proportions
    prop_age1_Im = 0.8,
    prop_age1_Is = 0.7,
    prop_age1_C = 0.75,
    prop_age1_R = 0.8,
    
    # Mosquito parameters
    Lambda_m = 1000,
    mu_m = 1/14,
    sigma_m = 1/7,
    
    # Transmission
    b = 0.5,
    beta_hm = 0.1,
    beta_mh = 0.1
  )
  
  # 4. Initial conditions
  N1_initial <- 1108333
  N2_initial <- 100000
  
  initial_conditions <- c(
    S1 = N1_initial - 100, E1 = 0,
    S2 = N2_initial - 50,  E2 = 0,
    Im = 500, Is = 173, C = 0, R = 0,
    Sm = 9000, Em = 500, Im_mosq = 500,
    Cum1 = 0, Cum2 = 0
  )
  
  # 5. Run simulation
  cat("Running modified chikungunya simulation...\n")
  results <- run_chikungunya_simulation(params, initial_conditions, time_points)
  
  # Calculate cumulative incidence
  results <- calculate_cumulative_incidence(results, params)
  
  # Calculate RO
  RO <- calculate_RO(params)
  cat(sprintf("Basic reproduction number (RO): %.3f\n", RO))
  
  # Plot results
  cat("Generating plots...\n")
  plot_chikungunya_results(results)
  
  return(list(results = results, params = params, RO = RO))
}


# Function for chronic disease burden disease analysis --------------------

analyze_chronic_burden <- function(results){
  
  # Peak Chronic cases
  peak_chronic <- max(results$C)
  peak_chronic_time <- results$time[which.max(results$C)]
  
  #Chronic prevalence over time 
  chronic_data <- data.frame(
    time = results$time,
    chronic_cases = results$C,
    chronic_prevalence = results$C / (results$S1 + results$E1 + results$S2 + results$E2 +
                                        results$Im + results$Is + results$C + results$R)
  )
  
  # Plot chronic burden 
  p_chronic <- ggplot(chronic_data, aes(x= time))+
    geom_line(aes(y= chronic_cases, color = "Chronic Cases"), linewidth = 1.2)+
    geom_line(aes(y = chronic_prevalence* 1000 , color = "Prevalence (per 1000)"),linewidth= 1.2)+
    scale_y_continuous(
      name = "Chronic Cases",
      sec.axis = sec_axis(~./1000, name = "Prevalence (per 1000)")
    )+ labs(title = "Chronic Chikungunya Burden Over Timme",
            x= "Time(days)", color = "Measure")+
    theme_minimal()+
    theme(legend.position = "bottom")
  
  print(p_chronic)
  
  cat(sprintf("Peak chronic cases: %.0f at day %.0f\n", peak_chronic, peak_chronic_time))
  
  return(chronic_data)
  
}



# Run the main simulation -------------------------------------------------

if(interactive()){
  simulation_results <- main_simulation()
  cat("Simulation completed successfully!\n")
  
  #Analyze chronic disease burden 
  cat("\nAnalyzing chronic disease burden...\n")
  chronic_analysis <- analyze_chronic_burden(simulation_results$results)
}






































