
pacman::p_load(tidyverse,
               dplyr)
# set parameters
parms=c(m=1e-4, b=0.02, v=0.1, r=0.3)
initial=c(S=50, I=1, R=0)
time.window=c(0, 100)

# initialize state and time variables
state <- initial
time <-  time.window[1]

# define output dataframe
output <- data.frame(t=time, 
                     S=state["S"], I=state["I"], R=state["R"])

# define how state variables S, I and R change for each process
# Each row is a transition: Birth, Death S, Infection, Death I, Recovery, Death R
processes <- matrix(c(
  1,  0,  0,  # Birth: S increases
  -1,  0,  0,  # Death S: S decreases
  -1,  1,  0,  # Infection: S -> I
  0, -1,  0,  # Death I: I decreases
  0, -1,  1,  # Recovery: I -> R
  0,  0, -1   # Death R: R decreases
), nrow=6, ncol=3, byrow=TRUE,
dimnames=list(c("birth", "death.S", "infection", "death.I", "recovery", "death.R"),
              c("S", "I", "R")))

# process probabilities (propensities)
probabilities <- function(state, parms){
  N <- sum(state)
  with(as.list(c(state, parms)), {
    return(c(
      birth     = m * N,
      death.S   = m * S,
      infection = b * S * I,
      death.I   = (m + v) * I,
      recovery  = r * I,
      death.R   = m * R
    ))
  })
}

while(time < time.window[2] & state["I"] > 0){
  
  # 1. Calculate process rates (propensities)
  rates <- probabilities(state, parms)
  total_rate <- sum(rates)
  
  if(total_rate == 0) break # Stop if no more events can happen
  
  # 2. WHEN does the next process happen? (Draw from exponential distribution)
  tau <- rexp(1, rate = total_rate)
  
  # 3. Update time
  time <- time + tau
  
  # 4. WHICH process happens? (Draw based on relative probabilities)
  event <- sample(1:nrow(processes), size = 1, prob = rates)
  
  # 5. Update states using the processes matrix
  state <- state + processes[event, ]
  
  # 6. Write into output
  output <- rbind(output, c(t=time, state))
}

# View results
print(head(output))



# Convert to long format for ggplot
output_long <- pivot_longer(output, 
                            cols = c("S", "I", "R"), 
                            names_to = "Compartment", 
                            values_to = "Count")

# Ensure the order of compartments is logical
output_long$Compartment <- factor(output_long$Compartment, levels = c("S", "I", "R"))



ggplot(output_long, aes(x = t, y = Count, color = Compartment)) +
  geom_step(linewidth = 1) +
  scale_color_manual(values = c("S" = "#2c7bb6", "I" = "#d7191c", "R" = "#1a9641")) +
  labs(title = "Stochastic SIR Model (Gillespie Algorithm)",
       subtitle = "Population dynamics with birth, death, and recovery",
       x = "Time",
       y = "Number of Individuals",
       color = "Status") +
  theme_minimal() +
  theme(legend.position = "bottom")
