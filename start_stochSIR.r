# SIMULATING A STOCHASTIC SIR MODEL USING A GILLESPIE ALGORITHM
# Description:
# 	parms: the parameters in the SIR model:
#   	dS/dt = m*(S+I+R) - m*S - b*S*I
#   	dI/dt = b*S*I - (m+v)*I - r*I
#   	dR/dt = r*I - m*R
#
# 	initial: are the initial values for the variables S, I and R
#
# 	time.window: specifies the start and end time of the simulation
#
# Output: a data.frame that start approximately like
#            t  S  I  R
#1   0.0000000 50  1  0
#2   0.4202398 49  2  0
#3   0.5046732 48  3  0
#    ...       ..  .  .

# set parameters
parms=c(m=1e-4,b=0.02,v=0.1,r=0.3)
initial=c(S=50, I=1, R=0)
time.window=c(0, 100)
    
# initialize state and time variables and write them into output
state <- initial
time <-  time.window[1]
  
# define output dataframe
output <- data.frame(t=time,
                     S=state["S"], I=state["I"], R=state["R"],
                     row.names=1)
  
# define how state variables S, I and R change for each process
processes <- matrix(0, nrow=6, ncol=3,
                    dimnames=list(c("birth",
                        "death.S",
                        "infection",
                        "death.I",
                        "recovery",
                        "death.R"),
                        c("dS","dI","dR")))
<...>
                    
# process probabilities
probabilities <- function(state){
  <...>
}


while(time < time.window[2] & state["I"]>0){

  # calculate process probabilities for current state
  <...>

  # WHEN does the next process happen?
  <...>

  # update time
  <...>

  # WHICH process happens after tau?
  <...>

  # update states
  <...>

  # write into output
  output <- rbind(output,c(time,state))
}

output