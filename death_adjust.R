dat <- read.csv("death_pop.csv")
dat <- dat[1:85,]

parse_num <- function(x){
  x <- as.character(x)
  x <- gsub(',', '', x)
  x <- gsub(' ', '', x)
  as.numeric(x)
}

deaths <-parse_num(dat$Deaths)
pop <- parse_num(dat$Population)

#probability that a person dies at each age 
prob_die <- deaths / pop
prob_die <- prob_die / sum(prob_die)
prob_die <- c(prob_die,rep(0,100))

#probability that a person is dead at each age
prob_dead <- cumsum(prob_die)


#probability that a person is each age
prob_age <- pop / sum(pop)
age <- 0:84

# diagnosis probability unadjusted for deaths over 85 years
diag_prob <- c(0.0144,0.0103,0.0073,0.0052,0.0037,0.0026,0.0019,0.0013)
diag_prob <- c(diag_prob,rep(.0009, 85-8))

# probability of diagnosis interval and not dead
adjust_lag_prob <- function(diag_prob, pop, prob_die){
  max_age <- length(pop)
  new_prob <- rep(0, max_age)
  for(i in 1:(max_age-1)){
    pd <- prob_die[(i+1):(i+max_age)]
    pd <- pd / sum(pd)
    pdead <- cumsum(pd)
    new_prob <- new_prob + prob_age[i] * (1 - pdead) * diag_prob
  }
  new_prob
}


max_age <- length(pop)
prob_die_at_interval <- rep(0, max_age)
for(i in 1:max_age){
  pd <- prob_die[(i+1):(i+max_age)]
  if(sum(pd) == 0) 
    pd <- 0
  else
    pd <- pd / sum(pd)
  prob_die_at_interval <- prob_die_at_interval + prob_age[i] * pd
}



lag_probs <- adjust_lag_prob(diag_prob, pop, prob_die)

#probability die before diagnosis
1 - sum(lag_probs)


