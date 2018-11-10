dat <-mortality

AreaCode2<-AreaCode
SubPop2<-SubPop

# dat <- read.csv("death_pop.csv")
dat <- dat[which(dat$AreaCode == AreaCode2 & dat$SubPop== SubPop2),]

deaths <-dat$Deaths
pop <- dat$Population

#probability that a person dies at each age 
prob_die <- deaths / pop
prob_die <- prob_die / sum(prob_die)
prob_die <- c(prob_die,rep(0,100))

#probability that a person is dead at each age
prob_dead <- cumsum(prob_die)


#probability that a person is each age
prob_age <- pop / sum(pop)
age <- 0:84

diag_prob <-c(0.0129, rep(0.0016,11), rep(0.0015,42), rep(0.0014,31))

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
  pd <- prob_die[(i):(i+max_age-1)]
  if(sum(pd) == 0) 
    pd <- 0
  else
    pd <- pd / sum(pd)
  prob_die_at_interval <- prob_die_at_interval + prob_age[i] * pd
}

lag_probs <- adjust_lag_prob(diag_prob, pop, prob_die)

#probability die before diagnosis
1 - sum(lag_probs)


