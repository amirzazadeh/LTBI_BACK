rm(list = ls()) #not great to do this

# Use "AreaCode<-3" for U.S. TB data / Use "AreaCode<-15" for CA TB data
# US 3,California	15, Alameda	8, Los Angeles	9, Orange	10, Sacramento	11, San Diego	12, San Francisco	13, Santa Clara	14
# California Non-US born 17, California US born 18
AreaCode <- 3

# Granularity of diagnosis counts in years.
diagInterval  <- 1

# Maximum number of years to active
lag <- 85

source("death_adjust.R")
# if (!require("devtools")) install.packages("devtools")
# devtools::install_github('hivbackcalc/package1.0/HIVBackCalc',
#                          build_vignettes=TRUE)
library(HIVBackCalc)
library(ggplot2)
library(xlsx)

# minor edit for nice plotting
plot.backproj <-function (x,
                          time,
                          showDiagCounts = TRUE,
                          case = "",
                          ...) {
  obs <- !is.na(x$y)
  if (missing(time))
    time <- as.numeric(names(x$y)[obs])
  else
    time <- seq(from = time[1],
                to = time[2],
                length.out = sum(obs))
  plot(
    time,
    x$lambda[obs],
    ylim = c(0, max(x$lambda[obs], x$y[obs]) * 1.1),
    type = "l",
    main = paste("Estimated Incidence",
                 case, sep = "\n"),
    ylab = "Count",
    ...
  )
  if (showDiagCounts)
    points(time, x$y[obs], col = "red")
}

# importing the distribution of time from (Latent TB) infection to diagnosis (Active TB)
pids_source <-
  read.xlsx2(
    "TB Data.xlsx",
    sheetIndex = "To R",
    startRow = 1,
    endRow = lag + 1,
    colIndex = 2:4
  )
pids_source[, 1] <- as.numeric(paste(pids_source[, 1]))
pids_source[, 2] <- as.numeric(paste(pids_source[, 2]))
pids_source[, 3] <- as.numeric(paste(pids_source[, 3]))

# Import Mortality Rate and Popualtion
PopDeaths <-
  read.xlsx2(
    "TB Data.xlsx",
    sheetIndex = "Pop Deaths",
    startRow = 1,
    endRow = 13,
    colIndex = 8:13
  )
DeathRate <-
  as.numeric(paste(PopDeaths$DeathRate[which(PopDeaths$AreaCode == AreaCode)]))
Pop <-
  as.numeric(paste(PopDeaths$Pop[which(PopDeaths$AreaCode == AreaCode)]))
Pop6plus <-
  as.numeric(paste(PopDeaths$Pop6plus[which(PopDeaths$AreaCode == AreaCode)]))

# Since only 8.60% (95%CI 5.18 to 12.01) of Latent TB infection progressed to Active TB in 50 years,
# we normalize the risk to 100% by dividing the risk for each year to the sum of 8.60% (and 5.18 12.01 for LL and UL).
# repeat using pids_point, pids_LL and pids_UL
pids <- pids_source[, 1] / sum(pids_source[, 1])  # import  Point Est.
pids <- pids_source[, 2] / sum(pids_source[, 2])  # use for LL Est.
pids <- pids_source[, 3] / sum(pids_source[, 3])  # use for UL Est.

counts_source <-
  read.xlsx2(
    "TB Data.xlsx",
    sheetIndex = "Reported TB Cases",
    startRow = 1,
    endRow = 25,
    colIndex = AreaCode
  )

# add 50 empty cells to the beginning of reported cases / required for calculation
counts <-
  c(rep(NA, 60), as.numeric(paste(as.vector(counts_source[, 1]))))

names(counts) <- 2017 - length(counts):1

# make empty data frames to save model output
Back_est <- data.frame(matrix(nrow = length(counts) + 50, ncol = 3))
Back_est$a <- NULL
Back_LL <- Back_est
Back_UL <- Back_est


pactive_before_die <- c(0,0,0)
# make a loop to do the calculation for Point, LL and UL LTBI incidence
for (K in 1:3) {
  
  # probability that individual becomes active i years after infection and hasn't died during that period.
  lagp <- pids_source[, K]
  lagp[(length(lagp)+1):85] <- 0
  lagp <- adjust_lag_prob(lagp, pop, prob_die)
  pactive_before_die[K] <- sum(lagp)
  
  pids <- lagp / sum(lagp)
  #pids <- pids_source[, K] / sum(pids_source[, K])
  pidFunc <- function (i) {
    ifelse(i > length(pids) - 1, 0, pids[i + 1])
  }
  
  incidence2 <- estimateIncidence(
    y = counts,
    pid = pidFunc,
    gamma = 0.001,
    verbose = TRUE
  )
  plot(incidence2)
  
  if (K == 1) {
    Back_est <-
      data.frame(
        year = names(incidence2$y),
        count = counts,
        incidence = incidence2$lambda
      )
  } else if (K == 2) {
    Back_UL <-
      data.frame(
        year = names(incidence2$y),
        count = counts,
        incidence_UL = incidence2$lambda
      )
  } else if (K == 3) {
    Back_LL <-
      data.frame(
        year = names(incidence2$y),
        count = counts,
        incidence_LL = incidence2$lambda
      )
  }
}

Back <-
  data.frame(
    year = Back_est[, 1],
    count = Back_est[, 2],
    incidence = Back_est[, 3],
    incidence_LL = Back_LL[, 3],
    incidence_UL = Back_UL[, 3]
  )

require(ggplot2)
BackLTBI <-
  ggplot(data = Back[c(51:74), ], aes(x = Back[c(51:74), 1], y = Back[c(51:74), 3]), group = 1) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = Back[c(51:74), 5], ymin = Back[c(51:74), 4])) +
  labs(x = "Year", y = "Incident LTBI #") +
  geom_line(
    size = 1,
    color = "RED",
    aes(x = Back[c(51:74), 1], y = Back[c(51:74), 2]),
    data = Back[c(51:74), ],
    group = 2
  )

# View the graph
BackLTBI

#### Calculation the # LTBI undiagnosed

# Number of total infected at each year
Back$TotalInfectedInYear <- Back$incidence / pactive_before_die[1]  # Point Est.
Back$TotalInfected <- 0
infected2016 <- 0
index2016 <- length(Back$TotalInfectedInYear)
pid <- pids_source[,1]
for(i in 1:index2016){
  interval <- index2016 - i
  ind <- min(interval,length(pid))
  if(interval != 0){
    infected2016 <- infected2016 + Back$TotalInfectedInYear[i] * (1 - sum(pid[1:ind])) *
      (1 - sum(prob_die_at_interval[1:interval]))
    infected2016 <- infected2016 + .5 * Back$TotalInfectedInYear[i] * (pid[ind] +
      prob_die_at_interval[interval] - pid[ind] * prob_die_at_interval[interval] )
  }
}



Back$TotalInfectedInYear_UL <- Back$incidence_UL / pactive_before_die[2]  # Point Est.
Back$TotalInfected <- 0
infected2016_UL <- 0
index2016 <- length(Back$TotalInfectedInYear_UL)
pid <- pids_source[,2]
for(i in 1:index2016){
  interval <- index2016 - i
  ind <- min(interval,length(pid))
  if(interval != 0){
    infected2016_UL <- infected2016_UL + Back$TotalInfectedInYear_UL[i] * (1 - sum(pid[1:ind])) *
      (1 - sum(prob_die_at_interval[1:interval]))
    infected2016_UL <- infected2016_UL + .5 * Back$TotalInfectedInYear_UL[i] * (pid[ind] +
                                                                         prob_die_at_interval[interval] - pid[ind] * prob_die_at_interval[interval] )
  }
}

infected2016_UL / Pop


Back$TotalInfectedInYear_LL <- Back$incidence_LL / pactive_before_die[3]  # Point Est.
Back$TotalInfected <- 0
infected2016_LL <- 0
index2016 <- length(Back$TotalInfectedInYear_LL)
pid <- pids_source[,3]
for(i in 1:index2016){
  interval <- index2016 - i
  ind <- min(interval,length(pid))
  if(interval != 0){
    infected2016_LL <- infected2016_LL + Back$TotalInfectedInYear_LL[i] * (1 - sum(pid[1:ind])) *
      (1 - sum(prob_die_at_interval[1:interval]))
    infected2016_LL <- infected2016_LL + .5 * Back$TotalInfectedInYear_LL[i] * (pid[ind] +
                                                                                  prob_die_at_interval[interval] - pid[ind] * prob_die_at_interval[interval] )
  }
}

infected2016_LL / Pop


df <- data.frame(Code=AreaCode,Est=infected2016/Pop, LL=infected2016_LL / Pop, UL=infected2016_UL / Pop)
df

write.csv(df,
          file=paste0("CurrentLatentInfections_",AreaCode,".csv"))


