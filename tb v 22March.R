# if (!require("devtools")) install.packages("devtools")
# devtools::install_github('hivbackcalc/package1.0/HIVBackCalc', 
#                          build_vignettes=TRUE)
# library(HIVBackCalc)
# library(ggplot2)

rm(list=ls())

diagInterval  <- 1
lag <- 50 

# importing the distribution of time from (Latent TB) infection to diagnosis (Active TB) 
library(xlsx)
pids_source <- read.xlsx2("~/Box Sync/CAPE - LTBI/Data/TB Data.xlsx", sheetIndex = "To R", startRow = 1,endRow = lag+1, colIndex = 2:4)
pids_source[,1]<-as.numeric(paste(pids_source[,1]))
pids_source[,2]<-as.numeric(paste(pids_source[,2]))
pids_source[,3]<-as.numeric(paste(pids_source[,3]))

# Since only 8.60% (95%CI 5.18 to 12.01) of Latent TB infection progressed to Active TB in 50 years, 
# we normalize the risk to 100% by dividing the risk for each year to the sum of 8.60% (and 5.18 12.01 for LL and UL).  
# repeat using pids_point, pids_LL and pids_UL
pids <- pids_source[,1] / sum(pids_source[,1])  # import  Point Est.
pids <- pids_source[,2] / sum(pids_source[,2])  # use for LL Est.
pids <- pids_source[,3] / sum(pids_source[,3])  # use for UL Est.

# Use "AreaCode<-3" for U.S. TB data / Use "AreaCode<-15" for CA TB data
# US 3,California	15, Alameda	8, Los Angeles	9, Orange	10, Sacramento	11, San Diego	12, San Francisco	13, Santa Clara	14
# California Non-US born 17, California US born 18
AreaCode<-3
counts_source <- read.xlsx2("~/Box Sync/CAPE - LTBI/Data/TB Data.xlsx", sheetIndex = "Reported TB Cases", startRow = 1,endRow = 25, colIndex = AreaCode)

# add 50 empty cells to the beginning of reported cases / required for calculation
counts <- c(rep(NA, 50), as.numeric(paste(as.vector(counts_source[,1]))))
 
names(counts) <- 2017 - length(counts):1

# make empty data frames to save model output
Back_est <-data.frame(matrix(nrow=length(counts)+50, ncol=3))
Back_est$a <- NULL
Back_LL<-Back_est
Back_UL<-Back_est

# make a loop to do the calculation for Point, LL and UL LTBI incidence 
for (K in 1:3) {
  
  pids <- pids_source[,K] / sum(pids_source[,K])  
  

pidFunc <- function (i) {
  ifelse(i > length(pids) - 1, 0, pids[i + 1])
}

# minor edit for nice plotting
plot.backproj <- function (x, time, showDiagCounts = TRUE, case = "", ...) {
  obs <- !is.na(x$y)
  if (missing(time)) 
    time <- as.numeric(names(x$y)[obs])
  else time <- seq(from = time[1], to = time[2], length.out = sum(obs))
  plot(time, x$lambda[obs], ylim = c(0, max(x$lambda[obs],x$y[obs])*1.1), type = "l", main = paste("Estimated Incidence", 
                                                                       case, sep = "\n"), ylab = "Count", ...)
  if (showDiagCounts) 
    points(time, x$y[obs], col = "red")
}

#slightly improved initial lambda
estimateIncidence <- function (y, pid, gamma = 0, tol = 10^-5, verbose = FALSE) 
{
  isna <- is.na(y)
  lambda <- y
  impute <- mean(y, na.rm = TRUE)
  for(i in length(y):1){
    if(!isna[i])
      impute <- y[i]
    else
      lambda[i] <- impute
  }
  #lambda <- rep(mean(y, na.rm = TRUE), length(y))
  
  ll <- lambda
  dev <- Inf
  while (dev > tol) {
    lambda <- meanEmUpdate(y, pid, lambda, gamma)
    dev <- sum((ll - lambda)^2/ll)
    if (verbose) {
      cat("lambda: ", paste(round(lambda, 1), collapse = " "), 
          "\n", "parameter change: ", dev, "\n", sep = "")
    }
    ll <- lambda
  }
  mod <- list(lambda = lambda, y = y, pid = pid, gamma = gamma, 
              tol = tol)
  class(mod) <- "backproj"
  mod
}

incidence2<- estimateIncidence(y=counts,
                               pid=pidFunc,
                               gamma=0.001,
                               verbose=TRUE)
plot(incidence2)

if (K==1) {
    Back_est<-data.frame(year = names(incidence2$y), count =counts, incidence=incidence2$lambda)
   } else if (K==2) {
     Back_UL<-data.frame(year = names(incidence2$y), count =counts, incidence_UL=incidence2$lambda)
  } else if (K==3) {
    Back_LL<-data.frame(year = names(incidence2$y), count =counts, incidence_LL=incidence2$lambda)
  }
}

Back<-data.frame(year = Back_est[,1], count = Back_est[,2], incidence = Back_est[,3], incidence_LL= Back_LL[,3], incidence_UL= Back_UL[,3])

require(ggplot2)
BackLTBI<- ggplot(data = Back[c(51:74),], aes(x = Back[c(51:74),1], y = Back[c(51:74),3]), group = 1) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = Back[c(51:74),5], ymin = Back[c(51:74),4])) +
  labs(x="Year", y="Incident LTBI #") +
  geom_line(size = 1, color = "RED",  aes(x= Back[c(51:74),1], y= Back[c(51:74),2]), data = Back[c(51:74),], group = 2)

# View the graph
BackLTBI

#### Calculation the # LTBI undiagnosed

# Import Mortality Rate and Popualtion
library(xlsx)
PopDeaths <- read.xlsx2("~/Box Sync/CAPE - LTBI/Data/TB Data.xlsx", sheetIndex = "Pop Deaths", startRow = 1,endRow = 13, colIndex = 8:13)
DeathRate<-as.numeric(paste(PopDeaths$DeathRate[which(PopDeaths$AreaCode==AreaCode)]))
Pop<-as.numeric(paste(PopDeaths$Pop[which(PopDeaths$AreaCode==AreaCode)]))
Pop6plus<-as.numeric(paste(PopDeaths$Pop6plus[which(PopDeaths$AreaCode==AreaCode)]))

# Number of total LTBI at each year
Back$TotalLTBI<-Back$incidence/sum(pids_source[,1])  # Point Est.
Back$TotalLTBI_LL<-Back$incidence_LL/sum(pids_source[,3])  # Lower limit 95%CI
Back$TotalLTBI_UL<-Back$incidence_UL/sum(pids_source[,2])  # Upper limit 95%CI

# Number of total LTBI alive at each year
Back$TotalLTBIalive<-Back$TotalLTBI-Back$TotalLTBI*DeathRate  # Point Est.
Back$TotalLTBIalive_LL<-Back$TotalLTBI_LL-Back$TotalLTBI_LL*DeathRate # Lower limit 95%CI
Back$TotalLTBIalive_UL<-Back$TotalLTBI_UL-Back$TotalLTBI_UL*DeathRate  # Upper limit 95%CI

# Number of total LTBI alive after 2016-Lag

LTBI<-sum(Back$TotalLTBIalive[which(as.numeric(as.character(Back$year))>=2016-lag+1)])
LTBI_LL<-sum(Back$TotalLTBIalive_LL[which(as.numeric(as.character(Back$year))>=2016-lag+1)])
LTBI_UL<-sum(Back$TotalLTBIalive_UL[which(as.numeric(as.character(Back$year))>=2016-lag+1)])


# Prevalence of LTBI in 2016 
Prv_LTBIalive<-LTBI/Pop
Prv_LTBIalive_LL<-LTBI_LL/Pop
Prv_LTBIalive_UL<-LTBI_UL/Pop

# Prevalence of LTBI in 2016 among Age 6+
Prv_LTBIalive6Plus<-LTBI/Pop6plus
Prv_LTBIalive6Plus_LL<-LTBI_LL/Pop6plus
Prv_LTBIalive6Plus_UL<-LTBI_UL/Pop6plus

Back_Output<-data.frame(year = Back[,1], count = Back[,2], "New LTBI ultimately diagnosed during the lag time" = Back[,3], "LL New LTBI ultimately diagnosed during the lag time"= Back[,4], "UL New LTBI ultimately diagnosed during the lag time"= Back[,5],"Total New LTBI alive"= Back[,9],"LL Total New LTBI Alive"= Back[,10],"UL Total New LTBI Alive"= Back[,11])
Prv_Output<-data.frame(Output="Total LTBI living in 2016","Area Code"=AreaCode,Point = LTBI, LowerLimit = LTBI_LL, UpperLimit= LTBI_UL)
Prv_Output2<-data.frame(Output="Prevalence of LTBI in 2016","Area Code"=AreaCode, Point = Prv_LTBIalive, LowerLimit = Prv_LTBIalive_LL, UpperLimit= Prv_LTBIalive_UL)
Prv_Output3<-data.frame(Output="Prevalence of LTBI among Age 6+ in 2016","Area Code"=AreaCode, Point = Prv_LTBIalive6Plus, LowerLimit = Prv_LTBIalive6Plus_LL, UpperLimit= Prv_LTBIalive6Plus_UL)

Prv_Output<-rbind(Prv_Output,Prv_Output2,Prv_Output3)

# change working directory
setwd("~/Box Sync/CAPE - LTBI/R")

# save the graph
pdf("BackLTBI.pdf")
plot(BackLTBI)
dev.off()

# Make Table for Incident LTBI for 1967 to 2016
library(xlsx)
write.xlsx(Back_Output[which(as.numeric(Back_Output$year)>=25),], "LTBI Total Output.xlsx")
write.xlsx(Prv_Output, "LTBI Prv Output.xlsx")

# Make Table for Incident LTBI for 1967 to 2016
library(xlsx)
write.xlsx(Back[which(as.numeric(Back$year)>=25),], "Back Output.xlsx")

