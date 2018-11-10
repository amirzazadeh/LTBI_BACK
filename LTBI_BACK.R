# change working directory to your choice
setwd("~/LTBI_Back")

# Code and installation instructions: https://github.com/hivbackcalc/package1.0
# install.packages("devtools")
# devtools::install_github('hivbackcalc/package1.0/HIVBackCalc')

library(HIVBackCalc)
library(ggplot2)

# Import and clean all data from "TB data" excel file
source("01_tidy.R")

# Choose the subpopualtion for your analysis: "Total"
SubPop <-"Total"

# Choose the first and last year of TB reproted cases 

firsty<-2000 # first year can not be lower than the first year of TB reported cases
lasty<-2016 # last year can not exceed the last year of TB reported cases
nyear<-lasty-firsty

# import the area codes for analysis
A<-totalpop$AreaCode[which(totalpop$SubPop==SubPop)]
AA<-0

## loop over the area code
for (S in A) {
AreaCode <- S

# Granularity of diagnosis counts in years.
diagInterval  <- 1

# Maximum number of years to active
lag <- 85

# function to adjust lag distribution for mortality (competing) risk
source("death_adjust.R")

# Import lag distribution (from Latent TB infection to Active TB)
pids_source <- react[,c(2:4)]

# we normalize the risk to 100% by dividing the risk for each year to the sum.
pids <- pids_source[, 1] / sum(pids_source[, 1])  # use for Point Est.
pids <- pids_source[, 2] / sum(pids_source[, 2])  # use for UL Est.
pids <- pids_source[, 3] / sum(pids_source[, 3])  # use for LL Est.

# Import Total Popualtion Size
TotPop <- totalpop$Pop[which(totalpop$AreaCode == AreaCode & totalpop$SubPop== SubPop)]

# Import reported annual active TB cases
counts_source <- cases$TB[which(cases$AreaCode == AreaCode & cases$SubPop== SubPop & (cases$Year>=firsty & cases$Year<=lasty))]

# add 100 empty cells to the beginning of reported cases / required for calculation
counts <- c(rep(NA, 100+nyear), counts_source[])
names(counts) <- (lasty+1) - length(counts):1

# make empty dataframes to save model's outputs
Back_est <- data.frame(matrix(nrow = length(counts) + 100+nyear, ncol = 3))
Back_est$a <- NULL
Back_LL <- Back_est
Back_UL <- Back_est


pactive_before_die <- c(0,0,0)
# make a loop to do the calculation for Point, LL and UL LTBI incidence
for (K in 1:3) {

   # probability that individual becomes active i years after infection and hasn't died during that period.
  lagp <- pids_source[, K]
  if(length(lagp) < 85)
    lagp[(length(lagp)+1):85] <- 0
  lagp <- adjust_lag_prob(lagp, pop, prob_die)
  pactive_before_die[K] <- sum(lagp)
  
  # store Adj Lag
  if (K == 1) {
    lagp_adj <-lagp
    } else if (K == 2) {
      lagp_adjLL <-lagp
    } else if (K == 3) {
      lagp_adjUL <-lagp
    }
  
  # normalize the risk to 100% 
  pids <- lagp / sum(lagp)
  
  # calculate LTBI incidence in each year 
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

  # store LTBI incidence in each year
  if (K == 1) {
    Back_est <-data.frame(year = names(incidence2$y),count = counts,incidence = incidence2$lambda)
    } else if (K == 2) {
      Back_UL <-data.frame(year = names(incidence2$y),count = counts,incidence = incidence2$lambda)
       } else if (K == 3) {
         Back_LL <-data.frame(year = names(incidence2$y),count = counts,incidence = incidence2$lambda)
   }
}


# store Active TB, LTBI Point, LL, and UL Cases 
Back <-
  data.frame(
    year = Back_est[, 1],
    count = Back_est[, 2],
    incidence = Back_est[, 3],
    incidence_LL = Back_LL[, 3], 
    incidence_UL = Back_UL[, 3] 
  )


#### Calculation the # LTBI undiagnosed

# Number of total infected at each year
Back$TotalInfectedInYear <- Back$incidence / pactive_before_die[1]  # Point Est.

infected2016 <- 0
index2016 <- length(Back$TotalInfectedInYear)
infectedAliveInYear<-c(0,rep(0,index2016-1))
pid <- pids_source[,1]

prob_die_at_interval<-c(prob_die_at_interval,rep(0,index2016-85))
 
for(i in 1:index2016){
  interval <- index2016 - i
  ind <- min(interval,length(pid))
  if(interval != 0){
    infected2016 <- infected2016 + Back$TotalInfectedInYear[i] * (1 - sum(pid[1:ind])) *
      (1 - sum(prob_die_at_interval[1:interval]))
    infected2016 <- infected2016 + .5 * Back$TotalInfectedInYear[i] * pid[ind] * (1 - sum(prob_die_at_interval[1:interval]) + .5 * prob_die_at_interval[interval] )
    }
  infectedAliveInYear[i]<-infected2016
  }

infected2016 / TotPop

# Store infected alive in each year till 2016
DiffinfectedAliveInYear<-c(0,rep(0,index2016-1))
for(i in 1:index2016) { 
  if(i==1) { DiffinfectedAliveInYear[i]<-infectedAliveInYear[i] } else {
    DiffinfectedAliveInYear[i]<-infectedAliveInYear[i]-infectedAliveInYear[i-1]
  }
}
DiffinfectedAliveInYear<-c(0,DiffinfectedAliveInYear[1:index2016-1])


Back$TotalInfectedInYear_LL <- Back$incidence_LL / pactive_before_die[3]  # LL Est.
infected2016_LL <- 0
index2016 <- length(Back$TotalInfectedInYear_LL)
infectedAliveInYear_LL<-c(0,rep(0,index2016-1))
pid <- pids_source[,3]
for(i in 1:index2016){
  interval <- index2016 - i
  ind <- min(interval,length(pid))
  if(interval != 0){
    infected2016_LL <- infected2016_LL + Back$TotalInfectedInYear_LL[i] * (1 - sum(pid[1:ind])) *
      (1 - sum(prob_die_at_interval[1:interval]))
    infected2016_LL <- infected2016_LL + .5 * Back$TotalInfectedInYear_LL[i] * pid[ind] * (1 - sum(prob_die_at_interval[1:interval]) + .5 * prob_die_at_interval[interval] )
  }
  infectedAliveInYear_LL[i]<-infected2016_LL
}

infected2016_LL / TotPop

# Store infected alive in each year till 2016
DiffinfectedAliveInYear_LL<-c(0,rep(0,index2016-1))
for(i in 1:index2016) { 
  if(i==1) { DiffinfectedAliveInYear_LL[i]<-infectedAliveInYear_LL[i] } else {
    DiffinfectedAliveInYear_LL[i]<-infectedAliveInYear_LL[i]-infectedAliveInYear_LL[i-1]
  }
}
DiffinfectedAliveInYear_LL<-c(0,DiffinfectedAliveInYear_LL[1:index2016-1])

Back$TotalInfectedInYear_UL <- Back$incidence_UL / pactive_before_die[2]  # UL Est.
infected2016_UL <- 0
index2016 <- length(Back$TotalInfectedInYear_UL)
infectedAliveInYear_UL<-c(0,rep(0,index2016-1))
pid <- pids_source[,2]
for(i in 1:index2016){
  interval <- index2016 - i
  ind <- min(interval,length(pid))
  if(interval != 0){
    infected2016_UL <- infected2016_UL + Back$TotalInfectedInYear_UL[i] * (1 - sum(pid[1:ind])) *
      (1 - sum(prob_die_at_interval[1:interval]))
    infected2016_UL <- infected2016_UL + .5 * Back$TotalInfectedInYear_UL[i] * pid[ind] * (1 - sum(prob_die_at_interval[1:interval]) + .5 * prob_die_at_interval[interval] )
    
    }
  infectedAliveInYear_UL[i]<-infected2016_UL
}

infected2016_UL / TotPop

# Store infected alive in each year till 2016
DiffinfectedAliveInYear_UL<-c(0,rep(0,index2016-1))
for(i in 1:index2016) { 
  if(i==1) { DiffinfectedAliveInYear_UL[i]<-infectedAliveInYear_UL[i] } else {
    DiffinfectedAliveInYear_UL[i]<-infectedAliveInYear_UL[i]-infectedAliveInYear_UL[i-1]
  }
}
DiffinfectedAliveInYear_UL<-c(0,DiffinfectedAliveInYear_UL[1:index2016-1])


# Store total new LTBI cases who are alive and not reactivated yet
TotNewAliveInYear <-
  data.frame(
    year = Back_est[, 1],
    TotNewAlive = DiffinfectedAliveInYear,
    TotNewAlive_LL = DiffinfectedAliveInYear_LL,
    TotNewAlive_UL = DiffinfectedAliveInYear_UL
  )

# Store Table 2: Estimated number of new latent TB infections in each year
Back_All <- merge(Back, TotNewAliveInYear ,by="year")

# Store New LTBI cases that will reactivate to TB & Total new LTBI cases & Total new LTBI Alive
write.csv(Back_All,
          file=paste0("BackAll_",SubPop,"_",AreaCode,".csv"))

# summerize current % of LTBI 
df <- data.frame(Code=AreaCode,Est=infected2016/TotPop, LL=infected2016_LL/TotPop, UL=infected2016_UL/TotPop)

AA<-AA+1
if(AA==1) {
 df_All<-df
 } else {
    df_All<-rbind(df_All,df)
}
  
write.csv(df,
          file=paste0("CurrentLatentInfections_",SubPop,"_",AreaCode,".csv"))

## store adjusted lag distribution for mortality (competing) risk  
LagP_Adj <- data.frame(PointAdj=lagp_adj,LLAdj=lagp_adjLL, ULAdj=lagp_adjUL,pids_source)

# calculate cumulative reactivation risks 
LagP_Adj[7]<-cumsum(LagP_Adj[1])
LagP_Adj[8]<-cumsum(LagP_Adj[2])
LagP_Adj[9]<-cumsum(LagP_Adj[3])
LagP_Adj[10]<-cumsum(LagP_Adj[4])
LagP_Adj[11]<-cumsum(LagP_Adj[5])
LagP_Adj[12]<-cumsum(LagP_Adj[6])
LagP_Adj[13] <- 1:85

write.csv(LagP_Adj,
          file=paste0("LagP_Adj_",SubPop,"_",AreaCode,".csv"))

}

#  for now, it uses TopPop 2016 for all calculations, even if you choose a differnt last year. 
# save all results in one csv file
df_All2<-merge(df_All,cases[which(cases$SubPop == SubPop & cases$Year==lasty),], by.x="Code", by.y="AreaCode", all.x = FALSE, all.y = FALSE) 
df_All2<-merge(df_All2,totalpop[which(totalpop$SubPop == SubPop & totalpop$Year==2016),c("AreaCode","SubPop","Pop","Year")], by.x="Code", by.y="AreaCode", all.x = FALSE, all.y = FALSE) 
df_All2$LTBICase<-round(df_All2$Est*df_All2$Pop)
df_All2$LTBICaseLL<-round(df_All2$LL*df_All2$Pop)
df_All2$LTBICaseUL<-round(df_All2$UL*df_All2$Pop)
df_All2<-df_All2[,c(1:5,7:9,11,13:15)]

df_All2
write.csv(df_All2,
          file=paste0("CurrentLatentInfections_",SubPop,"_",A[1],"_to_",S,".csv"))


### Making some useful graphs

require(ggplot2)
ylbl<-c("Adj. Risk % (Uncertainty Limit)","Risk % (Uncertainty Limit)","Adj. Cumulative Risk % (Uncertainty Limit)","Cumulative Risk % (Uncertainty Limit)")

# Graph New LTBI, Reported TB, and New LTBI alive (Uncertainty Limit)

Back_All2<-Back_All[which(as.numeric(paste(Back_All$year))>1931),]
LagP_Gr0 <-ggplot(Back_All2, aes(x=Back_All2$year, y=Back_All2$TotNewAlive, group=1)) +
  geom_line() +
  geom_point(shape=21, size=2, fill="black") +
  geom_errorbar(width=.4,size=.6, aes(ymin=Back_All2$TotNewAlive_LL, ymax=Back_All2$TotNewAlive_UL),colour="green") +
  aes(x=Back_All2$year, y=Back_All2$TotalInfectedInYear, group=1) +
  geom_errorbar(width=.4,size=.3, aes(ymin=Back_All2$TotalInfectedInYear_LL, ymax=Back_All2$TotalInfectedInYear_UL)) +
  geom_point(shape=21, size=2, fill="green") +
  aes(x=Back_All2$year, y=Back_All2$count, group=1) +
  geom_point(shape=21, size=3, fill="red")+
  ylab("Number of individuals with LTBI (95%UI)") +
  xlab("Year") +
  theme(axis.text.x = element_text(size=7,angle = 90, hjust = 1, vjust = 0.5)) +
  geom_point(data = Back_All2, aes(Back_All2$year, Back_All2$TotNewAlive), color="green") +
  geom_point(data = Back_All2, aes(Back_All2$year, Back_All2$TotalInfectedInYear), color="black")
LagP_Gr0  

fig3 <-ggplot(Back_All2, aes(x=Back_All2$year, y=Back_All2$TotNewAlive, group=1)) +
  geom_line() +
  geom_point(shape=21, size=2, fill="black") +
  geom_errorbar(width=.4,size=.6, aes(ymin=Back_All2$TotNewAlive_LL, ymax=Back_All2$TotNewAlive_UL),colour="green") +
  aes(x=Back_All2$year, y=Back_All2$TotalInfectedInYear, group=1) +
  geom_errorbar(width=.4,size=.3, aes(ymin=Back_All2$TotalInfectedInYear_LL, ymax=Back_All2$TotalInfectedInYear_UL)) +
  geom_point(shape=21, size=2, fill="green") +
  ylab("Number of individuals with LTBI (95%UI)") +
  xlab("Year") +
  theme(axis.text.x = element_text(size=7,angle = 90, hjust = 1, vjust = 0.5)) +
  geom_point(data = Back_All2, aes(Back_All2$year, Back_All2$TotNewAlive), color="green",shape = 15) +
  geom_point(data = Back_All2, aes(Back_All2$year, Back_All2$TotalInfectedInYear), color="black", shape = 16) 
fig3 



# Graph Lag : Adj. Risk % (Uncertainty Limit)
  DATA<-LagP_Adj[,c(1:3,13)]
  LagP_Gr1 <- ggplot(data = DATA, aes(x = DATA[, 4], y = DATA[, 1]*100), group = 1) + geom_point(size = 1) +
  geom_errorbar(aes(ymax = DATA[, 3]*100, ymin = DATA[, 2]*100)) + labs(x = "Years after latent TB infection", y = ylbl[1]) +
  geom_line(data = DATA, size = 0.5, color = "RED", aes(x = DATA[, 4], y = DATA[, 1]*100), group = 1) +
  theme_bw(base_size = 14) + scale_x_continuous(limits=c(0, lag))
LagP_Gr1

# Graph Lag : Risk % (Uncertainty Limit)
DATA<-LagP_Adj[,c(4:6,13)]
LagP_Gr2 <- ggplot(data = DATA, aes(x = DATA[, 4], y = DATA[, 1]*100), group = 1) + geom_point(size = 1) +
  geom_errorbar(aes(ymax = DATA[, 3]*100, ymin = DATA[, 2]*100)) + labs(x = "Years after latent TB infection", y = ylbl[2]) +
  geom_line(data = DATA, size = 0.5, color = "RED", aes(x = DATA[, 4], y = DATA[, 1]*100), group = 1) +
  theme_bw(base_size = 14) + scale_x_continuous(limits=c(0, lag))
LagP_Gr2

# Graph Lag : Adj. Cumulative Risk % (Uncertainty Limit)
DATA<-LagP_Adj[,c(7:9,13)]
LagP_Gr3 <- ggplot(data = DATA, aes(x = DATA[, 4], y = DATA[, 1]*100), group = 1) + geom_point(size = 1) +
  geom_errorbar(aes(ymax = DATA[, 3]*100, ymin = DATA[, 2]*100)) + labs(x = "Years after latent TB infection", y = ylbl[3]) +
  geom_line(data = DATA, size = 0.5, color = "RED", aes(x = DATA[, 4], y = DATA[, 1]*100), group = 1) +
  theme_bw(base_size = 14) + scale_x_continuous(limits=c(0, lag))
LagP_Gr3

# Graph Lag : Cumulative Risk % (Uncertainty Limit)
DATA<-LagP_Adj[,c(10:12,13)]
LagP_Gr4 <- ggplot(data = DATA, aes(x = DATA[, 4], y = DATA[, 1]*100), group = 1) + geom_point(size = 1) +
  geom_errorbar(aes(ymax = DATA[, 3]*100, ymin = DATA[, 2]*100)) + labs(x = "Years after latent TB infection", y = ylbl[4]) +
  geom_line(data = DATA, size = 0.5, color = "RED", aes(x = DATA[, 4], y = DATA[, 1]*100), group = 1) +
  theme_bw(base_size = 14) + scale_x_continuous(limits=c(0, lag))
LagP_Gr4



