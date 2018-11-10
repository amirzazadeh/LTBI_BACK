# set your own directory 

library(xlsx)

# TB Data.xlsx has 4 sheets 

# function to remove commas or spaces in 1,000 slots for numerals and converts to numeric, else will be factors
parse_num <- function(x){
  x <- as.character(x)
  x <- gsub(',', '', x)
  x <- gsub(' ', '', x)
  as.numeric(x)
}

d <- read.xlsx("TB Data.xlsx", 
                  sheetIndex = "ActiveTB")
d <- d[,-5] # drop Source
#str(d)
d$TB <- parse_num(d$TB) 

#fill out all state AreaCode to 2 digits 
levels(d$AreaCode)[levels(d$AreaCode) == "1"] <- "01"
levels(d$AreaCode)[levels(d$AreaCode) == "2"] <- "02"
levels(d$AreaCode)[levels(d$AreaCode) == "4"] <- "04"
levels(d$AreaCode)[levels(d$AreaCode) == "5"] <- "05"
levels(d$AreaCode)[levels(d$AreaCode) == "6"] <- "06"
levels(d$AreaCode)[levels(d$AreaCode) == "8"] <- "08"
levels(d$AreaCode)[levels(d$AreaCode) == "9"] <- "09"
# rename levels 991 (USB for all US) and 992 (FB for all US) to 99 
levels(d$AreaCode)[levels(d$AreaCode) %in% c("991", "992")] <- "99"
# make Name consistent with other sheets
levels(d$Name)[levels(d$Name) == "United States"] <- "US" 
cases <- d 

d <- read.xlsx("TB Data.xlsx", 
                sheetIndex = "TotPop")
# drop CA data from earlier project 
d <- d[-(1:5),]
# get rid of redundant CA county names and lower/upper case mismatch on Alameda
# likewise with AreaCode, in other words drop redundant factors from both
d$Name <- factor(d$Name)
d$AreaCode <- factor(d$AreaCode)
# make SubPop levels consistent 
levels(d$SubPop)[levels(d$SubPop) == "us-Born"] <- "US-Born" 
levels(d$SubPop)[levels(d$SubPop) == "total"] <- "Total"
d$Pop <- parse_num(d$Pop)
levels(d$AreaCode)[levels(d$AreaCode) == "1"] <- "01"
levels(d$AreaCode)[levels(d$AreaCode) == "2"] <- "02"
levels(d$AreaCode)[levels(d$AreaCode) == "4"] <- "04"
levels(d$AreaCode)[levels(d$AreaCode) == "5"] <- "05"
levels(d$AreaCode)[levels(d$AreaCode) == "6"] <- "06"
levels(d$AreaCode)[levels(d$AreaCode) == "8"] <- "08"
levels(d$AreaCode)[levels(d$AreaCode) == "9"] <- "09"
totalpop <- d

d <- read.xlsx("TB Data.xlsx", 
                sheetIndex = "death_pop")
d$Deaths <- parse_num(d$Deaths)
d$Population <- parse_num(d$Population)
levels(d$AreaCode)[levels(d$AreaCode) == "1"] <- "01"
levels(d$AreaCode)[levels(d$AreaCode) == "2"] <- "02"
levels(d$AreaCode)[levels(d$AreaCode) == "4"] <- "04"
levels(d$AreaCode)[levels(d$AreaCode) == "5"] <- "05"
levels(d$AreaCode)[levels(d$AreaCode) == "6"] <- "06"
levels(d$AreaCode)[levels(d$AreaCode) == "8"] <- "08"
levels(d$AreaCode)[levels(d$AreaCode) == "9"] <- "09"
names(d)[names(d) == "Single.Year.Ages.Code"] <- "Age"
mortality <- d 

d <- read.xlsx("TB Data.xlsx", 
                sheetIndex = "Lag")
d <- d[1:85, 1:4] 
d$Lag <- factor(d$Lag) 
d$Point <- parse_num(d$Point)
d$LL <- parse_num(d$LL)
d$UL <- parse_num(d$UL)
react <- d
rm(d)



