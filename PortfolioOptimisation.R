#======================================#
#**************************************#

# Portfolio Optimisation of Manually selected assets

#**************************************#
#======================================#

# Imports
library(GA)
library(dplyr)
library(ggplot2)
library(quantmod)
library(tbl2xts)
library(tidyquant)
library(yfR)
library(tidyr)

# Loading Assets
# The assets selected for the initial portfolio are from a diverse range of 
# industries during the financial trading period of *2021-2022*
assets <- c("SMPL", "WING", "NFLX", "LLY", "ETSY", 
            "CAT", "HAS", "SONY", "IBM", "LEVI")

# Prepare "training data" (2021-2022) to optimise for future "test data" (2022-2023).
trainReturns = lapply(assets, function(sym) {
  dailyReturn(na.omit(getSymbols(sym, from="2021-01-01", to="2022-01-01", auto.assign=FALSE)))
})

trainReturns<-do.call(merge.xts,trainReturns)
colnames(trainReturns) <- assets
nStocks <- ncol(trainReturns)

# Fitness function will be a multi-objective solution - max returns and min risk (covariance).
meanReturnsTrain <- apply(trainReturns,2,mean)
covMatTrain <- cov(trainReturns)

getReturnsTrain <- function(x) {
  (sum(x * meanReturnsTrain)*252) # multiply by 252 to annualise returns
}

getRiskTrain <- function(x) {
  (sqrt(t(x) %*%(covMatTrain*252) %*% x)) # multiply by 252 to annualise risk
}

# The Sharpe Ratio is used as a method for measuring risk adjusted returns.
# Optimising a portfolio is achieved by simply calculating the ratio of returns over risk.
sharpe <- function(x) {
  return(getReturnsTrain(x)/getRiskTrain(x))
}

# The vector of weights represents the proportion of an investment associated with each asset.
# The GA will be used to predict the optimal weights based on the Sharpe Ratio.
lower <- c(rep(0,nStocks))
upper <- c(rep(1,nStocks))

GA <- ga(type = "real-valued", fitness = sharpe, lower = lower, 
         upper = upper, maxiter = 1000, popSize = 100, seed=69)

summary(GA)

# Plot showing the maximisation of the fitness function suggests that 1000 
# generations might be overkill as it starts to stablise around 400 iterations.
plot(GA)

# The weights are required to be normalised so that they do not sum to > 1. 
# The ratio of each weight relative will remain the same so results should not be affected.
optimal_weights <-round(GA@solution / sum(GA@solution), 4)
colnames(optimal_weights) <- assets
optimal_weights

# It is clear that the LLY stock outperforms the majority of the other stocks in the training data. 
# The projected returns and risk from applying the GA's derived weights can also be calculated. 
ga_weights <- as.vector(optimal_weights)

gaRetTrain <- round(getReturnsTrain(ga_weights),4)
gaRiskTrain <- round(getRiskTrain(ga_weights),4)
gaSharpeTrain <- round(gaRetTrain/gaRiskTrain,4)

# Typically a Sharpe Ratio value greater than 1 is considered ideal so 2.1745 is considered very good in this case.
print(paste("Expected Returns from evolved GA weights:", gaRetTrain))
print(paste("Expected Risk from evolved GA weights:", as.vector(gaRiskTrain)))
print(paste("Sharpe Ratio from evolved GA weights:", gaSharpeTrain))

# Evaluation on future "test" data (returns from 2022-2023)
testReturns = lapply(assets, function(sym) {
  dailyReturn(na.omit(getSymbols(sym, from="2022-01-01", to="2023-01-01", auto.assign=FALSE)))
})

testReturns<-do.call(merge.xts,testReturns)
colnames(testReturns) <- assets

meanReturnsTest <- apply(testReturns,2,mean)
covMatTest <- cov(testReturns)

getReturnsTest <- function(x) {
  (sum(x * meanReturnsTest)*252)
}

getRiskTest <- function(x) {
  (sqrt(t(x) %*%(covMatTest*252) %*% x))
}

# Applying the Sharpe Ratio, the value is 0.0568 for the test data compared to 
# 2.175 for the train data. The initial concern of volatile markets after 
# the effect of COVID-19 is a possible explanation behind such a drastic drop.

gaRetTest <- round(getReturnsTest(ga_weights),4)
gaRiskTest <- round(getRiskTest(ga_weights),4)
gaSharpeTest <- round(gaRetTest/gaRiskTest,4)

print(paste("Expected Returns using Test Data:", gaRetTest))
print(paste("Expected Risk using Test Data:", as.vector(gaRiskTest)))
print(paste("Sharpe Ratio using Test Data:", gaSharpeTest))

#======================================#
#**************************************#

# Selection of best assets

#**************************************#
#======================================#

# Asset selection using GA
# Using the SP500 stocks would be a suitable starting point as 
# these typically are valued as high performing stocks. 
# A sample of 100 from the SP500 assets were taken but this process
# can be easily upscaled to using the entire pool of 500

sp500 <- yf_index_composition("SP500")
set.seed(1)
assets <- sample(sp500$ticker, 100)

newDailyReturns <- assets %>%
  tq_get(get  = "stock.prices",
         from = "2021-01-01",
         to   = "2022-01-01") %>%
  group_by(symbol) %>%
  tq_transmute(select     = adjusted, 
               mutate_fun = periodReturn, 
               period     = "daily", 
               col_rename = "returns")

dailyReturnsXTS <- tbl_xts(newDailyReturns, cols_to_xts = returns, spread_by = symbol)
symbols <- colnames(dailyReturnsXTS)

# Get mean returns and risk of each asset
meanReturns <- apply(dailyReturnsXTS,2,mean)
covMat <- cov(dailyReturnsXTS) 

getReturns <- function(x) {
  (sum(x * meanReturns)*252)
}

getRisk <- function(x) {
  (sqrt(t(x) %*%(covMat*252) %*% x))
}

sharpe <- function(x){
  return(getReturns(x)/getRisk(x))
}

# Objective function used for selection of asset
# based on sharpe ratio - threshold value of 1 for SR.
selection = function(x) {
  if (sharpe(x) <1) return(0)
  else return (sharpe(x)) 
}

# Custom population and mutation functions
# used to ensure a maximum of 10 assets selected.
# e.g. top 10 instead of top 50 for example 
customPopulation <- function(k){
  function(GA){
    m <- matrix(0, ncol = GA@nBits, nrow = GA@popSize)
    for(i in seq_len(GA@popSize))
      m[i, sample(GA@nBits, k)] <- 1 
    m
  }
}

customMutation <- function(GA, parent){
  index <- which(GA@population[parent,] == 1) 
  change <- sample(3, 1)
  index[sample(length(index), change)] <- sample(setdiff(seq_len(GA@nBits), index), change)
  parent <- integer(GA@nBits)
  parent[index] <- 1
  parent
}

gaControl("binary" = list(selection = "gabin_tourSelection", 
                          crossover = "gabin_uCrossover"))

# GA uses binary representation:
# 1 signals that an asset should be selected, 0 signals it should not.
selectionGA= ga(type='binary', fitness=function(x){selection(x)}, 
                nBits=ncol(dailyReturnsXTS), names=symbols, maxiter = 500,
                mutation=customMutation, population = customPopulation(10),
                popSize=100, keepBest=TRUE, seed=1)

best_assets = symbols[selectionGA@solution[1,]==1]
best_assets

#======================================#
#**************************************#

# Portfolio Optimisation of automatically selected assets

#**************************************#
#======================================#

newReturns = lapply(best_assets, function(sym) {
  dailyReturn(na.omit(getSymbols(sym, from="2021-01-01", to="2022-01-01", auto.assign=FALSE)))
})
newReturns<-do.call(merge.xts, newReturns)

colnames(newReturns) <- best_assets

meanReturnsNew <- apply(newReturns,2,mean)
covMatNew <- cov(newReturns)

getReturnsNew <- function(x) {
  (sum(x * meanReturnsNew)*252)
}

getRiskNew <- function(x) {
  (sqrt(t(x) %*%(covMatNew*252) %*% x))
}

sharpeNew <- function(x) {
  return(getReturnsNew(x)/getRiskNew(x))
}

newGA <- ga(type = "real-valued", fitness = sharpeNew, lower = lower, 
            upper = upper, maxiter = 1000, popSize = 100, seed=69)

new_optimal_weights <-round(newGA@solution / sum(newGA@solution), 4)
colnames(new_optimal_weights) <- best_assets
new_optimal_weights

testReturnsNew = lapply(best_assets, function(sym) {
  dailyReturn(na.omit(getSymbols(sym, from="2022-01-01", to="2023-01-01", auto.assign=FALSE)))
})

testReturnsNew<-do.call(merge.xts,testReturnsNew)
colnames(testReturnsNew) <- best_assets

meanReturnsNewTest <- apply(testReturnsNew,2,mean)
covMatNewTest <- cov(testReturnsNew)

getReturnsNewTest <- function(x) {
  (sum(x * meanReturnsNewTest)*252)
}

getRiskNewTest <- function(x) {
  (sqrt(t(x) %*%(covMatNewTest*252) %*% x))
}

sharpeNewTest <- function(x) {
  return(getReturnsNewTest(x)/getRiskNewTest(x))
}

new_ga_weights <- as.vector(new_optimal_weights)

gaSharpeNew <- round(sharpeNew(new_ga_weights),4)
gaSharpeNewT <- round(sharpeNewTest(new_ga_weights),4)

# Comparing GA selected assets vs. manually selected assets
print(paste("Sharpe using automatically selected assets (Train Data):", gaSharpeNew))
print(paste("Sharpe using automatically selected assets (Test Data)", gaSharpeNewT))

print(paste("Sharpe using manually selected assets (Train Data):", gaSharpeTrain))
print(paste("Sharpe using manually selected assets (Test Data):", gaSharpeTest))

#======================================#
#**************************************#

# Additional: Compare with balanced and random portfolios

#**************************************#
#======================================#

# Comparison with balanced portfolio (assets all have equal weightings: 0.1 in this case)
balanced_weights <- rep(0.1, length(assets))

bRetTrain <- round(getReturnsTrain(balanced_weights),4)
bRiskTrain <- round(getRiskTrain(balanced_weights),4)
bSharpeTrain <- round(bRetTrain/bRiskTrain,4)

bRetTest <- round(getReturnsTest(balanced_weights),4)
bRiskTest <- round(getRiskTest(balanced_weights),4)
bSharpeTest <- round(bRetTest/bRiskTest,4)

print(paste("Balanced Portfolio Returns (Train):", bRetTrain))
print(paste("Balanced Portfolio Risk (Train):", bRiskTrain))
print(paste("Balanced Portfolio Sharpe (Train):", bSharpeTrain))

print(paste("Balanced Portfolio Returns (Test):", bRetTest))
print(paste("Balanced Portfolio Risk (Test):", bRiskTest))
print(paste("Balanced Portfolio Sharpe (Test):", bSharpeTest))

generateRndWeights <- function() {
  r <- runif(length(assets))
  randomWeights <- sapply(r, function(x) x / sum(r))
}

# Comparison using weights obtained through RSA - random search algo (all must sum to 1)
# Either maximise returns or minimise risks
r <- list()
for(i in (1:1000)){
  r[[i]] <- runif(nStocks,0,1)
}

Rscaled <- lapply(r, function(r) r / sum(r))
Rreturns <- sapply(Rscaled, getReturnsTrain)
Rrisks <- sapply(Rscaled[], getRiskTrain)

bestRSAMax <- unlist(Rscaled[which.max(Rreturns)])
bestRSAMin <- unlist(Rscaled[which.min(Rrisks)])

# RSA maximising returns
# Train data
rsaRetTrainMax <- round(getReturnsTrain(unlist(bestRSAMax)),4)
rsaRiskTrainMax <- round(getRiskTrain(unlist(bestRSAMax)),4)
rsaTrainMaxSR <- round(rsaRetTrainMax/rsaRiskTrainMax,4)

print(paste("RSA(max returns) Returns:", rsaRetTrainMax))
print(paste("RSA(max returns) Risk:", rsaRiskTrainMax))
print(paste("RSA(max returns) Sharpe:", rsaTrainMaxSR))

# Test data
rsaRetTestMax <- round(getReturnsTest(unlist(bestRSAMax)),4)
rsaRiskTestMax <- round(getRiskTest(unlist(bestRSAMax)),4)
rsaTestMaxSR <- round(rsaRetTestMax/rsaRiskTestMax,4)

print(paste("RSA(max returns) Returns:", rsaRetTestMax))
print(paste("RSA(max returns) Risk:", rsaRiskTestMax))
print(paste("RSA(max returns) Sharpe:", rsaTestMaxSR))

# RSA minimising risk
# Train data
rsaRetTrainMin <- round(getReturnsTrain(unlist(bestRSAMin)),4)
rsaRiskTrainMin <- round(getRiskTrain(unlist(bestRSAMin)),4)
rsaTrainMinSR <- round(rsaRetTrainMin/rsaRiskTrainMin,4)

print(paste("RSA(min risk) Returns:", rsaRetTrainMin))
print(paste("RSA(min risk) Risk:", rsaRiskTrainMin))
print(paste("RSA(min risk) Sharpe:", rsaTrainMinSR))

# Test data
rsaRetTestMin <- round(getReturnsTest(unlist(bestRSAMin)),4)
rsaRiskTestMin <- round(getRiskTest(unlist(bestRSAMin)),4)
rsaTestMinSR <- round(rsaRetTestMin/rsaRiskTestMin,4)

print(paste("RSA(min risk) Returns:", rsaRetTestMin))
print(paste("RSA(min risk) Risk:", rsaRiskTestMin))
print(paste("RSA(min risk) Sharpe:", rsaTestMinSR))


