 
#Install Relevant Packages

library(tidyverse)
library(dplyr)
library(timeSeries)
library("data.table")
if(require("quantmod",quietly=TRUE))
  install.packages("quantmod",dependencies=TRUE,repos=c(CRAN="http://cran.rstudio.com"))
library("quantmod")
library("xts")
install.packages("fPortfolio")
library(fPortfolio)

library(ggplot2)
library(ggplotify)

# Small Cap

getSymbols("AOS")
barChart(AOS)

getSymbols("DX")
barChart(DX)

getSymbols("NNBR")
barChart(NNBR)

getSymbols("BLDR")
barChart(BLDR)

getSymbols("GME")
barChart(GME)

getSymbols("SAVA")
barChart(SAVA)

getSymbols("PACB")
barChart(PACB)

getSymbols("AZPN")
barChart(AZPN)

getSymbols("FDS")
barChart(FDS)

getSymbols("PVH")
barChart(PVH)

getSymbols("VMI")
barChart(VMI)

# Small Cap Returns
# AOS 
areturnsAOS<-annualReturn(AOS)
areturnsAOS <-as.data.frame(areturnsAOS)
areturnsAOS <- setDT(areturnsAOS, keep.rownames = TRUE)[]  
areturnsAOS$rn = as.Date(areturnsAOS$rn)
class(areturnsAOS$rn)

# DX 
areturnsDX<-annualReturn(DX)
areturnsDX <-as.data.frame(areturnsDX)
areturnsDX <- setDT(areturnsDX, keep.rownames = TRUE)[]  
areturnsDX$rn = as.Date(areturnsDX$rn)
class(areturnsDX$rn)

# NNBR
areturnsNNBR<-annualReturn(NNBR)
areturnsNNBR <-as.data.frame(areturnsNNBR)
areturnsNNBR <- setDT(areturnsNNBR, keep.rownames = TRUE)[]  
areturnsNNBR$rn = as.Date(areturnsNNBR$rn)
class(areturnsNNBR$rn)

# BLDR
areturnsBLDR<-annualReturn(BLDR)
areturnsBLDR <-as.data.frame(areturnsBLDR)
areturnsBLDR <- setDT(areturnsBLDR, keep.rownames = TRUE)[]  
areturnsBLDR$rn = as.Date(areturnsBLDR$rn)
class(areturnsBLDR$rn)

#GME
areturnsGME<-annualReturn(GME)
areturnsGME <-as.data.frame(areturnsGME)
areturnsGME <- setDT(areturnsGME, keep.rownames = TRUE)[]  
areturnsGME$rn = as.Date(areturnsGME$rn)
class(areturnsGME$rn)

#SAVA
areturnsSAVA<-annualReturn(SAVA)
areturnsSAVA <-as.data.frame(areturnsSAVA)
areturnsSAVA <- setDT(areturnsSAVA, keep.rownames = TRUE)[]  
areturnsSAVA$rn = as.Date(areturnsSAVA$rn)
class(areturnsSAVA$rn)

#PACB
areturnsPACB<-annualReturn(PACB)
areturnsPACB <-as.data.frame(areturnsPACB)
areturnsPACB <- setDT(areturnsPACB, keep.rownames = TRUE)[]  
areturnsPACB$rn = as.Date(areturnsPACB$rn)
class(areturnsPACB$rn)

##AZPN
areturnsAZPN<-annualReturn(AZPN)
areturnsAZPN <-as.data.frame(areturnsAZPN)
areturnsAZPN <- setDT(areturnsAZPN, keep.rownames = TRUE)[]  
areturnsAZPN$rn = as.Date(areturnsAZPN$rn)
class(areturnsAZPN$rn)

##FDS
areturnsFDS<-annualReturn(FDS)
areturnsFDS <-as.data.frame(areturnsFDS)
areturnsFDS <- setDT(areturnsFDS, keep.rownames = TRUE)[]  
areturnsFDS$rn = as.Date(areturnsFDS$rn)
class(areturnsFDS$rn)

#PVH
areturnsPVH<-annualReturn(PVH)
areturnsPVH <-as.data.frame(areturnsPVH)
areturnsPVH <- setDT(areturnsPVH, keep.rownames = TRUE)[]  
areturnsPVH$rn = as.Date(areturnsPVH$rn)
class(areturnsPVH$rn)

#VMI
areturnsVMI<-annualReturn(VMI)
areturnsVMI <-as.data.frame(areturnsVMI)
areturnsVMI <- setDT(areturnsVMI, keep.rownames = TRUE)[]  
areturnsVMI$rn = as.Date(areturnsVMI$rn)
class(areturnsVMI$rn)

scap<-cbind(areturnsVMI$rn, areturnsAOS$yearly.returns, areturnsDX$yearly.returns, areturnsNNBR$yearly.returns, areturnsBLDR$yearly.returns, areturnsGME$yearly.returns, areturnsSAVA$yearly.returns, areturnsAZPN$yearly.returns,areturnsFDS$yearly.returns,areturnsPVH$yearly.returns, areturnsVMI$yearly.returns)
scap<-as.data.frame(scap)
scap$V1 <- as.Date(scap$V1)
scap$V1 <- as.Date(scap$V1, format = "%Y-%m-%d")
class(scap$V2)

scap$V2 <- as.numeric(scap$V2)
scap$V3 <- as.numeric(scap$V3)
scap$V4 <- as.numeric(scap$V4)
scap$V5 <- as.numeric(scap$V5)
scap$V6 <- as.numeric(scap$V6)
scap$V7 <- as.numeric(scap$V7)
scap$V8 <- as.numeric(scap$V8)
scap$V9 <- as.numeric(scap$V9)
scap$V10 <- as.numeric(scap$V10)
scap$V11 <- as.numeric(scap$V11)
scap$V11 <- as.numeric(scap$V11)

colnames(scap) <- c("Date", "AOS", "DX", "NNBR", "BLDR", "GME", "SAVA", "AZPN", "FDS", "PVH", "VMI")


rownames(scap) <- c()

scap<-timeSeries(scap[,-1], charvec = scap$Date)
class(scap)
write.csv(scap, "C:/Data Science and Analytics/DSA 5303/Final Class Project/Small Cap/scap.csv")

# long portfolio with target return at 0.15

spec1<-portfolioSpec()
setRiskFreeRate(spec1)<-.05
setTargetReturn(spec1) <- 0.15

scaplong<-efficientPortfolio(scap, spec = spec1, constraints = "LongOnly")
 
scaplongfrontier<-portfolioFrontier(scap, spec1)

tailoredFrontierPlot(object = scaplongfrontier, mText = "MV Portfolio - LongOnly Constraints",
                     risk = "Cov")

par(mfrow = c(1, 2))
col <- rampPalette(ncol(scap), "purple2green")
names <- colnames(scap)
weights <- 100 * as.vector(getWeights(scaplong))
weightedReturns <- weights * getMean(scaplong)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")

# short portfolio with target return at 0.15

spec2<-portfolioSpec()
setRiskFreeRate(spec2)<-.05
setTargetReturn(spec2) <- 0.15

setSolver(spec2)<-"solveRshortExact"

scapshort<-efficientPortfolio(scap, spec = spec2, constraints = "Short")

scapshortfrontier<-portfolioFrontier(scap, spec2)

tailoredFrontierPlot(object = scapshortfrontier, mText = "MV Portfolio - Short Constraints",
                     risk = "Cov")

par(mfrow = c(1, 2))
col <- rampPalette(ncol(scap), "purple2green")
names <- colnames(scap)
weights <- 100 * as.vector(getWeights(scapshort))
weightedReturns <- weights * getMean(scapshort)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")


# short portfolio with target risk at 0.20

spec3<-portfolioSpec()
setRiskFreeRate(spec3)<-.05
setTargetRisk(spec3) <- 0.20

setOptimize(spec3) <- 'maxReturn'
setSolver(spec3)<-"solveRshortExact"

scapshortmaxReturn<-maxreturnPortfolio(scap, spec = spec3, constraints = "Short")
scapshortfrontiermax<-portfolioFrontier(scap, spec3)

tailoredFrontierPlot(object = scapshortfrontiermax, mText = "MV Portfolio - LongOnly Constraints",
                     risk = "Cov")

par(mfrow = c(1, 2))
col <- rampPalette(ncol(scap), "purple2green")
names <- colnames(scap)
weights <- 100 * as.vector(getWeights(scapshortmaxReturn))
weightedReturns <- weights * getMean(scapshortmaxReturn)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")

#Mid Cap Stocks

getSymbols("GATX")
barChart(GATX)

getSymbols("LEN")
barChart(LEN)

getSymbols("JBLU")
barChart(JBLU)

getSymbols("MLM")
barChart(MLM)

getSymbols("NFLX")
barChart(NFLX)

getSymbols("RS")
barChart(RS)

getSymbols("RJF")
barChart(RJF)

getSymbols("SLAB")
barChart(SLAB)

getSymbols("ROST")
barChart(ROST)

getSymbols("WRB")
barChart(WRB)

# Mid Cap Returns
# GATX 
areturnsGATX<-annualReturn(GATX)
areturnsGATX <-as.data.frame(areturnsGATX)
areturnsGATX <- setDT(areturnsGATX, keep.rownames = TRUE)[]  
areturnsGATX$rn = as.Date(areturnsGATX$rn)
class(areturnsGATX$rn)

# LEN 
areturnsLEN<-annualReturn(LEN)
areturnsLEN <-as.data.frame(areturnsLEN)
areturnsLEN <- setDT(areturnsLEN, keep.rownames = TRUE)[]  
areturnsLEN$rn = as.Date(areturnsLEN$rn)
class(areturnsLEN$rn)

# JBLU 
areturnsJBLU<-annualReturn(JBLU)
areturnsJBLU <-as.data.frame(areturnsJBLU)
areturnsJBLU <- setDT(areturnsJBLU, keep.rownames = TRUE)[]  
areturnsJBLU$rn = as.Date(areturnsJBLU$rn)
class(areturnsJBLU$rn)

# MLM
areturnsMLM<-annualReturn(MLM)
areturnsMLM <-as.data.frame(areturnsMLM)
areturnsMLM <- setDT(areturnsMLM, keep.rownames = TRUE)[]  
areturnsMLM$rn = as.Date(areturnsMLM$rn)
class(areturnsMLM$rn)

# NFLX 
areturnsNFLX<-annualReturn(NFLX)
areturnsNFLX <-as.data.frame(areturnsNFLX)
areturnsNFLX <- setDT(areturnsNFLX, keep.rownames = TRUE)[]  
areturnsNFLX$rn = as.Date(areturnsNFLX$rn)
class(areturnsNFLX$rn)


# RS
areturnsRS<-annualReturn(RS)
areturnsRS <-as.data.frame(areturnsRS)
areturnsRS <- setDT(areturnsRS, keep.rownames = TRUE)[]  
areturnsRS$rn = as.Date(areturnsRS$rn)
class(areturnsRS$rn)

# RJF 
areturnsRJF<-annualReturn(RJF)
areturnsRJF <-as.data.frame(areturnsRJF)
areturnsRJF <- setDT(areturnsRJF, keep.rownames = TRUE)[]  
areturnsRJF$rn = as.Date(areturnsRJF$rn)
class(areturnsRJF$rn)


# SLAB
areturnsSLAB<-annualReturn(SLAB)
areturnsSLAB <-as.data.frame(areturnsSLAB)
areturnsSLAB <- setDT(areturnsSLAB, keep.rownames = TRUE)[]  
areturnsSLAB$rn = as.Date(areturnsSLAB$rn)
class(areturnsSLAB$rn)

# ROST 
areturnsROST<-annualReturn(ROST)
areturnsROST <-as.data.frame(areturnsROST)
areturnsROST <- setDT(areturnsROST, keep.rownames = TRUE)[]  
areturnsROST$rn = as.Date(areturnsROST$rn)
class(areturnsROST$rn)


# WRB
areturnsWRB<-annualReturn(WRB)
areturnsWRB <-as.data.frame(areturnsWRB)
areturnsWRB <- setDT(areturnsWRB, keep.rownames = TRUE)[]  
areturnsWRB$rn = as.Date(areturnsWRB$rn)
class(areturnsWRB$rn)

mcap<-cbind(areturnsWRB$rn, areturnsGATX$yearly.returns, areturnsLEN$yearly.returns, areturnsJBLU$yearly.returns, areturnsMLM$yearly.returns, areturnsNFLX$yearly.returns, areturnsRS$yearly.returns, areturnsRJF$yearly.returns,areturnsSLAB$yearly.returns,areturnsROST$yearly.returns, areturnsWRB$yearly.returns)
mcap<-as.data.frame(mcap)
mcap$V1 <- as.Date(mcap$V1)
mcap$V1 <- as.Date(mcap$V1, format = "%Y-%m-%d")
 
mcap$V2 <- as.numeric(mcap$V2)
mcap$V3 <- as.numeric(mcap$V3)
mcap$V4 <- as.numeric(mcap$V4)
mcap$V5 <- as.numeric(mcap$V5)
mcap$V6 <- as.numeric(mcap$V6)
mcap$V7 <- as.numeric(mcap$V7)
mcap$V8 <- as.numeric(mcap$V8)
mcap$V9 <- as.numeric(mcap$V9)
mcap$V10 <- as.numeric(mcap$V10)
mcap$V11 <- as.numeric(mcap$V11)
 
colnames(mcap) <- c("Date", "GATX", "LEN", "JBLU", "MLM", "NFLX", "RS", "RJF", "SLAB", "ROST", "WRB")

rownames(mcap) <- c()

mcap<-timeSeries(mcap[,-1], charvec = mcap$Date)
class(mcap)

write.csv(mcap, "C:/Data Science and Analytics/DSA 5303/Final Class Project/Mid Cap/mcap.csv")


# long portfolio with target return at 0.15

spec4<-portfolioSpec()
setRiskFreeRate(spec4)<-0.05
setTargetReturn(spec4)<- 0.15
 
mcaplong<-efficientPortfolio(mcap, spec = spec4, constraints = "LongOnly")

mcaplongfrontier<-portfolioFrontier(mcap, spec4)

tailoredFrontierPlot(object = mcaplongfrontier, mText = "MV Portfolio - LongOnly Constraints",
                     risk = "Cov")


par(mfrow = c(1, 2))
col <- rampPalette(ncol(mcap), "purple2green")
names <- colnames(mcap)
weights <- 100 * as.vector(getWeights(mcaplong))
weightedReturns <- weights * getMean(mcaplong)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")

# short portfolio with target return at 0.15

spec5<-portfolioSpec()
setRiskFreeRate(spec5)<-0.05
setTargetReturn(spec5)<- 0.15
setSolver(spec5)<-"solveRshortExact"

mcapshort<-efficientPortfolio(mcap, spec = spec5, constraints = "Short")

mcapshortfrontier<-portfolioFrontier(mcap, spec5)

tailoredFrontierPlot(object = mcapshortfrontier, mText = "MV Portfolio - Short Constraints",
                     risk = "Cov")


par(mfrow = c(1, 2))
col <- rampPalette(ncol(mcap), "purple2green")
names <- colnames(mcap)
weights <- 100 * as.vector(getWeights(mcapshort))
weightedReturns <- weights * getMean(mcapshort)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")

# short portfolio with target risk at 0.20

spec6<-portfolioSpec()
setRiskFreeRate(spec6)<-0.05
setTargetRisk(spec6)<- 0.20
setSolver(spec6)<-"solveRshortExact"

setOptimize(spec6) <- 'maxReturn'
 
mcapshortmaxReturn<-maxreturnPortfolio(scap, spec = spec6, constraints = "Short")

mcapshortfrontiermax<-portfolioFrontier(mcap, spec6)

tailoredFrontierPlot(object = mcapshortfrontiermax, mText = "MV Portfolio - Short Constraints",
                     risk = "Cov")


par(mfrow = c(1, 2))
col <- rampPalette(ncol(mcap), "purple2green")
names <- colnames(mcap)
weights <- 100 * as.vector(getWeights(mcapshortmaxReturn))
weightedReturns <- weights * getMean(mcapshortmaxReturn)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")

#### large cap

getSymbols("AAPL")
barChart(AAPL)

getSymbols("MSFT")
barChart(MSFT)

getSymbols("GOOG")
barChart(GOOG)

getSymbols("IBM")
barChart(IBM)

getSymbols("BP")
barChart(BP)

getSymbols("XOM")
barChart(XOM)

getSymbols("CSCO")
barChart(CSCO)

getSymbols("DUK")
barChart(DUK)

getSymbols("JNJ")
barChart(JNJ)

getSymbols("INTC")
barChart(INTC)


## large cap returns

# AAPL 
areturnsAAPL<-annualReturn(AAPL)
areturnsAAPL <-as.data.frame(areturnsAAPL)
areturnsAAPL <- setDT(areturnsAAPL, keep.rownames = TRUE)[]  
areturnsAAPL$rn = as.Date(areturnsAAPL$rn)
class(areturnsAAPL$rn)

# MSFT 
areturnsMSFT<-annualReturn(MSFT)
areturnsMSFT <-as.data.frame(areturnsMSFT)
areturnsMSFT <- setDT(areturnsMSFT, keep.rownames = TRUE)[]  
areturnsMSFT$rn = as.Date(areturnsMSFT$rn)
class(areturnsMSFT$rn)

# IBM 
areturnsIBM<-annualReturn(IBM)
areturnsIBM <-as.data.frame(areturnsIBM)
areturnsIBM <- setDT(areturnsIBM, keep.rownames = TRUE)[]  
areturnsIBM$rn = as.Date(areturnsIBM$rn)
class(areturnsIBM$rn)

# GOOG 
areturnsGOOG<-annualReturn(GOOG)
areturnsGOOG <-as.data.frame(areturnsGOOG)
areturnsGOOG <- setDT(areturnsGOOG, keep.rownames = TRUE)[]  
areturnsGOOG$rn = as.Date(areturnsGOOG$rn)
class(areturnsGOOG$rn)

# BP 
areturnsBP<-annualReturn(BP)
areturnsBP <-as.data.frame(areturnsBP)
areturnsBP <- setDT(areturnsBP, keep.rownames = TRUE)[]  
areturnsBP$rn = as.Date(areturnsBP$rn)
class(areturnsBP$rn)

# XOM 
areturnsXOM<-annualReturn(XOM)
areturnsXOM <-as.data.frame(areturnsXOM)
areturnsXOM <- setDT(areturnsXOM, keep.rownames = TRUE)[]  
areturnsXOM$rn = as.Date(areturnsXOM$rn)
class(areturnsXOM$rn)

# CSCO 
areturnsCSCO<-annualReturn(CSCO)
areturnsCSCO <-as.data.frame(areturnsCSCO)
areturnsCSCO <- setDT(areturnsCSCO, keep.rownames = TRUE)[]  
areturnsCSCO$rn = as.Date(areturnsCSCO$rn)
class(areturnsCSCO$rn)

# DUK 
areturnsDUK<-annualReturn(DUK)
areturnsDUK <-as.data.frame(areturnsDUK)
areturnsDUK <- setDT(areturnsDUK, keep.rownames = TRUE)[]  
areturnsDUK$rn = as.Date(areturnsDUK$rn)
class(areturnsDUK$rn)

# JNJ 
areturnsJNJ<-annualReturn(JNJ)
areturnsJNJ <-as.data.frame(areturnsJNJ)
areturnsJNJ <- setDT(areturnsJNJ, keep.rownames = TRUE)[]  
areturnsJNJ$rn = as.Date(areturnsJNJ$rn)
class(areturnsJNJ$rn)

# INTC
areturnsINTC<-annualReturn(INTC)
areturnsINTC <-as.data.frame(areturnsINTC)
areturnsINTC <- setDT(areturnsINTC, keep.rownames = TRUE)[]  
areturnsINTC$rn = as.Date(areturnsINTC$rn)
class(areturnsINTC$rn)


lcap<-cbind(areturnsINTC$rn, areturnsAAPL$yearly.returns, areturnsMSFT$yearly.returns, areturnsIBM$yearly.returns, areturnsGOOG$yearly.returns, areturnsBP$yearly.returns, areturnsXOM$yearly.returns, areturnsCSCO$yearly.returns,areturnsDUK$yearly.returns,areturnsJNJ$yearly.returns, areturnsINTC$yearly.returns)
lcap<-as.data.frame(lcap)
lcap$V1 <- as.Date(lcap$V1)
lcap$V1 <- as.Date(lcap$V1, format = "%Y-%m-%d")
class(lcap$V2)

lcap$V2 <- as.numeric(lcap$V2)
lcap$V3 <- as.numeric(lcap$V3)
lcap$V4 <- as.numeric(lcap$V4)
lcap$V5 <- as.numeric(lcap$V5)
lcap$V6 <- as.numeric(lcap$V6)
lcap$V7 <- as.numeric(lcap$V7)
lcap$V8 <- as.numeric(lcap$V8)
lcap$V9 <- as.numeric(lcap$V9)
lcap$V10 <- as.numeric(lcap$V10)
lcap$V11 <- as.numeric(lcap$V11)
 
colnames(lcap) <- c("Date", "AAPL", "MSFT", "IBM", "GOOG", "BP", "XOM", "CSCO", "DUK", "JNJ", "INTC")


rownames(lcap) <- c()

lcap<-timeSeries(lcap[,-1], charvec = lcap$Date)
class(lcap)

write.csv(lcap, "C:/Data Science and Analytics/DSA 5303/Final Class Project/Large Cap/lcap.csv")



# long position with target return of 0.15

spec7<-portfolioSpec()
setRiskFreeRate(spec7)<-0.05
setTargetReturn(spec7)<- 0.15
 
lcaplong<-efficientPortfolio(lcap, spec = spec7, constraints = "LongOnly")
 
lcaplongfrontier<-portfolioFrontier(lcap, spec7)

tailoredFrontierPlot(object = lcaplongfrontier, mText = "MV Portfolio - LongOnly Constraints",
                     risk = "Cov")


par(mfrow = c(1, 2))
col <- rampPalette(ncol(lcap), "purple2green")
names <- colnames(lcap)
weights <- 100 * as.vector(getWeights(lcaplong))
weightedReturns <- weights * getMean(mcaplong)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")

# short position with target return of 0.15

spec8<-portfolioSpec()
setRiskFreeRate(spec8)<-0.05
setTargetReturn(spec8)<- 0.15
setSolver(spec8)<-"solveRshortExact"

lcapshort<-efficientPortfolio(lcap, spec = spec8, constraints = "Short")

lcapshortfrontier<-portfolioFrontier(lcap, spec8)

tailoredFrontierPlot(object = lcapshortfrontier, mText = "MV Portfolio - ShortOnly Constraints",
                     risk = "Cov")


par(mfrow = c(1, 2))
col <- rampPalette(ncol(lcap), "purple2green")
names <- colnames(lcap)
weights <- 100 * as.vector(getWeights(lcapshort))
weightedReturns <- weights * getMean(lcapshort)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")


# short portfolio with target risk at 0.20

spec9<-portfolioSpec()
setRiskFreeRate(spec9)<-0.05
setTargetRisk(spec9)<- 0.20
setSolver(spec9)<-"solveRshortExact"

setOptimize(spec9) <- 'maxReturn'

lcapshortmaxReturn<-maxreturnPortfolio(lcap, spec = spec9, constraints = "Short")

lcapshortfrontiermax<-portfolioFrontier(lcap, spec9)

tailoredFrontierPlot(object = lcapshortfrontiermax, mText = "MV Portfolio - Short Constraints",
                     risk = "Cov")


par(mfrow = c(2, 2))
col <- rampPalette(ncol(lcap), "purple2green")
names <- colnames(lcap)
weights <- 100 * as.vector(getWeights(lcapshortmaxReturn))
weightedReturns <- weights * getMean(lcapshortmaxReturn)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")

# combination of small, mid cap, and large cap stocks

mixedcap<-cbind(areturnsAAPL$rn, areturnsAAPL$yearly.returns, areturnsMSFT$yearly.returns, areturnsBP$yearly.returns, areturnsXOM$yearly.returns,  areturnsDUK$yearly.returns, 
                areturnsJNJ$yearly.returns, areturnsNFLX$yearly.returns, areturnsRS$yearly.returns, areturnsRJF$yearly.returns,areturnsSLAB$yearly.returns,
                areturnsWRB$yearly.returns, areturnsAOS$yearly.returns, areturnsFDS$yearly.returns, areturnsAZPN$yearly.returns)

mixedcap<-as.data.frame(mixedcap)
mixedcap$V1 <- as.Date(mixedcap$V1)
mixedcap$V1 <- as.Date(mixedcap$V1, format = "%Y-%m-%d")
class(mixedcap$VMI)

mixedcap$V2 <- as.numeric(mixedcap$V2)
mixedcap$V3 <- as.numeric(mixedcap$V3)
mixedcap$V4 <- as.numeric(mixedcap$V4)
mixedcap$V5 <- as.numeric(mixedcap$V5)
mixedcap$V6 <- as.numeric(mixedcap$V6)
mixedcap$V7 <- as.numeric(mixedcap$V7)
mixedcap$V8 <- as.numeric(mixedcap$V8)
mixedcap$V9 <- as.numeric(mixedcap$V9)
mixedcap$V10 <- as.numeric(mixedcap$V10)
mixedcap$V11 <- as.numeric(mixedcap$V11)
mixedcap$V12 <- as.numeric(mixedcap$V12)
mixedcap$V13 <- as.numeric(mixedcap$V13)
mixedcap$V14 <- as.numeric(mixedcap$V14)
mixedcap$V15 <- as.numeric(mixedcap$V15)
 

colnames(mixedcap) <- c("Date", "AAPL", "MSFT", "BP", "XOM", "DUK", "JNJ", "NFLX", "RS", "RJF", "SLAB", "WRB", "AOS", "FDS", "AZPN")


rownames(mixedcap) <- c()

mixedcap<-timeSeries(mixedcap[,-1], charvec = mixedcap$Date)
class(mixedcap)

write.csv(mixedcap, "C:/Data Science and Analytics/DSA 5303/Final Class Project/Mixed Cap/mixedcap.csv")

# long position with target return of 0.15

spec10<-portfolioSpec()
setRiskFreeRate(spec10)<-0.05
setTargetReturn(spec10)<- 0.15

mixedcaplong<-efficientPortfolio(mixedcap, spec = spec10, constraints = "LongOnly")

mixedcaplongfrontier<-portfolioFrontier(mixedcap, spec10)

tailoredFrontierPlot(object = mixedcaplongfrontier, mText = "MV Portfolio - LongOnly Constraints",
                     risk = "Cov")

par(mfrow = c(1, 1))

par(mfrow = c(1, 2))
col <- rampPalette(ncol(mixedcap), "purple2green")
names <- colnames(mixedcap)
weights <- 100 * as.vector(getWeights(mixedcaplong))
weightedReturns <- weights * getMean(mixedcaplong)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")

# short position with target return of 0.15

spec11<-portfolioSpec()
setRiskFreeRate(spec11)<-0.05
setTargetReturn(spec11)<- 0.15
setSolver(spec11)<-"solveRshortExact"

mixedcapshort<-efficientPortfolio(mixedcap, spec = spec11, constraints = "Short")

mixedcapshortfrontier<-portfolioFrontier(mixedcap, spec11)

tailoredFrontierPlot(object = mixedcapshortfrontier, mText = "MV Portfolio - ShortOnly Constraints",
                     risk = "Cov")


par(mfrow = c(1, 2))
col <- rampPalette(ncol(mixedcap), "purple2green")
names <- colnames(mixedcap)
weights <- 100 * as.vector(getWeights(mixedcapshort))
weightedReturns <- weights * getMean(mixedcapshort)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")


# short portfolio with target risk at 0.20

spec12<-portfolioSpec()
setRiskFreeRate(spec12)<-0.05
setTargetRisk(spec12)<- 0.20
setSolver(spec12)<-"solveRshortExact"

setOptimize(spec12) <- 'maxReturn'

mixedcapshortmaxReturn<-maxreturnPortfolio(mixedcap, spec = spec12, constraints = "Short")

mixedcapshortfrontiermax<-portfolioFrontier(mixedcap, spec12)

tailoredFrontierPlot(object = mixedcapshortfrontiermax, mText = "MV Portfolio - Short Constraints",
                     risk = "Cov")


par(mfrow = c(1, 2))
col <- rampPalette(ncol(mixedcap), "purple2green")
names <- colnames(mixedcap)
weights <- 100 * as.vector(getWeights(mixedcapshortmaxReturn))
weightedReturns <- weights * getMean(mixedcapshortmaxReturn)
barplot(height = weights, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Portfolio Weights", xlab = "Weights %")
barplot(height = weightedReturns, names.arg = names, horiz = TRUE, las = 1, col = col)
title(main = "Weighted Portfolio Returns", xlab = "Weighted Returns %")

# long frontier for all frontier points for mixed cap

weightsPlot(mixedcaplongfrontier)
text <- "Mean-Variance Portfolio - Long Only Constraints"
mtext(text, side = 3, line = 3, font = 2, cex = 0.9)
col <- rampPalette(ncol(mixedcap))
weightedReturnsPlot(mixedcaplongfrontier)



# short frontier for all frontier points for mixed cap

weightsPlot(mixedcapshortfrontier)
text <- "Mean-Variance Portfolio - Short Sell"
mtext(text, side = 3, line = 3, font = 2, cex = 0.9)
weightedReturnsPlot(mixedcapshortfrontier, col=col)
covRiskBudgetsPlot(mixedcapshortfrontier)


# shiny basic plots app

library(shiny)

shinyApp(
ui <- fluidPage(    
  radioButtons("radio", 
               label = HTML('<FONT color="red"><FONT size="5pt">Welcome to Portfolio Selection Program</FONT></FONT><br> <b>Which Stock do you want?</b>'),
               c("scaplong" = 1, "scapshort" = 2, "mcaplong"=3, "mcapshort" = 4, "lcaplong" = 5, "lcapshort" = 6)),      
  
  mainPanel(
    plotOutput("plot")
  )
),


server = function(input, output) {
  
  
  output$plot <- renderPlot({
    
    # Now let us print some stock data
    if(input$radio == 1){
   
      plot(scaplongfrontier, c(1,2,3,4,8))
    }
    else if(input$radio == 2){plot(scapshortfrontier, c(1,2,3,4,8))
    }
    else if  (input$radio == 3) {plot(mcaplongfrontier, c(1,2,3,4,8))
    } else if (input$radio ==4) {plot(mcapshortfrontier, c(1,2,3,4,8))}
      else if (input$radio ==5) {plot(lcaplongfrontier, c(1,2,3,4,8))}
      else {plot(lcapshortfrontier, c(1,2,3,4,8))}
    
    })
  }
)

# shiny basic output app

shinyApp(ui = ui, server = server)

library(shiny)

shinyApp(
  ui <- fluidPage(    
    radioButtons("radio", 
                 label = HTML('<FONT color="red"><FONT size="5pt">Welcome to Portfolio Selection Program</FONT></FONT><br> <b>Which Stock do you want?</b>'),
                 c("stats scaplong" = 1, " stats scapshort" = 2, " stats mcaplong"=3, "stats mcapshort" = 4, "stats lcaplong" = 5, "stats lcapshort" = 6)),      
    
    mainPanel(
      verbatimTextOutput("stats")
    )
  ),
  
  
  server = function(input, output) {
    
    
    output$stats <- renderPrint({
      
      # Now let us print some stock data
      if(input$radio == 1){
        
        scaplong
      }
      else if(input$radio == 2){scapshort
      }
      else if  (input$radio == 3) {mcaplong
      } else if (input$radio ==4) {mcapshort}
      else if (input$radio ==5) {lcaplong}
      else {lcapshort}
      
    })
  }
)


shinyApp(ui = ui, server = server)
