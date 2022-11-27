# R-package Change Point: https://cran.r-project.org/web/packages/changepoint/changepoint.pdf
library(changepoint)
# change in mean
y=c(rnorm(100,0,1),rnorm(100,5,1)) # continuous data.
ansmean=cpt.mean(y,penalty="MBIC",pen.value=0,method="PELT",Q=5,test.stat="Normal",class=TRUE,
                 param.estimates=TRUE,minseglen=1)
plot(ansmean,cpt.col='blue')
print(ansmean)

#Created Using changepoint version 2.2.3
#Changepoint type      : Change in mean
#Method of analysis    : PELT
#Test Statistic  : Normal
#Type of penalty       : MBIC with value, 15.89495
#Minimum Segment Length : 1
#Maximum no. of cpts   : Inf
#Changepoint Locations : 100

# Page 6. For change points in mean.
ansmean=cpt.mean(cp_Y$T_k,penalty="MBIC",pen.value=0.05,method="PELT",test.stat="Poisson",
            class=TRUE,param.estimates=TRUE,shape=1,minseglen=1)
plot(cp_Y$T_k) # observed data
plot(ansmean,cpt.col='blue') # change points
print(ansmean)

# Page 6. For change points in mean.
ansmean=cpt.mean(cp_Y$T_k,penalty="MBIC",pen.value=0.05,method="BinSeg",test.stat="Poisson",
                 class=TRUE,param.estimates=TRUE,shape=1,minseglen=1)
plot(cp_Y$T_k) # observed data
plot(ansmean,cpt.col='blue') # change points
print(ansmean)
#Changepoint type      : Change in mean and variance
#Method of analysis    : PELT
#Test Statistic  : Poisson
#Type of penalty       : MBIC with value, 6.907755
#Minimum Segment Length : 1
#Maximum no. of cpts   : Inf
#Changepoint Locations : 1 2 4 8 9

# Using two common change point detection method to detect a change in the mean of our infection data,
# PELT and BinSeg detected change point locations in the data at 1 2 4 8 9. Segmenting the data
# in this way would result in very biased transmission rate estimates.


