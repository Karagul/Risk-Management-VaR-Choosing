#409 Risk Managment HW3
#Author: Ming-Hao Yu
#Date: 2018/04/21

library(data.table)
library(lubridate)
library(tseries)
setwd("C:\\Users\\Ming-Hao\\Desktop\\MFE\\409-Financial Risk Measurement and Management")

if(!exists("hw3")){
    hw3 = fread(file="homework3_data.csv")
    hw3[, date:= as.Date(date, "%Y/%m/%d")]
}

#Q1
start = which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))[1]
VaR = VaRexp = VaRb = VaRexpb = vector()
upper = lower = upperb = lowerb = upperexpb = lowerexpb = vector()
j = 1

for(i in start:length(hw3$date)){
    data = data.table(hw3[1:(i-1), ])
    lambda = 0.995
    n = length(data$date)
    ii = seq(1,n)
    weight = lambda^(n-ii)*(1-lambda)/(1-lambda^n)
    data[, Weight := weight]
    VaR[j] = abs(sort(data$gain)[ceiling(0.01*i)])
    
    data = data[order(data$gain), ]
    ptr = which(cumsum(data$Weight)>=0.01)[1]
    VaRexp[j] = abs(data$gain[ptr])
    
    mu = mean(data$gain)
    sigma = sd(data$gain)
    x = mu+qnorm(0.01)*sigma
    f = exp(-(x-mu)^2/2/sigma^2)/sigma/sqrt(2*pi)
    StdDev = 1/f*sqrt(0.99*0.01/i)
    
    upper[j] = VaR[j] + 2*StdDev
    lower[j] = VaR[j] - 2*StdDev
    j = j+1
}

ylim=c(min(VaR, VaRexp), max(VaR,VaRexp))
plot(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))], 
     y=VaR, xlab="date", ylab="VaR", main="1-day 99%-VaR", type="l", col="red", ylim=ylim)
points(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))], 
       y=VaRexp, col="blue", type="l")
legend("topleft", c("historical VaR", "exponential weighted VaR"), col=c("red", "blue"), cex=0.8, lwd=1)

ylim=c(min(VaR, VaRexp, hw3$gain), max(VaR,VaRexp, hw3$gain))
plot(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))], 
     y=-VaR, xlab="date", ylab="Gain", main="1-day 99%-VaR", type="l", col="red", ylim=ylim)
points(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))], 
       y=-VaRexp, col="blue", type="l")
points(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))],
       y=hw3$gain[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))],
       type="l", col="black")
legend("topleft", c("historical VaR", "exponential weighted VaR", "Gain"), col=c("red", "blue", "black"), cex=0.8, lwd=1)

#calculate exception of VaR/VaRexp
exception1 = which(hw3$gain[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))] < -VaR)+(length(hw3$gain)-length(VaR))
exception2 = which(hw3$gain[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))] < -VaRexp)+(length(hw3$gain)-length(VaR))

#Historical VaR
plot(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))], 
     y=-VaR, xlab="date", ylab="Gain", main="Historical 1-day 99%-VaR", type="l", col="red", ylim=ylim)
points(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))],
       y=hw3$gain[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))],
       type="l", col="black")
points(x=hw3$date[exception1], y=hw3$gain[exception1], col="dark green", type="p", pch=16)
legend("topleft", c("historical VaR", "Gain"), col=c("red", "black"), cex=0.8, lwd=1)
print(c("exception:", length(exception1)))

#exponential VaR
plot(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))], 
     y=-VaRexp, xlab="date", ylab="Gain", main="Exponential weighted 1-day 99%-VaR", type="l", col="blue", ylim=ylim)
points(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))],
       y=hw3$gain[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))],
       type="l", col="black")
points(x=hw3$date[exception2], y=hw3$gain[exception2], col="dark green", type="p", pch=16)
legend("topleft", c("exponential weighted VaR", "Gain"), col=c("blue", "black"), cex=0.8, lwd=1)
print(c("exception:", length(exception2)))

#Q2
j = 1

for(i in start:length(hw3$date)) {
    data = data.table(hw3[1:(i-1), ])
    lambda = 0.995
    n = length(data$date)
    ii = seq(1,n)
    weight = lambda^(n-ii)*(1-lambda)/(1-lambda^n)
    data[, Weight := weight]
    VaR[j] = abs(sort(data$gain)[ceiling(0.01*i)])
    
    data = data[order(data$gain), ]
    v1 = v2 = vector()
    
    #bootstrap historical VaR
    sample1 = matrix(sample(data$gain, size=i*1000, replace=T), nrow=i, ncol=1000)
    v1 = apply(sample1, MARGIN=2, FUN=function(x){return(abs(sort(x)[ceiling(0.01*i)]))})
    VaRb[j] = mean(v1)
    upperb[j] = quantile(v1, 0.975)
    lowerb[j] = quantile(v1, 0.025)
    
    #bootstrap exponential weighted VaR
    sample2 = matrix(sample(data$gain, size=i*1000, replace=T, prob=data$Weight), nrow=i, ncol=1000)
    v2 = apply(sample2, MARGIN=2, FUN=function(x){
        return(abs(sort(x)[ceiling(0.01*i)]))
    })
    VaRexpb[j] = mean(v2)
    upperexpb[j] = quantile(v2, 0.975)
    lowerexpb[j] = quantile(v2, 0.025)
    j = j+1
}

Q2 = matrix(nrow=length(VaRb), ncol=3, c(lowerb, VaRb, upperb))
colnames(Q2) = c("lwr (2.5%)", "VaR", "upr (97.5%)")

plot(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))], 
     y=-VaRb, xlab="date", ylab="Gain", main="Bootstrap Historical 1-day 99%-VaR", type="l", col="red", ylim=ylim)
arrows(x0=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))], 
       y0=-Q2[, 1], y1=-Q2[, 3], lwd=1.5, angle=90, code=3, length=0.05, col="orange")
points(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))],
       y=hw3$gain[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))],
       type="l", col="black")
points(x=hw3$date[exception1], y=hw3$gain[exception1], col="dark green", type="p", pch=16)
legend("topleft", c("historical VaR", "Gain", "95% Confidence Interval"), col=c("red", "black", "orange"), cex=0.8, lwd=1)
head(Q2, 20)

Q2exp = matrix(nrow=length(VaRexpb), ncol=3, c(lowerexpb, VaRexpb, upperexpb))
colnames(Q2exp) = c("lwr (2.5%)", "VaR", "upr (97.5%)")

plot(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))], 
     y=-VaRexpb, xlab="date", ylab="Gain", main="Bootstrap exponential weighted 1-day 99%-VaR", type="l", col="blue", ylim=ylim)
arrows(x0=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))], 
       y0=-Q2exp[, 1], y1=-Q2exp[, 3], lwd=1.5, angle=90, code=3, length=0.05, col="orange")
points(x=hw3$date[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))],
       y=hw3$gain[which(hw3$date >= as.Date("2007-01-01", "%Y-%m-%d"))],
       type="l", col="black")
points(x=hw3$date[exception1], y=hw3$gain[exception1], col="dark green", type="p", pch=16)
legend("topleft", c("historical VaR", "Gain", "95% Confidence Interval"), col=c("blue", "black", "orange"), cex=0.8, lwd=1)
head(Q2exp, 20)

#Q3
jarque.bera.test(hw3$gain)
# p-value is much small than 0.05 which implies that we should reject the Null-hypothesis, therefore, the gain is NOT normal distribution.

#Q4
j = 1
gainNormal = vector()

for(i in start:length(hw3$date)) {
    volatility = sd(hw3$gain[(i-30):(i-1)])
    gainNormal[j] = hw3$gain[i]/volatility
    j = j+1
}

#Normalized gain
jarque.bera.test(gainNormal)

#Original gain
jarque.bera.test(hw3$gain[start:(length(hw3$date))])

#Normalized gain has smaller p-value (0.0001221) than the orighinal one (2.2e-16) so that the normailzed gain is closer to the normal distribution.

#Q5
j = 1
gainAdvancedNormal = vector()

for(i in start:(length(hw3$date)-30)) {
    volatility = sd(hw3$gain[(i+1):(i+30)])
    gainAdvancedNormal[j] = hw3$gain[i]/volatility
    j = j+1
}
jarque.bera.test(gainAdvancedNormal)

#The p-value of gain normalized by the future volatility is similar to the gain normalized by the previous one month volatility, which is also much close to normal distribution compared to the original gain.

#Q6

