rm(list=objects())
setwd("/home/achraf/Desktop/")

library(dplyr)
library(xts)
library(mgcv)
library(magrittr)
library(forecast)


df = read.csv("sample_superstore.csv")

# On ne garde que les ventes de fourniture
df = df[df$Category == "Furniture", ]

start_date = min(as.Date(df$Order.Date, "%d/%m/%Y"))
end_date = max(as.Date(df$Order.Date, "%d/%m/%Y"))

start_date = strptime(start_date, "%Y-%m-%d")
end_date = strptime(end_date, "%Y-%m-%d")
Date = seq.POSIXt(from=start_date, to=end_date, by="day")

Order.Date = strptime(df$Order.Date, "%d/%m/%Y")
X = xts(df$Quantity, order.by=Order.Date)

X.sum = apply.daily(X, sum)
plot(X.sum, type="l")
X.mean = apply.monthly(X.sum, mean)
plot(X.mean, type="l")

df_new = data.frame("Date"=index(X.sum), "Sum"=as.numeric(X.sum))

df_new %<>%
  mutate(Date = as.Date(df_new$Date)) %>%
  complete(Date=seq.Date(as.Date(min(df_new$Date)), as.Date(max(df_new$Date)), by="day" ))

for(i in c(1:length(df_new$Sum))) {
  if(is.na(df_new$Sum[i])) {
      df_new$Sum[i] = 0
  }
}

X.sum = xts(df_new$Sum, order.by=df_new$Date)
plot(X.sum)

for(i in c(1:length(df_new$Sum))) {
  if(df_new$Sum[i] == 0 ) {
    if(i > 1) {
      df_new$Sum[i] = df_new$Sum[i-1]
    } else {
      df_new$Sum[i] = 0 
    }
  }
}

X.before = xts(df_new$Sum, order.by=df_new$Date)

### TRAIN TEST SPLIT
split_date = strptime("01/01/2017", format = "%d/%m/%Y")
train_set = window(X.before, start=start_date, end=split_date)
test_set = window(X.before, start=split_date, end=end_date)


Date = seq.POSIXt(from=start_date, to=split_date, by="day")
X = xts(train_set, order.by=Date)
n = length(X)
t = c(1:n)

# REGRESSION
### Linéaire
reg = lm(X ~ t)
summary(reg)
ychap.lm = xts(reg$fitted, order.by=Date)
plot(X, type="l", col="gray")
lines(reg$fitted.values, col="red")

### SPLINES
g = gam(X~s(t, k=15))
summary(g)
ychap.gam = xts(g$fitted, order.by=Date)
plot(X, type='l')
lines(xts(ychap.gam, order.by=Date), col='red', lwd= 6)

res = seasonality(data=X , trend=ychap.gam, K=10, period=365)
plot(X, type="l")
lines(res$season, col="red", lwd=2)
plot(X - ychap.gam - res$season)

### GAUSSIEN
h = 365
x = seq(1,max(t),length=n)
W = matrix(unlist(lapply(x,
  function(x){
   dnorm(x-t,0,sd=sqrt(h/2))/sum(dnorm(x-t,0,sd=sqrt(h/2)))
  })
),ncol=n,nrow=n,byrow=F)

ychap.kernel = colSums(as.numeric(X) * W)
ychap.kernel = xts(ychap.kernel,order.by=Date)

plot(X,type='l', col="gray")
lines(ychap.kernel, col='red')

### POLYNOME LOCAUX
lo = loess(X~t, degree=2, span=0.9)
ychap.lo = xts(lo$fitted, order.by=Date)
plot(X, type="l")
lines(ychap.lo,col="red", lwd=5)


# Saisonalité
seasonality = function(data, trend=NA, period=365, K=10) {
  if(is.na(trend)) {
    x = data
  } else {
    x = data - trend
  }
  
  acf(x, lag.max=5000)
  w = 2*pi/period
  fourier = cbind(cos(w*t), sin(w*t))
  
  for(i in c(2:K)){
    fourier = cbind(fourier, cos(i*w*t), sin(i*w*t))
  }
  reg = lm(x ~ fourier)
  summary(reg)
  season = xts(reg$fitted, order.by=Date)
  res = list()
  res$season = season
  res$reg = reg
  return(res)
}

res = seasonality(data=X, K=10)
plot(X, type="l")
lines(res$season, col="red", lwd=2)

# LISSAGE
### SIMPLE
simple_smooth = function(x, alpha){
  s=x
  for(i in c(2:length(x))){
    s[i]= alpha*x[i]+(1-alpha)*s[i-1]
  }
  return(s)
}
res = xts(simple_smooth(X, 0.1), order.by=Date)

plot(X,type="l", col="gray")
lines(res, col='red')

### DOUBLE
double_smooth = function(x, alpha)
{
  xsmooth = x
  l = array(x[1], dim=length(x))
  b = array(x[2]-x[1], dim=length(x))
  
  for(i in c(2:length(x)))
  {
    l[i] = xsmooth[i-1] + (1-(1-alpha)^2)*(x[i]-xsmooth[i-1])
    b[i] = b[i-1]+alpha^2 * (x[i]-xsmooth[i-1])
    xsmooth[i] = l[i]+b[i]
  }
  
  res = list()
  res$smooth = xsmooth
  res$l = l
  res$b = b
  return(res)
}


alpha = seq(0.05, 0.95, length=100)
forecast = lapply(alpha, double_smooth, x=as.numeric(X))
erreur = unlist(
  lapply(forecast,
         function(x){mean((tail(X,n-1)-head(x$smooth,n-1))^2)}))
plot(alpha,erreur,type='l')

X.smooth = double_smooth(as.numeric(X), alpha[which.min(erreur)])
plot(X,type='l')
lines(xts(X.smooth$smooth, order.by=Date), col='red')

plot(X.smooth$l, type='l',ylim=range(X.smooth$l,X.smooth$b),col='blue')
lines(X.smooth$b, col='red')

res = seasonality(data=X.smooth$smooth, K=3)
plot(X, type="l")
lines(xts(X.smooth$smooth, order.by=Date), type="l", col="green")
lines(res$season, col="red", lwd=2)


### LISSAGE DOUBLE DE HOLT-WINTERS
double_holt_winters = function(x,alpha,beta,delta,T)
{
  xsmooth=x
  l<-array(x[2],dim=length(x))
  b<-array(x[2]-x[1],dim=length(x))
  s<-array(x[1],dim=length(x))
  
  for(i in c(2:length(x)))
  {
    l[i]<-alpha*(x[i]-s[max(i-T,1)])+(1-alpha)*(l[i-1]+b[i-1])
    b[i]<-beta*(l[i]-l[i-1])+(1-beta)*b[i-1]
    s[i]<-delta*(x[i]-l[i])+(1-delta)*s[max(i-T,1)]
    xsmooth[i]<-l[i]+b[i]+s[i]
  }
  
  res<-list()
  res$smooth<-xsmooth
  res$l=l
  res$b<-b
  res$s<-s
  return(res)
}


alpha = 0.2
beta = 0.2
delta = 0.2
period = 10
X.seas.exp.mooHW = double_holt_winters(as.numeric(X), alpha,
                                        beta, delta, period)


plot(X, type="l")
lines(xts(X.seas.exp.mooHW$smooth - X.seas.exp.mooHW$s, order.by=Date), col="red")

seasonality(data = X.seas.exp.mooHW$smooth)
plot(X, type="l")
lines(res$season, col="red", lwd=2)


lines(xts(X.seas.exp.mooHW$s, order.by=Date), col="green")



### PREDICTIONS
predict.smooth = function(Xsmooth,inst,horizon,smooth.type="double")
{
  
  if(smooth.type=="simple")
  {
    n<-length(Xsmooth)
    prev<-c(Xsmooth[1:inst],rep(Xsmooth[inst],horizon))
  }
  
  if(smooth.type=="double")
  {
    n<-length(Xsmooth$smooth)
    prev<-c(Xsmooth$smooth[1:inst],Xsmooth$l[inst]+Xsmooth$b[inst]*c(1:horizon))
  }
  
  if(smooth.type=="holt-winter")
  {
      
  }
  return(prev)
}

smooth.double.forecast = predict.smooth(X.smooth, inst=length(X), horizon=365,
                                        smooth.type="double")
smooth.simple.forecast = predict.smooth(X.smooth, inst=length(X), horizon=365,
                                        smooth.type="simple")

plot(as.numeric(X), col="gray", type="l")
lines(smooth.double.forecast, col="red")
lines(smooth.simple.forecast, col="green")
abline(v=length(X),lty='dashed')


### SAISONALITE CYCLIQUE
cycle = c(rep(c(1:(length(x)/3)),3),1)
cycle = cycle[1:length(cycle)-1]
plot(cycle)

g = gam(X~s(cycle,k=20, bs='cc'))
summary(g)
ychap.gam.season = xts(g$fitted, order.by=Date)
plot(X, type='l')
lines(ychap.gam.season, col='red', lwd = 2)

acf(X - ychap.gam.season)
acf(X - res$season)

plot(ychap.gam.season + ychap.gam)
plot(res$season, col="green")
lines(ychap.gam.season + ychap.gam, col="red")
