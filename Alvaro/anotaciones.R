https://www.statmethods.net/graphs/density.html 
x <- data.frame(exprseset)[, c(12)]
d <- density(x)
plot(d)
h<- hist(x, lty = 1, col = mycols,
         main = "Intensity distributions", density = TRUE)
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col=mycols[8], lwd=5)
lines(xfit, yfit, col=mycols[2], lwd=1)

for (i in 1:length(CELfiles)) {
  x <- data.frame(exprseset)[, i]
  h<- hist(x, lty = 1, col = "white",
           main = "Intensity distributions", plot = FALSE)
  xfit<-seq(min(x),max(x),length=40)
  yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
  yfit <- yfit*diff(h$mids[1:2])*length(x)
  lines(xfit, yfit, col=mycols[i], lwd=2)
}

for (i in 1:length(CELfiles)) {
  x <- data.frame(exprseset)[, i]
  d <- density(x)
  plot(d)
}
