### This file generates Figure 21.1a
### Regression Discontinuity Design


##### show outcome ####
require(here)

x1 <- as.matrix(seq(0,1,by=0.01))
x2 <- as.matrix(seq(1,2.5,by=0.01))

b1 <- (x1/3)^(1/6)
b2 <- (x2/3)^(1/6)

t1 <- (x1/3)^(1/3)
t2 <- (x2/3)^(1/3)

wd <- 2.4
size_point <- 3
size_text <- 1.5

pic_path <- here("pic/HANSEN21-1a-outcome-expect.png")
#pdf(here("pic/HANSEN21-1a-outcome.pdf"),family="Helvetica",colormodel='cmyk',width=4.5,height=4.5)
#postscript(here("pic/HANSEN21-1a-outcome.eps"),paper="special",width=4.5,height=4.5,horizontal=FALSE,onefile=FALSE,family="Helvetica",colormodel='cmyk')
png(pic_path, width = 450, height = 450)

plot(x1,b1,type="l",lty=1,col="blue",xaxs="i",yaxs="i",,ylab="Output Variable Y",xlab="Running Variable X",bty="n",xlim=c(0,2.5),ylim=c(0,1.1),yaxt="n",xaxt="n",lwd=wd,yaxt="n")
axis(side=1,lwd=wd)
axis(side=2,seq(0,1.2,.2),lwd=wd)
lines(x1,t1,lty=2,lwd=wd, col='green')
lines(x2,b2,lty=2,lwd=wd, col='blue')
lines(x2,t2,lwd=wd, col = "green")
lines(c(1,1),c(b2[1],t2[1]),lwd=3, col = "red")
abline(h=0,lwd=wd)
abline(v=0,lwd=wd)
lines(c(1,1),c(0,1.2),lty=2,lwd=wd)
text(1.05,.1,"c", cex = size_text)
text(1.75,.1,"X>c", cex = size_text, col = "green")
text(0.5,.1,"X<c", cex = size_text, col = "blue")
text(.95,0.76,expression(bar(tau)), cex =  size_text)
text(1.8,.76,expression(Y=m[1](X)), cex = size_text)
text(.45,.82,expression(Y=m[0](X)), cex = size_text)

dev.off()

##### show conditional expect of outcome  ####
x1 <- as.matrix(seq(0,1,by=0.02))
x2 <- as.matrix(seq(1,2.5,by=0.02))

b1 <- (x1/3)^(1/6) +rnorm(51,sd=0.01)
b2 <- (x2/3)^(1/6) + rnorm(76,sd =0.02)

t1 <- (x1/3)^(1/3) +rnorm(51,sd=0.01)
t2 <- (x2/3)^(1/3) + rnorm(76,sd =0.02)

wd <- 2.4
size_point <- 3
size_text <- 1.5

pic_path <- here("pic/HANSEN21-1a-outcome.png")
#pdf(here("pic/HANSEN21-1a-outcome-expect.pdf"),family="Helvetica",colormodel='cmyk',width=4.5,height=4.5)
#postscript(here("pic/HANSEN21-1a-outcome-expect.eps"),paper="special",width=4.5,height=4.5,horizontal=FALSE,onefile=FALSE,family="Helvetica",colormodel='cmyk')
png(pic_path, width = 600, height = 450)

plot(x1,b1,type="p",pch=21,col="blue",xaxs="i",yaxs="i",ylab="Output Variable Y",xlab="Running Variable X",bty="n",xlim=c(0,2.5),ylim=c(0,1.1),yaxt="n",xaxt="n",lwd=wd,yaxt="n")
axis(side=1,lwd=wd)
axis(side=2,seq(0,1.2,.2),lwd=wd)
points(x1,t1,pch=21,col=rgb(red = 0,green =  255,blue = 0, alpha = 80,maxColorValue=255))
points(x2,b2,pch=21,col=rgb(red = 0,green = 0, blue = 255, alpha = 80,maxColorValue=255))
points(x2,t2,pch =21, col = "green")
lines(c(1,1),c(b2[1],t2[1]),lwd=3, col = "red")
abline(h=0,lwd=wd)
abline(v=0,lwd=wd)
lines(c(1,1),c(0,1.2),lty=2,lwd=wd)
text(1.05,.1,"c", cex = size_text)
text(1.75,.1,"X>c", cex = size_text, col = "green")
text(0.5,.1,"X<c", cex = size_text, col = "blue")
#text(.95,0.76,expression(bar(tau)), cex =  size_text)
text(1.8,.76,expression(Y[1](X>c[0])), cex = size_text)
text(.45,.82,expression(Y[0](X<c[0])), cex = size_text)

dev.off()