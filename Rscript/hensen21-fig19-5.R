##  This file generates Figures 19.5 and 19.6
##  Wage/Experience

#### attention 1/2 ####
## this is the part 1/2 R script of the whole analysis
## so you will run `hensen21-fig19-6.R` forthcoming 

#### Uses data file cps09mar.txt ####

dt_file <- here::here("data/cps09mar.dta")
#dat <- read.table(dt_file)
dat <- haven::read_dta(dt_file)

####data set 1 ####
## subsample of black men with 12 years of education
bf <- (dat[,11]==2)&(dat[,2]==1)&(dat[,4]==12)
dat1 <- dat[bf,]
## earning per hours
y <- as.matrix(log(dat1[,5]/(dat1[,6]*dat1[,7])))
x <- as.matrix(dat1[,1]-dat1[,4]-6)
n <- length(y)

dt_cps1 <-  tibble(X=x, Y=y)

#### basic plot ####
fsize <- 16
p0 <- ggplot() +
  geom_point(aes(X, Y),data = dt_cps1, pch=21) +
  labs(x= "X职业年数", y ="log(Y)时均工资的对数") +
  scale_x_continuous(breaks = seq(0,70,10), limits = c(0,65)) +
  scale_y_continuous(breaks = seq(0,4,1), limits = c(0,4)) +
  theme_bw() +
  theme(text = element_text(size = fsize))


#### Reference Rule ####
## window range
sx1 <- 0
sx2 <- 40
x1 <- matrix(1,n,1)
zz <- cbind(x1,x,x^2,x^3,x^4)
beta <- solve((t(zz)%*%zz),(t(zz)%*%y))
xtrim <- (x<=sx2)*(x>=sx1)
b <- mean(((beta[3]+x*3*beta[4]+(x^2)*6*beta[5])^2)*xtrim)
e <- y - zz%*%beta
sig <- (sum(e^2))/(n-5)
hrot <- 0.58*(((sx2-sx1)*sig/n/b)^.2)


#### CV Bandwidth Selection ####
g <- 200
h1 <- hrot/3
h2 <- 3*hrot
hh <- seq(h1,h2,(h2-h1)/g)
hn <- length(hh)

LL <- matrix(0,n,hn)
for (i in 1:hn){
  hi <- hh[i]
  for (j in 1:n){
    xj <- x-x[j]
    k <- dnorm(xj/hi)
    k[j] <- 0
    z <- cbind(x1,xj)
    zk <- z*(k%*%cbind(1,1))
    beta <- solve(t(zk)%*%z,t(zk)%*%y)
    LL[j,i] <- (y[j]-beta[1])^2
  }
}

LL2 <- LL*(xtrim%*%matrix(1,1,hn))
cvLL <- colMeans(LL2)
i <- which.min(cvLL)
hLL <- hh[i]
hLL_tex <-paste0("$h_{CV}=$",number(hLL,0.0001),"")
LLmin <- min(cvLL)

# combine result as data.frame
tb_cvc <- tibble(
  h_tune = hh, 
  cv_LL = cvLL)


#### plot 19.5-a Cross-Validation Criterion ####

lwd <- 0.8
lwadd <- 0.2

p_cvc <- ggplot(aes(x = h_tune),data = tb_cvc) +
  geom_line(aes(y = cv_LL),
            lty = "dashed", color = "blue",
            lwd = lwd) +
  labs(x= "谱宽h", y ="交叉验证准则函数值CV(h)") +
  scale_x_continuous(#expand = expansion(add = c(0, .06)),
    breaks = seq(1,12,1), 
    limits = c(1,12)) +
  scale_y_continuous(
    expand = expansion(add = c(0,0)),
    breaks = seq(0.2290,0.2310,0.0005), 
    limits = c(0.2290,0.2310),
    labels = scales::number_format(accuracy = 0.0001)
                     ) +
  geom_segment(aes(x=hLL, xend = hLL,
                   y= LLmin, yend = c(0.2290)),
               arrow = arrow(length = unit(0.1,"cm"),
                             type = "closed"),
               lty = "dashed", color= "gray", lwd=lwd) +
  geom_text(aes(x = hLL+0.5, y=LLmin-0.0004), 
            label =latex2exp::TeX(hLL_tex),
            #parse=TRUE,
            color= "blue", size=4)  +
  theme_bw() +
  theme(text = element_text(size = fsize))


#### CEF LL Regression Estimation ####
g <- 201
xg <- seq(sx1,sx2,(sx2-sx1)/(g-1))
m1 <- matrix(0,g,1)
m2 <- matrix(0,g,1)
for (j in 1:g){
  xj <- x-xg[j]
  z <- cbind(x1,xj)
  z1 <- z*(dnorm(xj/hrot)%*%cbind(1,1))
  z2 <- z*(dnorm(xj/hLL)%*%cbind(1,1))
  beta1 <- solve(t(z1)%*%z,t(z1)%*%y)
  beta2 <- solve(t(z2)%*%z,t(z2)%*%y)
  m1[j,1] <- beta1[1]
  m2[j,1] <- beta2[1]
}

# combine result as data.frame
tb_mxh <- tibble(xg = xg,
                 mx1=as.vector(m1), 
                 mx2=as.vector(m2))

#### plot 19.5-b Local Linear Regression ####

# basic plot
p00 <- ggplot() +
  geom_point(aes(X, Y),data = dt_cps1, pch=21,alpha =0.1) +
  labs(x= "X职业年数", y ="log(Y)时均工资的对数") +
  scale_x_continuous(breaks = seq(0,40,5), limits = c(0,40)) +
  scale_y_continuous(breaks = seq(2,2.7,0.1), limits = c(2,2.7)) +
  theme_bw() +
  theme(text = element_text(size = fsize))

# plot mx with h_rot
p_mxh1 <- p00 +
  geom_line(aes(x = xg, y = m1, 
                color="m1", lty="m1"),
            lwd = lwd+lwadd,
            data = tb_mxh) +
  theme_bw() +
  scale_color_manual(
    name="",
    breaks = c("m1"),
    labels = c(expression(m(x):h[ROT])),
    values=c("green"))+
  scale_linetype_manual(
    name="",
    breaks = c("m1"),
    labels = c(expression(m(x):h[ROT])),
    values=c("dotted"))+
  theme(legend.position = "right",
        text = element_text(size = fsize))

# plot mx with h_CV
p_mxh2 <- p00 +
 geom_line(aes(x = xg, y = m2, 
                color="m2", lty="m2"),
            lwd = lwd+lwadd,
            data = tb_mxh) +
  scale_color_manual(
    name="",
    breaks = c("m2"),
    labels = c(expression(m(x):h[CV])),
    values=c("blue"))+
  scale_linetype_manual(
    name="",
    breaks = c("m2"),
    labels = c(expression(m(x):h[CV])),
    values=c("dashed"))+
  theme(legend.position = "right",
        text = element_text(size = fsize))

# plot mx with both h selections
p_mxh <- p00 +
  geom_line(aes(x = xg, y = m1, 
                color="m1", lty="m1"),
            lwd = lwd+lwadd,
            data = tb_mxh) +
  geom_line(aes(x = xg, y = m2, 
                color="m2", lty="m2"), 
            lwd = lwd+lwadd,
            data = tb_mxh) +
  scale_color_manual(
    name="",
    breaks = c("m1", "m2"),
    labels = c(expression(m(x):h[ROT]),expression(m(x):h[CV])),
    values=c("green", "blue"))+
  scale_linetype_manual(
    name="",
    breaks = c("m1", "m2"),
    labels = c(expression(m(x):h[ROT]),expression(m(x):h[CV])),
    values=c("dotted", "dashed"))+
  theme(legend.position = "right",
        text = element_text(size = fsize))

