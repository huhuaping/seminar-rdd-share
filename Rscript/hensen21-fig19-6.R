##  This file generates Figures 19.5 and 19.6
##  Wage/Experience

#### attention 2/2 ####

## this is the part 2/2 R script
## so you should run `hensen21-fig19-5.R` firstly 

#### ROT Residual Estimation ####
h <- hrot
e2 <- matrix(0,n,1)
for (j in 1:n){
  xj <- x-x[j]
  k <- dnorm(xj/h)
  k[j] <- 0
  z <- cbind(x1,xj)
  zk <- z*(k%*%cbind(1,1))
  beta <- solve(t(zk)%*%z,t(zk)%*%y)
  e2[j,1] <- (y[j]-beta[1])^2
}

####  Variance ROT Bandwidth ####
beta <- solve((t(zz)%*%zz),(t(zz)%*%e2))
b <- mean(((beta[3]+x*3*beta[4]+(x^2)*6*beta[5])^2)*xtrim)
e <- e2 - zz%*%beta
sig <- (sum(e^2))/(n-5)
hrotv <- 0.58*(((sx2-sx1)*sig/n/b)^.2)


#### Variance CV Bandwidth Selection ####
h1 <- 2
h2 <- 40
hh <- seq(h1,h2,(h2-h1)/g)
hn <- length(hh)
LLV <- matrix(0,n,hn)
NWV <- matrix(0,n,hn)
for (i in 1:hn){
  hi <- hh[i]
  for (j in 1:n){
    xj <- x-x[j]
    k <- dnorm(xj/hi)
    k[j] <- 0
    mu <- (t(k)%*%e2)/sum(k)
    NWV[j,i] <- (e2[j]-mu)^2
    z <- cbind(x1,xj)
    zk <- z*(k%*%cbind(1,1))
    beta <- solve(t(zk)%*%z,t(zk)%*%e2)
    LLV[j,i] <- (e2[j]-beta[1])^2
  }
}
cvNWV <- colMeans(NWV*(xtrim%*%matrix(1,1,hn)))
cvLLV <- colMeans(LLV*(xtrim%*%matrix(1,1,hn)))
i1 <- which.min(cvNWV)
i2 <- which.min(cvLLV)
hNWV <- hh[i1]
hLLV <- hh[i2]
hNWV_tex <-paste0("$h_{NW}=$",number(hNWV,0.1),"")
hLLV_tex <-paste0("$h_{CV}=$",number(hLLV,0.1),"")

NWminV <- min(cvNWV)
LLminV <- min(cvLLV)

# combine all result
tbl_cvcv <- tibble(
  h_tune = hh, 
  cv_NW = cvNWV,
  cv_LL = cvLLV)

#### plot 19.6-a Cross-Validation Criterion ####
## Cross-Validation for Conditional Variance

p_cvcv <- ggplot(aes(x = h_tune),data = tbl_cvcv) +
  geom_line(aes(y = cv_NW, lty = "CV_NW", color = "CV_NW"),
            lwd = lwd) +
  geom_line(aes(y = cv_LL, lty = "CV_LL", color = "CV_LL"),
            lwd = lwd) +
  labs(x= "谱宽h", y ="交叉验证准则函数值CV(h)") +
  scale_x_continuous(#expand = expansion(add = c(0, .06)),
    breaks = seq(0,40,5), 
    limits = c(0,40)) +
  scale_y_continuous(expand = expansion(add = c(0,0)),
                     breaks = seq(1.200, 1.216,0.004), 
                     limits = c(1.200, 1.216),
                     labels = scales::number_format(accuracy = 0.001)
                     ) +
  theme_bw() +
  geom_segment(aes(x=hNWV, xend = hNWV,
                   y= NWminV, yend = 1.2),
               arrow = arrow(length = unit(0.1,"cm"),
                             type = "closed"),
               lty = "dotted", color= "gray", lwd=lwd)  +
  geom_text(aes(x = hNWV-2, y=NWminV-0.002), 
            label =latex2exp::TeX(hNWV_tex),
            color= "orange", size=3) +
  geom_segment(aes(x=hLLV, xend = hLLV,
                   y= LLminV, yend = 1.2),
               arrow = arrow(length = unit(0.1,"cm"),
                             type = "closed"),
               lty = "dashed", color= "gray", lwd=lwd) +
  geom_text(aes(x = hLLV-2, y=LLminV-0.002), 
            label =latex2exp::TeX(hLLV_tex),
            color= "blue", size=3)  +
  scale_color_manual(
    name="",
    breaks = c("CV_NW", "CV_LL"),
    values=c("orange", "blue"))+
  scale_linetype_manual(
    name="",
    breaks = c("CV_NW", "CV_LL"),
    values=c("dotted", "dashed"))+
  theme(legend.position = "right")

#### Variance Estimation ####
se <- matrix(0,g,1)
u <- cbind(e2,e2)
for (j in 1:g){
  xj <- x-xg[j]
  z <- cbind(x1,xj)
  K <- dnorm(xj/hrot)
  z1 <- z*cbind(K,K)
  ZKZ <- solve(t(z1)%*%z)
  ZK2 <- t(z1)%*%(z1*u)
  V1 <- ZKZ%*%ZK2%*%ZKZ
  se[j,1] = sqrt(V1[1,1])
}

L1 <- m1-1.96*se
U1 <- m1+1.96*se

# combine all result as data.frame
tbl_band <- tibble(
  xg =xg,
  mx = m1,
  s=as.vector(se),
  s2=se^2,
  lwr = L1,
  upr = U1
)

#### Plot 19.6-b Confidence Bands ####

p_band <- p00 +
  geom_line(aes(x = xg, y = mx),
            data = tbl_band,
            color = "blue", lty = "solid", lwd = lwd) +
  geom_ribbon(aes(x=xg,ymin = lwr, ymax = upr),
              data = tbl_band,
              alpha = 0.2)


