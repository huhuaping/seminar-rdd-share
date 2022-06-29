
##  This file generates Figure 21.1b, Table 21.1 Baseline, and equation (21.5)
##  Head Start Regression Discontinuity

####  The data file LM2007.dta is used ####

library(haven)
library(here)
dt_path <- here::here("data/LM2007.dta")
dat <- read_dta(dt_path)
y <- dat$mort_age59_related_postHS
x <- dat$povrate60

c <- 59.1984

dt_head <- tibble(X= x,
                  Y = y) %>%
  mutate(D= ifelse(X>=c,1,0))

#### basic plot ####
lwd <- 0.6
fsize <- 16
p0 <- ggplot() +
  geom_point(aes(X, Y),data = dt_head, pch=21,alpha=0.6) +
  geom_vline(aes(xintercept=c),
             color = "orange", lty = "dotted", lwd=lwd) +
  labs(x= "X贫困线", y ="Y儿童死亡率") +
  scale_y_continuous(breaks = seq(0,80,20), limits = c(0,80)) +
  scale_x_continuous(breaks = seq(10,90,10), limits = c(10,90)) +
  theme_bw() +
  theme(text = element_text(size = fsize))

#### RDD Helper function ####

# specify parameters and basic calculation
h <- 8
T <- as.numeric(x >= c)
y1 <- y[T==0]
y2 <- y[T==1]
x1 <- x[T==0]
x2 <- x[T==1]
n1 <- length(y1)
n2 <- length(y2)
n <- n1+n2

# set bins for both side
g1 <- seq(15,59.2,.2)
g2 <- seq(59.2,82,.2)
G1 <- length(g1)
G2 <- length(g2)

# Helper function for  triangle kernel function
TriKernel <- function(x) {
  s6 <- sqrt(6)
  ax <- abs(x)
  f <- (1-ax/s6)*(ax < s6)/s6
  return(f)
}

# Helper function for LLR estimation
LL_EST <- function(y,x,g,h) {
  G <- length(g)
  m <- matrix(0,G,1)
  z <- matrix(1,length(y),1)
  for (j in 1:G){
    xj <- x-g[j]
    K <- TriKernel(xj/h)
    zj <- cbind(z,xj)
    zz <- t(cbind(K,xj*K))
    beta <- solve(zz%*%zj,zz%*%y)
    m[j] <- beta[1]
  }
  return(m)
}

# Helper function for forecast residuals
LL_Residual <- function(y,x,h) {
  n <- length(y)
  e <- matrix(0,n,1)
  z <- matrix(1,n,1)
  for (j in 1:n){
    xj <- x-x[j]
    K <- TriKernel(xj/h)
    K[j] <- 0
    zj <- cbind(z,xj)
    zz <- t(cbind(K,xj*K))
    beta <- solve(zz%*%zj,zz%*%y)
    e[j] <- y[j] - beta[1]
  }
  return(e)
}

# Helper function for standard error
LL_SE <- function(y,x,g,h) {
  G <- length(g)
  s <- matrix(0,G,1)
  z <- matrix(1,length(y),1)
  e <- LL_Residual(y,x,h)     # used inner
  for (j in 1:G){
    xj <- x-g[j]
    K <- TriKernel(xj/h)
    zj <- cbind(z,xj)
    zz <- cbind(K,xj*K)
    ZKZ <- solve(t(zz)%*%zj)
    Ke <- K*e
    ze <- cbind(Ke,xj*Ke)
    V <- ZKZ%*%crossprod(ze)%*%ZKZ
    s[j] = sqrt(V[1,1])
  }
  return(s)
}

#### RDD Estimation ####

# estimate m(x) for both side
m1 <- LL_EST(y1,x1,g1,h)
m2 <- LL_EST(y2,x2,g2,h)

# obtain standard error and interval for both side
s1 <- LL_SE(y1,x1,g1,h)
s2 <- LL_SE(y2,x2,g2,h)
L1 <- m1-1.96*s1
U1 <- m1+1.96*s1
L2 <- m2-1.96*s2
U2 <- m2+1.96*s2

# rdd effect 
theta_lft <- m1[G1]    # m(x) for left cut-off
theta_rgt <- m2[1]     # m(x) for right cut-off
theta <- m2[1]-m1[G1]  # the ATE rdd effect 
v_g1 <- s1[G1]^2
v_g2 <- s2[1]^2

# t test
setheta <- sqrt(s1[G1]^2 + s2[1]^2)
tstat <- theta/setheta
pvalue <- 2*(1-pnorm(abs(tstat)))


#### tibble of mx estimate related ####
tbl_mx1 <- tibble(group = "control",
                  xg = g1,
                  mx=as.vector(m1),
                  s = as.vector(s1),
                  lwr = as.vector(L1), 
                  upr = as.vector(U1))

tbl_mx2 <- tibble(group = "treat",
                  xg = g2,
                  mx=as.vector(m2),
                  s = as.vector(s2),
                  lwr = as.vector(L2), 
                  upr = as.vector(U2))


tb_mxh <- bind_rows(tbl_mx1, tbl_mx2)

#### tibble of estimate result ####
tbl_theta01 <- tibble(
  model = "baseline",
  pars = c("theta"),
  est = c(theta),
  se  = c(setheta)
)


#### plot CEF mx ####

lwd <- 0.6
lwadd <- 0.2

# basic plot
p00 <- ggplot() +
  geom_point(aes(X, Y),data = dt_head, pch=21,alpha=0.3) +
  geom_vline(aes(xintercept=c),
             color = "red", lty = "dotted", lwd=lwd) +
  labs(x= "X贫困线", y ="Y儿童死亡率") +
  scale_y_continuous(breaks = seq(0,5,1), limits = c(0,5)) +
  scale_x_continuous(breaks = seq(10,85,10), limits = c(10,85)) +
  theme_bw() +
  theme(text = element_text(size = fsize))

# mx plot for left part
p_mxh1 <- p00 +
  geom_line(aes(x = xg, y = mx, 
                color="m1", lty="m1"),
            lwd = lwd+lwadd,
            data = base::subset(tb_mxh,group=="control")) +
  #theme_bw() +
  scale_color_manual(
    name="",
    breaks = c("m1"),
    labels = c(expression(m(x):control)),
    values=c("purple"))+
  scale_linetype_manual(
    name="",
    breaks = c("m1"),
    labels = c(expression(m(x):control)),
    values=c("solid"))+
  theme(legend.position = "right",
        text = element_text(size = fsize))


# mx plot for right part
p_mxh2 <- p00 +
  geom_line(aes(x = xg, y = mx, 
                color="m2", lty="m2"),
            lwd = lwd+lwadd,
            data = base::subset(tb_mxh,group=="treat")) +
  scale_color_manual(
    name="",
    breaks = c("m2"),
    labels = c(expression(m(x):treat)),
    values=c("blue"))+
  scale_linetype_manual(
    name="",
    breaks = c("m2"),
    labels = c(expression(m(x):treat)),
    values=c("solid"))+
  theme(legend.position = "right",
        text = element_text(size = fsize))

# mx plot for both part
p_mxh <- p00 +
  geom_line(aes(x = xg, y = mx, 
                color="m1", lty="m1"),
            lwd = lwd+lwadd,
            data = base::subset(tb_mxh,group=="control")) +
  geom_line(aes(x = xg, y = mx, 
                color="m2", lty="m2"), 
            lwd = lwd+lwadd,
            data = base::subset(tb_mxh,group=="treat")) +
  scale_color_manual(
    name="",
    breaks = c("m1", "m2"),
    labels = c(expression(m(x):control),expression(m(x):treat)),
    values=c("purple", "blue"))+
  scale_linetype_manual(
    name="",
    breaks = c("m1", "m2"),
    labels = c(expression(m(x):control),expression(m(x):treat)),
    values=c("solid", "solid"))+
  theme(legend.position = "right")

#### Plot Confidence Bands ####

# band plot for left side
p_band_left <- p_mxh1 +
  geom_ribbon(aes(x=xg,ymin = lwr, ymax = upr),
              data = base::subset(tb_mxh,group=="control"),
              alpha = 0.2) 

# band plot for right side
p_band_right <- p_mxh2 +
  geom_ribbon(aes(x=xg,ymin = lwr, ymax = upr),
              data = base::subset(tb_mxh,group=="treat"),
              alpha = 0.2) 

# band plot for both side
p_band <- p_mxh +
  geom_ribbon(aes(x=xg,ymin = lwr, ymax = upr),
              data = base::subset(tb_mxh,group=="control"),
              alpha = 0.2) +
  geom_ribbon(aes(x=xg,ymin = lwr, ymax = upr),
              data = base::subset(tb_mxh,group=="treat"),
              alpha = 0.2) +
  geom_text(aes(x = c-3, y = theta_lft, 
                label=number(theta_lft, 0.0001)),
            color = "purple")+
  geom_text(aes(x = c+3, y = theta_rgt, 
                label=number(theta_rgt, 0.0001)),
            color = "blue") +
  geom_segment(aes(x=c, xend = c,
                   y= theta_lft, yend = theta_rgt),
               arrow = arrow(length = unit(0.1,"cm"),
                             type = "closed"),
               lty = "solid", color= "orange", lwd=lwd)  +
  geom_text(aes(x = c+5, y=theta_lft+0.5*theta), 
            label = TeX(paste0("\\hat{\\theta}=",number(theta,0.0001))),
            color= "orange", size=4)


# wrapper function to avoid two .R file's global variable conflict
## when use `include_graphics()` after a new global env
## see [url](https://www.r-bloggers.com/2020/08/why-i-dont-use-r-markdowns-ref-label/)

draw_band <- function(p_base=p_mxh, df=tb_mxh,
                      theta=theta,theta_lft = theta_lft,
                      theta_rgt = theta_rgt){
  p <- p_base +
  geom_ribbon(aes(x=xg,ymin = lwr, ymax = upr),
              data = base::subset(df,group=="control"),
              alpha = 0.2) +
  geom_ribbon(aes(x=xg,ymin = lwr, ymax = upr),
              data = base::subset(df,group=="treat"),
              alpha = 0.2) +
  geom_text(aes(x = c-3, y = theta_lft, 
                label=number(theta_lft, 0.0001)),
            color = "purple")+
  geom_text(aes(x = c+3, y = theta_rgt, 
                label=number(theta_rgt, 0.0001)),
            color = "blue") +
  geom_segment(aes(x=c, xend = c,
                   y= theta_lft, yend = theta_rgt),
               arrow = arrow(length = unit(0.1,"cm"),
                             type = "closed"),
               lty = "solid", color= "orange", lwd=lwd)  +
  geom_text(aes(x = c+5, y=theta_lft+0.5*theta), 
            label = TeX(paste0("\\hat{\\theta}=",number(theta,0.0001))),
            color= "orange", size=4) 
  return(p)
}

theta_lft1 <- theta_lft
theta_rgt1 <- theta_rgt
theta1 <- theta

p_band_wrapper <- draw_band(p_base=p_mxh, df=tb_mxh,
                            theta=theta1,theta_lft = theta_lft1,
                            theta_rgt = theta_rgt1)


#### RDD equivalent simple OLS regression ####

h0 <- h*sqrt(3)
w <- abs(x-c)<=h
y0 <- y[w]
x0 <- x[w]
T0 <- T[w]
n0 <- length(y0)
Z <- cbind(matrix(1,n0,1),x0,(x0-c)*T0,T0)
ZZ <- solve(crossprod(Z))
beta <- ZZ%*%crossprod(Z,y0)
e <- y0 - Z%*%beta
u <- Z*(e%*%matrix(1,1,4))
V <- ZZ%*%crossprod(u)%*%ZZ
s <- sqrt(diag(V))




