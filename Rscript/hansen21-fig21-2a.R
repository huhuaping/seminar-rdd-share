##  This file generates Figure 21.2a and Table 21.1 Covariates
##  Head Start RDD with Covariates
##  This file uses the package haven
##  This file uses the data file LM2007.dta

####  The data file LM2007.dta is used ####
library(haven)

dt_path <- here::here("data/LM2007.dta")
dat <- read_dta(dt_path)
y <- dat$mort_age59_related_postHS
x <- dat$povrate60
Za <- dat$census1960_pctblack
Zb <- dat$census1960_pcturban

c <- 59.1984

dt_head <- tibble(X= x,
                  Y = y,
                  Za = Za,
                  Zb = Zb) %>%
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

# specify bins
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
  e <- LL_Residual(y,x,h)
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

#### step 1: LL estimate Y on X ####
## obtain residual of stage 1 
e1 <- LL_Residual(y1,x1,h)
e2 <- LL_Residual(y2,x2,h)
e <- rbind(e1,e2)

#### step 2: LL estimate Z on X ####
## obtain residual of stage 2 
Ra <- LL_Residual(Za,x,h)
Rb <- LL_Residual(Zb,x,h)
Z <- cbind(Za,Zb)
R <- cbind(Ra,Rb)

# step 3: OLS  ####
## regress estimate residual of stage 1 on residual of stage 2
## obtain the coefﬁcient and standard errors
## no intercept !
XXR <- solve(crossprod(R))     # inverse matrix
beta <- XXR%*%(t(R)%*%e)
u <- e - R%*%beta
Ru <- R*(u%*%matrix(1,1,ncol(R)))
V <- XXR%*%crossprod(Ru)%*%XXR    # robust se
sbeta <- sqrt(diag(V))

tstat <- beta/sbeta
pvalue <- 2*(1-pnorm(abs(tstat)))
betas <- cbind(beta,sbeta,tstat,pvalue)

#### step 4: Construct the residual ####
## obtain constructed residual

zm <- colMeans(Z)%*%beta

Z1 <- Z[T==0,]
Z2 <- Z[T==1,]
yZ <- y - Z%*%beta
yZ1 <- y1 - Z1%*%beta
yZ2 <- y2 - Z2%*%beta

#### step 5: LL ATE estimate ####
## LL estimate of constructed residual on X
## obtain  ATE and standard errors

m1 <- LL_EST(yZ1,x1,g1,h) + matrix(zm,G1,1)
m2 <- LL_EST(yZ2,x2,g2,h) + matrix(zm,G2,1)
s1 <- LL_SE(yZ1,x1,g1,h)
s2 <- LL_SE(yZ2,x2,g2,h)
L1 <- m1-1.96*s1
U1 <- m1+1.96*s1
L2 <- m2-1.96*s2
U2 <- m2+1.96*s2

# rdd effect
theta_lft <- m1[G1]
theta_rgt <- m2[1]
theta <- m2[1]-m1[G1]
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


#### tibble of sample result ####
tbl_residual <- dt_head %>%
  mutate(e = as.vector(e),   # residual of LL 1
         Ra = as.vector(Ra),   # residual of LL 2
         Rb = as.vector(Rb),
         RZ = as.vector(yZ) # residual of construct
  )

#### tibble of estimate result ####
tbl_theta02 <- tibble(
  model = "covariate",
  pars = c("theta", "black","urban"),
  est = c(theta, as.vector(beta)),
  se  = c(setheta, as.vector(sbeta))
)

#### plot CEF mx ####

lwd <- 0.6
lwadd <- 0.2

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

# mx plot for both side
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
  theme(legend.position = "right",
        text = element_text(size = fsize))

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
p_band_cov <- p_mxh +
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


