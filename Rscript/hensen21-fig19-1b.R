##  This file generates Figure 19.1b
##  Local Linear Regression

#### prepare ####

n <- 100
xm = 10
x <- seq(0,xm,.01)
xn <- length(x)
set.seed(180)
xdat <- runif(n,0,xm)
a <- pi/4
m <- sin(xdat*a)/(xdat*a)
m <- sin((xdat-2)*a)/((xdat-2)*a)
ydat <- m + rnorm(n)/4

dt <- tibble(index = 1:length(xdat), X = xdat, Y = ydat)

#### calculate ####

# set bandwidth
h1 <- 1
h2 <- h1/sqrt(3)

# Rolling LL and Kernel LL estimation
mg <- matrix(0,xn,1)
mr <- matrix(0,xn,1)
for (j in 1:xn){
  xj <- xdat-x[j]
  z <- cbind(matrix(1,n,1),xj)
  k1 <- abs(xj) < h1      # Rectangle Kernel
  k2 <- dnorm(xj/h2)      # Gaussian Kernel
  zk1 <- z*k1
  zk2 <- z*k2
  betar <- solve(t(zk1)%*%z,t(zk1)%*%ydat)
  betag <- solve(t(zk2)%*%z,t(zk2)%*%ydat)
  mr[j] <- betar[1]   # for Rolling LL
  mg[j] <- betag[1]   # for Kernel LL
}


# Binned LL estimation
x0 <- c(1,3,5,7,9)
## for bins 1
x1 <- as.matrix(subset(xdat,xdat<2))
y1 <- as.matrix(subset(ydat,xdat<2))
z1 <- cbind(matrix(1,length(y1),1),x1)
beta1 <- solve(t(z1)%*%z1,t(z1)%*%y1)   # get OLS  coefficients
m1 <- beta1[1]+beta1[2]*x0[1]    # mean centers
x1 <- as.matrix(subset(x,x<2))
f1 <- beta1[1]+x1*beta1[2]     # forecast with new data

## for bins 2
x2 <- as.matrix(subset(xdat,(xdat>=2)&(xdat<4)))
y2 <- as.matrix(subset(ydat,(xdat>=2)&(xdat<4)))
z2 <- cbind(matrix(1,length(y2),1),x2)
beta2 <- solve(t(z2)%*%z2,t(z2)%*%y2)
m2 <- beta2[1]+beta2[2]*x0[2]          
x2 <- as.matrix(subset(x,(x>=2)&(x<4)))
f2 <- beta2[1]+x2*beta2[2]

## for bins 3
x3 <- as.matrix(subset(xdat,(xdat>=4)&(xdat<6)))
y3 <- as.matrix(subset(ydat,(xdat>=4)&(xdat<6)))
z3 <- cbind(matrix(1,length(y3),1),x3)
beta3 <- solve(t(z3)%*%z3,t(z3)%*%y3)
m3 <- beta3[1]+beta3[2]*x0[3]
x3 <- as.matrix(subset(x,(x>=4)&(x<6)))
f3 <- beta3[1]+x3*beta3[2]

## for bins 4
x4 <- as.matrix(subset(xdat,(xdat>=6)&(xdat<8)))
y4 <- as.matrix(subset(ydat,(xdat>=6)&(xdat<8)))
z4 <- cbind(matrix(1,length(y4),1),x4)
beta4 <- solve(t(z4)%*%z4,t(z4)%*%y4)
m4 <- beta4[1]+beta4[2]*x0[4]
x4 <- as.matrix(subset(x,(x>=6)&(x<8)))
f4 <- beta4[1]+x4*beta4[2]

## for bins 5
x5 <- as.matrix(subset(xdat,xdat>=8))
y5 <- as.matrix(subset(ydat,xdat>=8))
z5 <- cbind(matrix(1,length(y5),1),x5)
beta5 <- solve(t(z5)%*%z5,t(z5)%*%y5)
m5 <- beta5[1]+beta5[2]*x0[5]
x5 <- as.matrix(subset(x,(x>=8)))
f5 <- beta5[1]+x5*beta5[2]

## combine as data.frame
tbl_a0 <- cbind(beta1, beta2,beta3, beta4,beta5) %>%
  as_tibble(.name_repair = "unique") %>%
  rename_all(~paste0("bd",1:5)) %>%
  add_column(par= c("intercept", "slope"), .before = "bd1")

##  binned LL CEF estimation
m0 <- c(m1,m2,m3,m4,m5)        ## mean centers
f0 <- c(f1, f2,f3, f4, f5)     ## LL estimation

#### combine all result ####

tbl_result <- tibble(id = 1:length(x),
                     x = x,
                     x_upr = x +0.01) %>%
  mutate(bd = cut_interval(x,length = 2, right=FALSE)) %>%
  mutate(
    lwr = as.numeric(str_extract(bd, "(\\d{1})(?=\\,)")),
    upr = as.numeric(str_extract(bd, "(?<=\\,)(\\d{1,2})")),
    mid = 0.5*(lwr+ upr)
  ) %>%
  mutate(lft = x-h1,
         rgt = x +h1,
         bins = str_c('[',number(lft,0.01),',',
                      number(rgt,0.01), ')')) %>%
  mutate(my = rep(m0, times=c(200,200,200,200,201)),
         m0 = (as.vector(f0)),
         m1 = (as.vector(mr)),
         m2 = as.vector(mg))


tbl_match <- tbl_result %>%
  mutate(data = map(x,~dt)) %>%
  unnest(data) %>%
  mutate(
    isbins = ifelse(
      ((X>=x-h1)& (X<(x+h1))),
      1, 0)
  )  %>%
  mutate(bd = as.character(bd)) %>%
  mutate(
    isbd = ifelse(
      ((X>=lwr)& (X<(upr))),
      1, 0)
  )


# binned LL
tbl_m0 <- tbl_match %>%
  filter(isbins==1) %>%
  select(index, X, Y, bd,bins,lwr, upr,mid, x,my,m0, isbins) %>%
  unique() %>%
  arrange(lwr,index) %>%
  group_by(x) %>%
  mutate(sum_k = sum(isbins)) %>%
  ungroup() 


# rolling LL
tbl_m1 <- tbl_match %>%
  filter(isbins==1) %>%
  select(index, X, Y, bd,bins,lft, rgt,x,m1, isbins) %>%
  unique() %>%
  arrange(lft,index) %>%
  group_by(x) %>%
  mutate(sum_k = sum(isbins)) %>%
  ungroup() %>%
  rename("lwr"="lft", "upr"="rgt")

# kernel LL
tbl_m2 <- tbl_match %>%
  filter(isbins==1) %>%
  select(index, X, Y, bd,bins,lft, rgt,x,m2, isbins) %>%
  unique() %>%
  arrange(lft,index) %>%
  group_by(x) %>%
  mutate(sum_k = sum(isbins)) %>%
  ungroup() %>%
  rename("lwr"="lft", "upr"="rgt")

#### draw plot ####

# basic plot
fsize <- 16
p0 <- ggplot() +
  geom_point(aes(X, Y),data = dt, pch=21) +
  labs(x= "自变量X", y ="因变量Y") +
  scale_x_continuous(breaks = seq(0,10,1), limits = c(0,10)) +
  scale_y_continuous(breaks = seq(-1,2,0.5), limits = c(-1,2)) +
  theme_bw() +
  theme(text = element_text(size = fsize))

p00 <- p0 +
  geom_vline(xintercept = x0, lty="dashed") +
  geom_point(aes(X, Y, color=as.factor(bd)),
             data = tbl_m0 %>% select(X, Y, bd) %>% unique(),
             pch=21) +
  geom_rect(aes(ymin=-1,ymax=2,
                xmin= lwr, xmax = upr, 
                fill = as.factor(bd)),
            data = tbl_m0 %>% select(bd, lwr, upr) %>% unique(),
            alpha = 0.05, inherit.aes = FALSE) +
  theme(legend.position = "none",
        text = element_text(size = fsize))

### binned LL plot
p1 <- p00 +
  geom_point(aes(mid, my),
             data = tbl_m0 %>% select(mid, my) %>% unique(), 
             pch=15, color="black") +
  geom_line(aes(x =x, y = m0,
                color=as.factor(bd)),
               data = tbl_m0 %>% select(x,m0, bd) %>% unique(),
               lty="solid", lwd=0.8) +
  theme(legend.position = "none",
        text = element_text(size = fsize))

## rolling LL plot

p2 <- p0 +
  geom_point(aes(mid, my),
             data = tbl_m0 %>% select(mid, my) %>% unique(), 
             pch=15, color="black")  +
  geom_line(aes(x =x, y = m1,
                color=as.factor(bd)),
            data = tbl_m1 %>% select(x,m1, bd) %>% unique(),
            lty="solid", lwd=0.8) +
  theme(legend.position = "none",
        text = element_text(size = fsize))

## Kernel LL  plot
p3 <- p0 +
  geom_point(aes(mid, my),
             data = tbl_m0 %>% select(mid, my) %>% unique(), 
             pch=15, color="black")  +
  geom_line(aes(x =x, y = m2,
                color=as.factor(bd)),
            data = tbl_m2 %>% select(x,m2, bd) %>% unique(),
            lty="solid", lwd=0.8)+
  theme(legend.position = "none",
        text = element_text(size = fsize))

## all three plot
p_all <- p0 +
  geom_vline(xintercept = x0, lty="dashed")  +
  geom_point(aes(mid, my),
             data = tbl_m0 %>% select(mid, my) %>% unique(), 
             pch=15, color="black", size =2)  +
  geom_line(aes(x =x, y = m0, 
                color="binned reg"),
            data = tbl_m0 %>% select(x,m0, bd) %>% unique(),
            lty="solid", lwd=0.8) +
  geom_line(aes(x, m1, color="rolling reg"),
            data = tbl_m1 %>% select(x,m1, bd) %>% unique(),
            lty="solid", lwd=0.8) +
  geom_line(aes(x, m2, color = "LL reg"),
            data = tbl_m2 %>% select(x,m2, bd) %>% unique(),
            lty="solid", lwd=0.8) +
  scale_color_manual(
    name="mx reg",
    breaks = c("binned reg", "rolling reg","LL reg"),
    values=c("green", "red","blue"))+
  theme(legend.position = "right",
        text = element_text(size = fsize))

