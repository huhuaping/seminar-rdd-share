#  This file generates Figure 19.1a
#  Nadaraya-Watson Regression

#### prepare ####
n <- 100
xm <- 10

# generate X and Y data
set.seed(180)
xdat <- runif(n,0,xm)
a <- pi/4
m <- sin((xdat-2)*a)/((xdat-2)*a)
ydat <- m + rnorm(n)/4

dt <- tibble(index = 1:length(xdat), X = xdat, Y = ydat)



#### calculate ####

h1 <- 1          # for binned
h2 <- h1/sqrt(3) # for NW

x0 <- seq(1,9,2) # bins break
x0n <- length(x0)
# index of Logic value in each bins
k0 <- (abs(x0%*%matrix(1,1,n)-matrix(1,x0n,1)%*%t(xdat)) <= h1)
# binned means
m0 <- (k0%*%ydat)/rowSums(k0)

# create new x with n=1001
x <- seq(0,xm,.01)
xn <- length(x)      # 1001

# for Rolling, dim(k1)=c(1001, 100)
k1 <- (abs(x%*%matrix(1,1,n)-matrix(1,xn,1)%*%t(xdat)) <= h1)
m1 <- (k1%*%ydat)/rowSums(k1)

# for NW, dim(k2)=c(1001, 100)
k2 <- dnorm(abs(x%*%matrix(1,1,n)-matrix(1,xn,1)%*%t(xdat))/h2)
m2 <- (k2%*%ydat)/rowSums(k2)

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
         bins = str_c('[',lft,',',rgt, ')')) %>%
  mutate(m0 = rep(m0, times=c(200,200,200,200,201)),
         m1 = as.vector(m1),
         m2 = as.vector(m2))


tbl_match <- tbl_result %>%
  mutate(data = map(x,~dt)) %>%
  unnest(data) %>%
  mutate(
    isbins = ifelse(
      ((X>=x-h1)& (X<(x+h1))),
      1, 0)
  )  %>%
  mutate(
    k0=ifelse(
      ((X>=lwr)& (X<upr)),
      1, 0),
    k1=as.numeric(as.vector(t(k1))),
    k2 = as.vector(t(k2))
    ) %>%
  mutate(bd = as.character(bd))


# binned mean
tbl_m0 <- tbl_match %>%
  filter(k0==1) %>%
  select(index, X, Y, bd,lwr, upr,mid,k0,m0) %>%
  unique() %>%
  arrange(index) %>%
  group_by(bd) %>%
  mutate(sum_k = sum(k0),
         sum_ky = sum(k0*Y)) %>%
  ungroup() %>%
  mutate(bins = bd) %>%
  rename("x"="mid")

# scrolling mean
tbl_m1 <- tbl_match %>%
  filter(k1==1) %>%
  select(index, X, Y, bd,bins,lft, rgt,x,k1,m1) %>%
  unique() %>%
  arrange(lft,index) %>%
  group_by(x) %>%
  mutate(sum_k = sum(k1),
         sum_ky = sum(k1*Y)) %>%
  ungroup() %>%
  rename("lwr"="lft", "upr"="rgt")

# kernel mean
tbl_m2 <- tbl_match %>%
  #filter(k1==1) %>%
  select(index, X, Y, bd,bins,lft, rgt,x,k2,m2) %>%
  unique() %>%
  arrange(lft,index) %>%
  group_by(x) %>%
  mutate(sum_k = sum(k2),
         sum_ky = sum(k2*Y)) %>%
  ungroup() %>%
  rename("lwr"="lft", "upr"="rgt")

#### extract result for different method ####

# help function
gen_bins <- function(df){
  out <- df %>% 
    select(#index, X, Y,
      x, bins,sum_ky, sum_k,
      starts_with("m")) %>%
    unique() %>%
    arrange(x) %>%
    mutate(x=number(x, 0.01))
}

bin1 <- gen_bins(tbl_m0)
n1 <- nrow(bin1)

bin2 <- gen_bins(tbl_m1)
n2 <- nrow(bin2)

bin3 <- gen_bins(tbl_m2)
n3 <- nrow(bin3)


#### draw plot ####

# basic scatter plot
fsize <- 16
p0 <- ggplot() +
  geom_point(aes(X, Y),data = dt, pch=21) +
  labs(x= "自变量X", y ="因变量Y") +
  scale_x_continuous(breaks = seq(0,10,1), limits = c(0,10)) +
  scale_y_continuous(breaks = seq(-1,2,0.5), limits = c(-1,2)) +
  theme_bw() +
  theme(text = element_text(size = 16))

### basic with bins
p00 <- p0 +
  geom_vline(xintercept = x0, lty="dashed") +
  geom_point(aes(X, Y, color=as.factor(bd)),
             data = tbl_m0, pch=21) +
  geom_rect(aes(ymin=-1,ymax=2,
                xmin= lwr, xmax = upr, 
                fill = as.factor(bd)),
            data = tbl_m0 %>% select(bd, lwr, upr) %>% unique,
            alpha = 0.05, inherit.aes = FALSE) +
  theme(legend.position = "none",
        text = element_text(size = 16))

### basic with means
p000 <- p00 +
  geom_point(aes(x, m0),
             data = tbl_m0 %>% select(x, m0) %>% unique(), 
             pch=15, color="black", size=3)

### binned plot
p1 <- p00 +
  geom_point(aes(x, m0),
             data = tbl_m0 %>% select(x, m0) %>% unique(), 
             pch=15, color="black", size=3) +
  geom_segment(aes(x =lwr, xend=upr,
                   y =m0, yend = m0,
                   color=as.factor(bd)),
            data = tbl_m0 %>% select(lwr,bd,upr, m0) %>% unique(),
            lty="solid", lwd=0.8) +
  theme(legend.position = "none",
        text = element_text(size = fsize))

## scrolling  plot

p20 <- p0 +
  geom_point(aes(x, m1, color=as.factor(bd)),
             data = tbl_m1 %>% select(x,m1,bins,bd) %>% unique(),
             pch=20) +
  theme(legend.position = "none",
        text = element_text(size = fsize))

p2 <- p0 +
  geom_point(aes(x, m0),
             data = tbl_m0 %>% select(x, m0) %>% unique(), 
             pch=15, color="black") +
  geom_line(aes(x, m1, color=as.factor(bd)),
            data = tbl_m1 %>% select(x,m1,bins,bd) %>% unique(), 
            lty="solid", lwd=0.8) +
  theme(legend.position = "none",
        text = element_text(size = fsize))

## Nadaraya-Watson  plot
p3 <- p0 +
  geom_point(aes(x, m0),
             data = tbl_m0 %>% select(x, m0) %>% unique(), 
             pch=15, color="black") +
  geom_line(aes(x, m2, color=as.factor(bd)),
            data = tbl_m2 %>% select(x,m2,bins,bd) %>% unique(), 
            lty="solid", lwd=0.8) +
  theme(legend.position = "none",
        text = element_text(size = fsize))

## all three plot

p_all <- p0 +
  geom_vline(xintercept = x0, lty="dashed")  +
  geom_point(aes(x, m0),
             data = tbl_m0 %>% select(x, m0) %>% unique(), 
             pch=15, size =3, color="black")  +
  geom_segment(aes(x =lwr, xend=upr,
                   y =m0, yend = m0,
                   color= "binned"),
               data = tbl_m0 %>% select(lwr,bd,upr, m0) %>% unique(),
               lty="solid", lwd=0.8)  +
  geom_line(aes(x, m1, color="rolling"),
            data = tbl_m1 %>% select(x,m1) %>% unique(), 
            lty="solid", lwd=0.8) +
  geom_line(aes(x, m2, color = "NW"),
            data = tbl_m2 %>% select(x,m2) %>% unique(), 
            lty="solid", lwd=0.8) +
  scale_color_manual(
    name="mx",
    breaks = c("binned", "rolling","NW"),
    values=c("green", "red","blue"))+
  theme(legend.position = "right",
        text = element_text(size = fsize))
