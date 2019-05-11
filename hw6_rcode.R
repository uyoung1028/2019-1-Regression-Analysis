library(ggplot2)
library(MASS)
#install.packages('lmtest')
library(lmtest)
library(dplyr)

getwd()
setwd('/Users/uyoung/Documents/Ewha/Graduate/data')

############### PART 1 ###############
## a)
data <- read.csv("Haskins-Default-100.csv",header=T)
bankruptcy <- data[,-c(1,3)]

f.m <- glm(Y~., data=bankruptcy, family=binomial(link=logit))
#r.m <- stepAIC(f.m, direction="backward")

AICHp <- data.frame(p=24:1,AIC=0)
variable<-list()
bankruptcy2 <- bankruptcy

for(i in 1:24){
  variable[[i]]<-colnames(bankruptcy2)
  H <- glm(Y~., data=bankruptcy2, family=binomial(link=logit))
  AICHp[i,2] <- H$aic
  k <- which.max(abs(summary(H)$coef[-1,4]))
  bankruptcy2 <- bankruptcy2[,-(k[[1]]+1)]
}

#data <- read.csv("Haskins-Default-100.csv",header=T)
#bankruptcy <- data[,-c(1,3)]

AICHp[25,2]<-glm(Y~1,data=bankruptcy2,family=binomial(link=logit))$aic
AICHp
variable[[which.min(AICHp$AIC)]]

dev.off()

ggplot(AICHp, aes(p,AIC)) + 
  geom_point() + 
  geom_line() 

## b)
r.m <- glm(Y~R5+R6+R16+R18+R22,data=bankruptcy,family=binomial(link=logit))
summary(r.m)

## c)
lrtest(r.m, f.m)

## d)
pchs <- nrow(bankruptcy)
pchs[] <- 1

pchs[bankruptcy$Y == 1] <- 1
pchs[bankruptcy$Y == 0] <- 4
pairs(Y~R5+R6+R16+R18+R22, bankruptcy, pch=pchs)

#data <- read.csv("Haskins-Default-100.csv",header=T)
#bankruptcy <- data[,-c(1,3)]

myglm <- summary(glm(Y~R5+R6+R16+R18+R22, data=bankruptcy, family=binomial(link=logit)))$coefficients
designX <- cbind(rep(1, 100), bankruptcy[,c(6,7,17,19,23)])
betaX <- as.matrix(designX)%*%as.matrix(summary(r.m)$coefficients[,1])
prob <- 1/(1+exp(-betaX))
bankruptcy3 <- data.frame(y=bankruptcy[,1],logit=log(prob/(1-prob)))

ggplot(bankruptcy3,aes(x=as.factor(y),y=logit,color=as.factor(y)))+geom_boxplot()
ggplot(bankruptcy3,aes(logit,group=as.factor(y),colour=as.factor(y)))+geom_density()


## e)
#data <- read.csv("Haskins-Default-100.csv",header=T)
#bankruptcy <- data[,-c(1,3)]

loglik.logis <- function(beta){
  z <- beta[1] + beta[2]*bankruptcy$R5 + beta[3]*bankruptcy$R6 + beta[4]*bankruptcy$R16 + beta[5]*bankruptcy$R18 + beta[6]*bankruptcy$R22
  p <- 1 / (1+exp(-z))
  -sum(log(p^bankruptcy$Y*(1-p)^(1-bankruptcy$Y)))
}
myglm <- summary(r.m)$coefficients

result <- optim(par=myglm[,1], fn=loglik.logis, method="BFGS", hessian=T)

H <- result$hessian
U <- solve(H)
U

mle <- c(H[1,1], H[1,2], H[1,3], H[1,4], H[1,5], H[1,6])
mle

se1 <- sqrt(diag(U))

ci.h <- function(par,hessian){
  se <- sqrt(diag(U))
  ci.U <- 0; ci.L <- 0
  for(i in 1:length(par)){
    ci.U[i] <- par[i] + 1.645*se[i]
    ci.L[i] <- par[i] - 1.645*se[i]
  }
  ci <- data.frame(par.hat=par, ci.L=ci.L, ci.U=ci.U)
  ci
}
ci.h(result$par,result$hessian)


## f)
prediction <- ifelse(prob >= 0.5, 1, 0)
err.table <- table(bankruptcy$Y, prediction)
err <- (err.table[1,2] + err.table[2,1]) / 100
err

fit2 <- glm(Y~.,data=bankruptcy, family=binomial(link=logit))
prob2 <- predict(fit2, type="response")
prediction2 <- ifelse(prob2 >= 0.5, 1, 0)
err.table2 <- table(bankruptcy$Y, prediction2)
err2 <- (err.table2[1,2] + err.table2[2,1]) / 100
err2


############### PART 2 ###############
## a)
profile0<-function(a) {
  loglike<-function(b) {
    z<-a+b[1]*bankruptcy$R5+b[2]*bankruptcy$R6+b[3]*bankruptcy$R16+b[4]*bankruptcy$R18+b[5]*bankruptcy$R22
    Gz<-1/(1+exp(-z))
    return(sum(log(Gz^bankruptcy$Y*(1-Gz)^(1-bankruptcy$Y))))
  }
  return(optim(par=rep(mle[-1]),loglike,hessian=T,control=list(fnscale=-1))$value)
}

profile1<-function(a) {
  loglike<-function(b) {
    z<-b[1]+a*bankruptcy$R5+b[2]*bankruptcy$R6+b[3]*bankruptcy$R16+b[4]*bankruptcy$R18+b[5]*bankruptcy$R22
    Gz<-1/(1+exp(-z))
    return(sum(log(Gz^bankruptcy$Y*(1-Gz)^(1-bankruptcy$Y))))
  }
  return(optim(par=mle[-2],loglike,hessian=T,control=list(fnscale=-1))$value)
}

profile2<-function(a) {
  loglike<-function(b) {
    z<-b[1]+b[2]*bankruptcy$R5+a*bankruptcy$R6+b[3]*bankruptcy$R16+b[4]*bankruptcy$R18+b[5]*bankruptcy$R22
    Gz<-1/(1+exp(-z))
    return(sum(log(Gz^bankruptcy$Y*(1-Gz)^(1-bankruptcy$Y))))
  }
  return(optim(par=mle[-3],loglike,hessian=T,control=list(fnscale=-1))$value)
}

profile3<-function(a) {
  loglike<-function(b) {
    z<-b[1]+b[2]*bankruptcy$R5+b[3]*bankruptcy$R6+a*bankruptcy$R16+b[4]*bankruptcy$R18+b[5]*bankruptcy$R22
    Gz<-1/(1+exp(-z))
    return(sum(log(Gz^bankruptcy$Y*(1-Gz)^(1-bankruptcy$Y))))
  }
  return(optim(par=mle[-4],loglike,hessian=T,control=list(fnscale=-1))$value)
}

profile4<-function(a) {
  loglike<-function(b) {
    z<-b[1]+b[2]*bankruptcy$R5+b[3]*bankruptcy$R6+b[4]*bankruptcy$R16+a*bankruptcy$R18+b[5]*bankruptcy$R22
    Gz<-1/(1+exp(-z))
    return(sum(log(Gz^bankruptcy$Y*(1-Gz)^(1-bankruptcy$Y))))
  }
  return(optim(par=mle[-5],loglike,hessian=T,control=list(fnscale=-1))$value)
}

profile5<-function(a) {
  loglike<-function(b) {
    z<-b[1]+b[2]*bankruptcy$R5+b[3]*bankruptcy$R6+b[4]*bankruptcy$R16+b[5]*bankruptcy$R18+a*bankruptcy$R22
    Gz<-1/(1+exp(-z))
    return(sum(log(Gz^bankruptcy$Y*(1-Gz)^(1-bankruptcy$Y))))
  }
  return(optim(par=mle[-6],loglike,hessian=T,control=list(fnscale=-1))$value)
}


b0<-seq(-5,5,length.out=100)
prof.loglike0<-c()
for (i in 1:length(b0)) {
  prof.loglike0[i]<-profile0(b0[i])
}

b1<-seq(10,50,length.out=100)
prof.loglike1<-c()
for (i in 1:length(b1)) {
  prof.loglike1[i]<-profile1(b1[i])
}

b2<-seq(-50,-10,length.out=100)
prof.loglike2<-c()
for (i in 1:length(b2)) {
  prof.loglike2[i]<-profile2(b2[i])
}

b3<-seq(-70,0,length.out=100)
prof.loglike3<-c()
for (i in 1:length(b3)) {
  prof.loglike3[i]<-profile3(b3[i])
}

b4<-seq(15,55,length.out=100)
prof.loglike4<-c()
for (i in 1:length(b4)) {
  prof.loglike4[i]<-profile4(b4[i])
}

b5<-seq(-15,45,length.out=100)
prof.loglike5<-c()
for (i in 1:length(b5)) {
  prof.loglike5[i]<-profile5(b5[i])
}

ggplot(data.frame(b0,prof.loglike0),aes(b0,prof.loglike0))+geom_point()+geom_line()+ggtitle("Profile likelihood - beta0")
ggplot(data.frame(b1,prof.loglike1),aes(b1,prof.loglike1))+geom_point()+geom_line()+ggtitle("Profile likelihood - beta1")
ggplot(data.frame(b2,prof.loglike2),aes(b2,prof.loglike2))+geom_point()+geom_line()+ggtitle("Profile likelihood - beta2")
ggplot(data.frame(b3,prof.loglike3),aes(b3,prof.loglike3))+geom_point()+geom_line()+ggtitle("Profile likelihood - beta3")
ggplot(data.frame(b4,prof.loglike4),aes(b4,prof.loglike4))+geom_point()+geom_line()+ggtitle("Profile likelihood - beta4")
ggplot(data.frame(b5,prof.loglike5),aes(b5,prof.loglike5))+geom_point()+geom_line()+ggtitle("Profile likelihood - beta5")

mydata<-data.frame(b0=b0,b1=b1,b2=b2,b3=b3,b4=b4,b5=b5,newy0=NA,newy1=NA,newy2=NA,newy3=NA,newy4=NA,newy5=NA)

for (i in 1:length(b0)) {
  mydata[i,7]<-sqrt(2*(profile0(mle[1])-profile0(mydata[i,1])))*sign(mydata[i,1]-mle[1])
  mydata[i,8]<-sqrt(2*(profile1(mle[2])-profile1(mydata[i,2])))*sign(mydata[i,2]-mle[2])
  mydata[i,9]<-sqrt(2*(profile2(mle[3])-profile2(mydata[i,3])))*sign(mydata[i,3]-mle[3])
  mydata[i,10]<-sqrt(2*(profile3(mle[4])-profile3(mydata[i,4])))*sign(mydata[i,4]-mle[4])
  mydata[i,11]<-sqrt(2*(profile4(mle[5])-profile4(mydata[i,5])))*sign(mydata[i,5]-mle[5])
  mydata[i,12]<-sqrt(2*(profile5(mle[6])-profile5(mydata[i,6])))*sign(mydata[i,6]-mle[6])
}

CIb0<-b0[c(nrow(as_tibble(mydata) %>% filter(newy0<=-1.645))+1,100-nrow(as_tibble(mydata) %>% filter(newy0>=1.645)))]
CIb1<-b1[c(nrow(as_tibble(mydata) %>% filter(newy1<=-1.645))+1,100-nrow(as_tibble(mydata) %>% filter(newy1>=1.645)))]
CIb2<-b2[c(nrow(as_tibble(mydata) %>% filter(newy2<=-1.645))+1,100-nrow(as_tibble(mydata) %>% filter(newy2>=1.645)))]
CIb3<-b3[c(nrow(as_tibble(mydata) %>% filter(newy3<=-1.645))+1,100-nrow(as_tibble(mydata) %>% filter(newy3>=1.645)))]
CIb4<-b4[c(nrow(as_tibble(mydata) %>% filter(newy4<=-1.645))+1,100-nrow(as_tibble(mydata) %>% filter(newy4>=1.645)))]
CIb5<-b5[c(nrow(as_tibble(mydata) %>% filter(newy5<=-1.645))+1,100-nrow(as_tibble(mydata) %>% filter(newy5>=1.645)))]

CIb0; CIb1; CIb2; CIb3; CIb4; CIb5;



