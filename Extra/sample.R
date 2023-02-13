library(qcc)
library(MASS)
library(ggplot2)
library(gridExtra)

## DATA IMPORT
dados <- read.csv("https://raw.githubusercontent.com/ProfNascimento/mCpk/main/sample.csv", sep=",")

Achoc=qcc.groups(data=dados$weight,dados$sample)
cartaX= qcc(Achoc,type = "xbar")

# Shewhart rule
shewhart.rules(cartaX)
violating.runs(cartaX) #elemento de numero 10, removê-lo

AchocCorrigido = Achoc[-10,]

cartaXCorrigido= qcc(AchocCorrigido,type = "xbar")
cartaR = qcc(AchocCorrigido,type="R")

################
# Control Charts
################
n=5;m=29-1
Subamostra=1:m

medias= apply(AchocCorrigido, 1,mean)
ampli = apply(AchocCorrigido, 1,range)

amplitude=c()

for(j in 1:m){
  amplitude[j]<- (ampli[2,j]) - (ampli[1,j])
}
d2=2.3248;d3=0.8674
SD= mean(amplitude)/d2

######################
## Range Control Chart
######################
R= data.frame(cbind(amplitude,Subamostra))
Rbar=mean(amplitude)
LSCr = Rbar + 3*(d3*Rbar/(d2))
LICr = Rbar - 3*(d3*Rbar/d2) ; if(LICr<0){LICr=0}

ggplot(R, aes(x=Subamostra,y=amplitude))+
  geom_line(size=1,col="darkslategray4") +  geom_point(size=3,col="darkslategray4") +
  geom_line(aes(y = Rbar), linetype=1,size=1,col="firebrick")+
  geom_line(aes(y = LICr), linetype=2,size=1,col="firebrick")+
  geom_line(aes(y = LSCr), linetype=2,size=1,col="firebrick")+
  scale_y_continuous("Amplitude") + scale_x_continuous("Amostra")

#####################
## Mean Control Chart
#####################
n=5
Xbar= data.frame(cbind(medias,Subamostra))
LM= mean(medias)
LSC= LM + (3*SD)/sqrt(n)
LIC= LM - 3*SD/sqrt(n)

ggplot(Xbar, aes(x=Subamostra,y=medias))+
  geom_line(size=1,col="darkslategray4") +  geom_point(size=3,col="darkslategray4") +
  geom_line(aes(y = LM), linetype=1,size=1,col="firebrick")+
  geom_line(aes(y = LIC), linetype=2,size=1,col="firebrick")+
  geom_line(aes(y = LSC), linetype=2,size=1,col="firebrick")+
  scale_y_continuous("Média") + scale_x_continuous("Amostra")

#############################
##  Interval Cpk closed-form
#############################
N=length(AchocCorrigido)
x=sort(as.numeric(AchocCorrigido))

#MLE Normal
fit = fitdistr(x,"normal")
media =  fit$estimate[1]
desvio = fit$estimate[2]

#MMLE Gamma
escG= try((sum(x)/N))
formG = try((N*sum(x)) / (N*sum(x*log(x)) - (sum(x)*sum(log(x))) ))

#LMM Weibull
pos = (1:length(x))
m2 = try(((2/(length(x)^2 - length(x))) * (sum((pos - 1)*x  ))) - mean(x))
formW =  try(( - log(2) / ( log( 1 -( m2 / mean(x)) ) )))
escW = try((mean(x) /gamma( 1 / formW+1)))

ks.test(x,"pgamma",scale=(escG/formG),shape=formG)
ks.test(x,"pweibull",shape=formW,scale=escW)
ks.test(x,"pnorm",mean=fit$estimate[1],sd=fit$estimate[2])

## CHECKING DIST. ADJUSTMENT
ggplot(data.frame(x),aes(x)) +
  geom_histogram(aes(y = ..density..),bins=8,col="azure4",fill="white") +
  geom_line(aes(y = ..density.., colour = "Smooth"), size=1,stat = "density")+
  stat_function(fun = dweibull, size=1,args=list(shape=formW,scale=escW),aes(colour = 'Weibull'))+
  stat_function(fun = dgamma,size=1, args=list(shape=formG,scale=(escG/formG)),aes(colour = 'Gamma'))+
  stat_function(fun = dnorm, args=list(mean=media,sd=desvio),aes(colour = 'Normal'))+
  scale_colour_manual(name="PROB.", values=c("#9999CC","Azure3","lightpink2","Cadetblue"))+
  scale_y_continuous("Density") + scale_x_continuous("Weight (g)")

##########################
## modified Cpk (Weibull)
# LSE=412; LIE=388
##########################
# devtools::install_github("ProfNascimento/mCpk")
library(mCpk)

## RAW DATA
Test = cpk.plot(Achoc,n=5,m=29,LIE=388,LSE=412,dist="Weibull")
Test

Test2 = cpk.plot(Achoc,n=5,m=29,LIE=388,LSE=412,dist="Normal")

## CLEANED DATA
Test3 = cpk.plot(x,n=5,m=28,LIE=388,LSE=412,dist="Weibull")
Test3

##-----------------------------------------------------##
# ----------------Visualizing Test3-------------------- #
# Process Specification vs Empirical displacement (trend)
x=rweibull(10000,shape=formW,scale=escW)
p1 = ggplot(as.data.frame(x), aes(x=x)) + geom_density(alpha=.3) + xlab("Process Specification - Grams (g)")+
     geom_vline(aes(xintercept = 400), linetype=1,size=1,col="firebrick") +
     geom_vline(aes(xintercept = 388), linetype=2,size=1,col="firebrick") +
     geom_vline(aes(xintercept = 412), linetype=2,size=1,col="firebrick") + xlim(c(380,420))+ coord_flip()

p2 = ggplot(Xbar, aes(x=Subamostra,y=medias))+
  geom_line(size=1,col="darkslategray4") +  geom_point(size=3,col="darkslategray4") +
  geom_line(aes(y = LM), linetype=1,size=.7,col="firebrick")+
  geom_line(aes(y = LIC), linetype=2,size=.7,col="firebrick")+
  geom_line(aes(y = LSC), linetype=2,size=.7,col="firebrick")+
  ylab("Control Chart (x-bar)") + xlab("Sample (m)")+
  ylim(c(380,420)) 

cowplot::plot_grid(p1,p2,nrow=1, rel_widths = c(1, 3))
