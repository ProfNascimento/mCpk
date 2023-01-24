#' Plot Cpk
#'
#' This function plot the modified interval Cpk.
#'
#' @import ggplot2
#' @importFrom MASS fitdistr
#' @param data,n,m,LSE,LIE,dist,B0,CI Path to the input file
#' @return tabela
#' @export
cpk.plot = function(data,n,m,LSE,LIE,dist,B0=1000,CI=0.95){
  Cpk<-c();  CpkSup<-c();
  CpkInf<-c();  CpkB<-c()
  SAMPLE = 1:m;             #PROCESS WINDOW SIZE
  B = B0                    #BOOTSTRAP SAMPLE SIZE

  #DATA WRANGLING
  data <- qcc::qcc.groups(as.numeric(dados[(length(dados)-n*m+1):length(dados)]),
                          rep(SAMPLE,n))

  # NORMAL DISTRIBUTION
  if(dist == "Normal"){
    for(i in 1:m){
      w<-c()

      for(j in 1:n){
        w[j]<-data[i,j]
      }

      nfit = fitdistr(w,"normal")
      media = nfit$estimate[1]
      dp = nfit$estimate[2]
      cpknorm = min((LSE - media)/(3*dp),(media-LIE)/(3*dp))
      Cpk[i]<-cpknorm

      for(k in 1:B){
        l = rnorm(n,media ,dp)
        nfitB = fitdistr(l,"normal")
        mediaB = nfitB$estimate[1]
        dpB = nfitB$estimate[2]
        CpkB[k]= min(abs(LSE - mediaB)/(3*dpB),(mediaB-LIE)/(3*dpB))
      }

      CpkSup[i]<-quantile(CpkB, probs = c(CI+(1-CI)/2), na.rm = FALSE,names = FALSE,type = 7)
      CpkInf[i]<- quantile(CpkB, probs = c((1-CI)/2), na.rm = FALSE,names = FALSE,type = 7)
    }
  }
  # GAMMA DISTRIBUTION
  if (dist=="Gamma"){
    for(i in 1:m){
      x<-c()
      for(j in 1:n){
        x[j]<-data[i,j]
      }

      s2<-sum(x)
      halpha<-n/(((n*sum((x)*log(x))/s2)-sum(log(x))))
      hbet<-(1/(n*(n-1)))*(n*sum((x)*log(x))-s2*sum(log(x)))

      Upv <- qgamma(c(CI+(1-CI)/2),shape=halpha,scale=hbet)
      Lpv <- qgamma(c((1-CI)/2),shape=halpha,scale=hbet)
      Mv <-  qgamma(0.5,shape=halpha,scale=hbet)
      cpkgama <- min((LSE - Mv)/(Upv-Mv), (Mv - LIE)/(Mv-Lpv))
      Cpk[i]<-cpkgama

      for(k in 1:B){
        z = rgamma(n,shape=halpha ,scale=hbet) ;

        balphaM = (n*sum(z)) / (n*sum(z*log(z)) - (sum(z)*sum(log(z))) )
        bbetM =   (1/(n^2))*(n*sum((z)*log(z))- (sum(z)*sum(log(z))) )

        UM <- qgamma(c(CI+(1-CI)/2),balphaM,scale=bbetM)
        LM <- qgamma(c((1-CI)/2),balphaM,scale=bbetM)
        MM <- qgamma(0.5,shape=balphaM,scale=bbetM)
        cpkbootM<- min((LSE - MM)/(UM-MM), (MM - LIE)/(MM-LM))
        CpkB[k] <- cpkbootM
      }

      CpkSup[i]<-quantile(CpkB, probs = c(CI+(1-CI)/2), na.rm = FALSE,names = FALSE,type = 7)
      CpkInf[i]<- quantile(CpkB, probs = c((1-CI)/2), na.rm = FALSE,names = FALSE,type = 7)
    }
  }
  # WEIBULL DISTRIBUTION
  if (dist=="Weibull"){

    for(i in 1:m){
      y<-c()

      for(j in 1:n){
        y[j]<-data[i,j]
      }

      y = sort(y);pos=c(1:length(y))

      m1=mean(y)
      m2 = try(((2/(n^2 - n)) * (sum((pos - 1)*y  ))) - mean(y))
      formW =  try(( - log(2) / ( log( 1 -( m2 / m1) ) )))
      escW = try((m1 /gamma( 1 / formW+1)))

      UM <- qweibull(c(CI+(1-CI)/2),formW,scale=escW)
      LM <- qweibull(c((1-CI)/2),formW,scale=escW)
      MM <- qweibull(0.5,formW,scale=escW)
      cpkweibull<- min((LSE - MM)/(UM-MM), (MM - LIE)/(MM-LM))
      Cpk[i]<-cpkweibull

      for(k in 1:B){

        w = sort(rweibull(n,shape=formW ,scale=escW))
        m1y = mean(w)
        posy = c(1:length(w))
        m2y = try(((2/(n^2 - n)) * (sum((pos - 1)*w ))) - mean(w))
        balpha =  try(( - log(2) / ( log( 1 -( m2y / m1y) ) )))
        bbeta = try((m1y /gamma(1 / balpha+1)))

        UM <- qweibull(c(CI+(1-CI)/2),balpha,scale=bbeta)
        LM <- qweibull(c((1-CI)/2),balpha,scale=bbeta)
        MM <- qweibull(0.5,balpha,scale=bbeta)
        cpkbootM<- min((LSE - MM)/(UM-MM), (MM - LIE)/(MM-LM))
        CpkB[k] <- cpkbootM
      }

      CpkSup[i]<-quantile(CpkB, probs = c(CI+(1-CI)/2), na.rm = FALSE,names = FALSE,type = 7)
      CpkInf[i]<- quantile(CpkB, probs = c((1-CI)/2), na.rm = FALSE,names = FALSE,type = 7)

    }
  }

  if(dist == "Normal"|dist=="Gamma"|dist=="Weibull"){
    CPK = as.data.frame(cbind(SAMPLE,Cpk,CpkSup,CpkInf))
    inferior<-min(CpkInf) ; superior<-max(CpkSup)

    gr = ggplot(CPK, aes(x=SAMPLE,y=Cpk)) +
      scale_y_continuous(limits = c(min(Cpk,inferior),max(1.6,superior))) +
      geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1,fill = "tomato", alpha = 0.02) +
      geom_rect(xmin = -Inf, xmax = Inf, ymin = 1.33, ymax = Inf, fill = "cadetblue3", alpha = 0.02) +
      geom_rect(xmin = -Inf, xmax = Inf, ymin = 1, ymax = 1.33, fill = "azure4", alpha = 0.02) +
      geom_line(size=1) +  geom_point(size=3) +
      geom_line(aes(y = CpkInf), linetype=2,size=1,col="dimgray") +
      geom_line(aes(y = CpkSup), linetype=2,size=1,col="dimgray")

    op = print(gr)

    tab = CPK[,-1]
    tabela=cbind(tab)
    colnames(tabela)=c("Cpk point est.","UI","LI")

    print( paste("Process Cpk Upper-Interv.(UI): ", round(mean(CpkSup),2)) )
    print( paste("Process point Cpk: ", round(mean(Cpk),2)) )
    print("SUMMARY")
    print(
      matrix(
        c(sum(CPK[,3]<1),sum(CPK[,3]<1.33 & CPK[,3]>1),sum(CPK[,3]>1.33),
          sum(CPK[,2]<1),sum(CPK[,2]<1.33 & CPK[,2]>1),sum(CPK[,2]>1.33),
          sum(CPK[,4]<1),sum(CPK[,4]<1.33 & CPK[,4]>1),sum(CPK[,4]>1.33)), ncol=3,byrow = TRUE,
        dimnames = list(c("Cpk UI","Cpk", "Cpk LI"),
                        c("Not","Barely","Capable"))
      )
    )
    return(tabela)
    par(op)
  }
}
