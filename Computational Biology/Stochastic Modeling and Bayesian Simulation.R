#########################
# exemplul 1 random walk
Walk1d <- function() {
  n <- 100
  y <- vector(length=n)
  y[1] <- 0
  for (i in 2:n) y[i] <- y[(i-1)] + sample(c(-1,1),1) #random sample
  plot(1:n,y,type='l',ylim=c(-20,20))
  }
Walk1d()

par(mfrow=c(2,2))
Walk2d<-function(){
  xstart <- 0
  ystart <- 0
  xmove <- sample(c(-1,1),500,repl=T)
  ymove <- sample(c(-1,1),500,repl=T)
  xmove <- xstart + cumsum(xmove)
  ymove <- ystart + cumsum(ymove)
  plot(xmove,ymove,xlim=c(-40,40),ylim=c(-40,40),xlab='x',ylab='y',type='l')
 }
for (i in 1:4) Walk2d()

#########################
# exemplul 2 Modelarea unei secvente de ADN

 DNA<- matrix(c(0.3,0.1,0.2,0.1,0.2,0.2,0.2,0.8,0.2,0.4,0.2,0.1,0.3,0.3,0.4,0),
          nrow=4,ncol=4)
 DNA2=DNA%*%DNA
 DNA4<-DNA2%*%DNA2
 DNA8<-DNA4%*%DNA4
 DNA16<-DNA8%*%DNA8
 
#given our initial distribution the posterior distribution
#of nucleotides is 15% A, 35%T, 25%C and 25%G, regardless of position.
 
 
 ####################################
 # exemplul 3 Lant Markov ascuns
 states<- c("Begin", "Exon", "Donor", "Intron") 
 observations<- c("A", "C", "G", "T") 
 
 #transition matrix 
 Begin_p<- c(0,1,0,0) 
 Exon_p<- c(0,0.9,0.1,0) 
 Donor_p<-c(0,0,0,1)
 Intron_p<-c(0,0,0.05,0.95)
 transition_m<-matrix(c(Begin_p,Exon_p, Donor_p,Intron_p), 4, 4, byrow = TRUE) 
 rownames(transition_m) <- states
 colnames(transition_m) <- states
 transition_m
 
 # emission matrix
 Begin_states_p<- c(0,0,0,0) 
 Exon_states_p<- c(0.25,0.25,0.25,0.25) 
 Donor_states_p<-c(0.05,0,0.95,0)
 Intron_states_p<-c(0.4,0.1,0.1,0.4)
 
 emission_m <- matrix(c(Begin_states_p, Exon_states_p,
                        Donor_states_p,Intron_states_p), 4, 4, byrow = TRUE) # Create a 2 x 4 matrix
 rownames(emission_m) <- states
 colnames(emission_m) <- observations
 emission_m 
 
 library(HMM)
 hmm=initHMM(states,observations,startProbs=NULL,transition_m,
             emission_m)
 
 myseq <- c("C", "T", "T", "C", "A" , "T", "G", "T", "G", "A", "A", "A","G",
            "C", "A", "G", "A", "C", "G", "T", "A", "A", "G", "T", "C", "A")
 viterbi(hmm,myseq)
 
 ############################################
 
 secventa_exon=c("A","C","G","C","C","T","T","A","A",
                 "A","T","A", "C","C","A","A","A","T","A","G",
                 "A","T","G","G","C","T","T","C","T","A","C")
 secventa_intron=c("A","T","G","G","C","A","T","C","C","T","T",
                   "A","C","A","A","T","G","G","A","G","T","A","T",
                   "C","T","A","G")
 
 
 proportii_exoni=c()
 for (i in 1:4)
   proportii_exoni[i]=sum(secventa_exon==observations[i])/length(secventa_exon)
 proportii_exoni
 
 proportii_introni=c()
 for (i in 1:4)
   proportii_introni[i]=sum(secventa_intron==observations[i])/length(secventa_intron)
 proportii_introni
 
 
 #########################
 # exemplul 4 Monte Carlo - simulate an exponential distribution
simExp=function(n,lambda)
   {u<-runif(n)
   x<-(-1/lambda)*log(1-u)}
par(mfrow=c(1,1))
plot(simExp(1000,2))

x1<-simExp(1000,2)
x2<-rexp(1000,2)
par(mfrow=c(1,2))
hist(x1,30,main="Inverse CDF")
hist(x2,30,main="Direct Simulation")



########################
#exemplul 5 Gibbs pentru bivariata
#Gibbs sampler to simulate a biviariate normal distribution by iteratively sampling from these conditionals.
gibbsBVN=function(x,y, n, rho){
    #create a matrix to store values
    m<-matrix(ncol=2,nrow=n)
    #store initial values in matrix
    m[1,]<-c(x,y)
      #sampling iteration loop
      for (i in 2:n){
      #rnorm takes sd not variance
      #update x conditional on y
      x<-rnorm(1,rho*y,sqrt(1-rho^2))
      #update y conditional on x from above
      y<-rnorm(1,rho*x,sqrt(1-rho^2))
      #store values in matrix
      m[i,]<-c(x,y)}
    m}

par(mfrow=c(2,2))
startX0Y0<-gibbsBVN(0,0,200,0); plot(gibbsBVN(0,0,200,0))
startX5Y5<-gibbsBVN(5,5,200,0); plot(gibbsBVN(5,5,200,0))
startXn5Y5<-gibbsBVN(-5,5,200,0); plot(gibbsBVN(-5,5,200,0))
startXn5Yn5<-gibbsBVN(-5,-5,200,0); plot(gibbsBVN(-5,-5,200,0))


corr0<-gibbsBVN(0,0,1000,0)
corr3<-gibbsBVN(0,0,1000,0.3)
corr5<-gibbsBVN(0,0,1000,0.5)
corr98<-gibbsBVN(0,0,1000,0.98)
par(mfrow=c(2,2))
plot(corr0[,1],corr0[,2],main="XYCorr=0")
plot(corr3[,1],corr3[,2],main="XYCorr=0.3")
plot(corr5[,1],corr5[,2],main="XYCorr=0.5")
plot(corr98[,1],corr98[,2],main="XYCorr=0.98")

par(mfrow=c(2,2))
hist(corr3[,1],nclass=20,main="X Marg, Corr=0.3")
hist(corr3[,2],nclass=20,main="Y Marg, Corr=0.3")
hist(corr98[,1],nclass=20,main="X Marg, Corr=0.98")
hist(corr98[,1],nclass=20,main="Y Marg, Corr=0.98")


##################################
#exemplul 6 Gibbs MCMC grupe de sange
library(BRugs)
setwd("C:/Users/Floarea/Desktop/UVT/MCB/proiect MCB")
modelCheck("bt.txt")
modelData("btData.txt")
modelCompile(numChains=3)
modelInits("btInits.txt")
modelGenInits()


library(boot); library(coda); library(R2WinBUGS)
data=list(n=c(750,250,75,925),total=2000)
inits=list(a=c(1,1),b=c(1,1),o=c(1,1))
bugs(data,inits,model.file="bt.bugs",
     parameters=c("a","b","o"),
     n.chains=2, n.iter=1000,
     bugs.directory="C:/Users/Floarea/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14/")

samplesSet(c("p","q","r"))
modelUpdate(1000)



#####################
# alte exemple

# metrop algorithm
nchain<-1000
y<-3
n<-5
theta<-vector(length=nchain)
theta[1]<-0.04

for(i in 2:nchain){
  thetastar<-runif(1)#new value from the proposal distribution
  new_theta=thetastar^y*(1-thetastar)^(n-y)
  old_theta=(theta[i-1]^y*(1-theta[i-1])^(n-y))
  r<-new_theta/old_theta
  u<-runif(1) #draw a unif val
  if(u<r){
    theta[i]<-thetastar}
  else theta[i]<-theta[i-1]
  }

par(mfrow=c(1,1))
hist(theta[101:1000],nclass=25,prob=T)
xx <- (1:100)/100
lines(xx,dbeta(xx,4,3))

par(mfrow=c(1,3))
plot(1:100,theta[1:100],main="First 100 Runs")
lines(1:100,theta[1:100])
plot(1:1000,theta[1:1000],main="All Runs")
lines(1:1000,theta[1:1000])
plot(901:1000,theta[901:1000],main="Last 100 Runs")
lines(901:1000,theta[901:1000])


# Function to draw histograms of sample results
draw.hist <- function(data, text="") {
  tmp <- hist(data, breaks=0:(max(data)+1), xaxt="n", right=FALSE, freq=TRUE, main=text)
  axis(1, at=tmp$mids, labels=0:max(data))
}
