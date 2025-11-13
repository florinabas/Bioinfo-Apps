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



