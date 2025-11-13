library(readxl) #sessionInfo()

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
setwd("C:/Users/Floarea/Desktop/UVT/sem 1/biostat si R/proiect")
data <- read_excel("date_migrena.xlsx")
str(data) #informatii despre continut
names(data) #numele variabilelor
dim(data)

# pregatirea datelor
str(data)
#pre-freq - chr desi ar trebui numar

# redenumire coloane cu spatii
colnames(data)[2:8] = c("age","gender", "aura","pre_freq", "pre_intens", "post_freq", "post_intens")

# modificare codificare gen din 1/2 in F/M si aura din 1/2/3 in A/N/AN
data$gender[data$gender==1]="F";data$gender[data$gender==2]="M"
data$aura[data$aura==1]="A"; data$aura[data$aura==2]="N";data$aura[data$aura==3]="AN"

# eliminare NAN (pre-tx freq contine o variabila .)
control_patient=data[which(data$pre_freq=="."),]
data=data[-which(data$pre_freq==".") ,]

# transformare formate
class(data$pre_freq) #character
data$pre_freq=as.numeric(data$pre_freq) #class(data$pre_freq) #numeric
data$gender=as.factor(data$gender)
data$aura=as.factor(data$aura)

str(data)
attach(data)

# grafic varsta in functie de gen
install.packages("ggplot2")
library(ggplot2)
ggplot(aes(x=gender, y=age), data=data) + 
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=4)

table(factor(gender))

summary(subset(age,gender=="F"))
summary(subset(age,gender=="M"))

# grafic varsta in functie de aura
ggplot(aes(x=aura, y=age), data=data) + 
  geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=4)
table(factor(aura))



# analiza grafica eficienta tratament - histograme
par(mfrow=c(2, 2)); par(mar=c(5, 2, 2, 2))
hist(pre_freq);hist(post_freq);hist(pre_intens);hist(post_intens);

# intensitate pre vs post tratament
length(post_intens[post_intens==0])
pre_intens[post_intens==0]

#comparatii in functie de gen
data_f=subset(data,gender=="F") #nrow(data_f)#23
data_m=subset(data,gender=="M") #nrow(data_m)#5
c(mean(data_f$pre_freq), mean(data_f$pre_intens)) #5.74, 8.87
c(mean(data_m$pre_freq),mean(data_m$pre_intens)) #6.4, 8.6
rm(data_f,data_m)

boxplot(pre_freq~gender, main="pre_freq") #posibil outlier (19)
boxplot(post_freq~gender, main="post_freq")

boxplot(pre_intens~gender, main="pre_intens")
boxplot(post_intens~gender, main="post_intens")

# scapam de outlier
data_fo=data[-which.max(data$pre_freq),]
cor(data_fo$age,data_fo$pre_freq) #0.21

#vedem daca ramane aceeasi distributie si post tratament
qqplot(pre_freq,post_freq)
qqplot(data_fo$pre_freq,data_fo$post_freq)


##############################
#comparare medii pre post tratament
c(mean(pre_freq), mean(post_freq)) 
c(mean(pre_intens),mean(post_intens)) 

# verificam distributia frecventelor
par(mfrow=c(2, 2)); par(mar=c(4, 2, 3, 2))
plot(density(pre_freq))
plot(density(post_freq))
plot(density(pre_intens))
plot(density(post_intens))
# asimetrie f pronuntata - n ar trebui sa merg pe parametric

# verificam normaltiatea
# H0 - esantionul vine dintr-o distributie normala

# Testul Shapiro-Wilk (volum esantion <5000)
shapiro.test(pre_freq);shapiro.test(post_freq);shapiro.test(pre_intens);shapiro.test(post_intens);

#Testul Anderson-Darling (volum esantion>7)
library(nortest)
ad.test(pre_freq);ad.test(post_freq);ad.test(pre_intens);ad.test(post_intens)
#pval<0.05 ->resping H0 -> nu sunt din distributie normala


# distributie non-normala si esantioane mici - teste neparametrice 
# esantioane dependente (perechi) - testul Wilcoxon 
# H0: miu1=miu2, H1: miu1>miu2
wilcox.test(pre_freq,post_freq,alternative="greater",paired=TRUE)
# p val < 0.05 ->  resping ipoteza nula -> frecventa este mai ridicata inainte de tratament
wilcox.test(pre_intens,post_intens,alternative="greater",paired=TRUE)
# p val < 0.05 ->  resping ipoteza nula -> intensitatea este mai ridicata inainte de tratament

#anova vs kruskal
par(mfrow=c(1, 2)); 
plot(pre_freq,pre_intens)
plot(data_fo$pre_freq,data_fo$pre_intens)

#anaizam variantele celor 4 grupuri - cu outlier
x=numeric(4)
for (i in 7:10){
  x[i-6]=sd(pre_freq[pre_intens==i])}

# fara outlier
x=numeric(4)
for (i in 7:10){
  x[i-6]=sd(data_fo$pre_freq[data_fo$pre_intens==i])}

rm(i,x)
#anova - comp medii >2 populatii
#F distribution of the ratio of variances -can be used to analyze signal to noise ratios
#signal- between variation (changes in the mean response due to treatments)
# noise - within variability (error)
#for a test statistic the ratio of the MSB/MSE - signal to noise ratio


#conditii - obs indep, variante aprox egale, distr reziduuri aprox normala
# obs sunt indep - provind de la pacienti diferiti
#compar dispersiile
library(car)
# testul Levene - populatii non-normale si mai mult de doua grupuri
leveneTest(pre_freq~as.factor(pre_intens))
#pval>0.05 -> nu resping H0 conform careia dispersiile sunt egale
#pot aplica ANOVA
leveneTest(data_fo$pre_freq~as.factor(data_fo$pre_intens))


summary(aov(pre_freq~as.factor(pre_intens)))
#pval<0.05 -> resping H0 conform mediile sunt egale
pairwise.t.test(pre_freq,as.factor(pre_intens))
# grupa de intensitate 7 este diferita de restul

summary(aov(data_fo$pre_freq~as.factor(data_fo$pre_intens)))
#pval>0.05 - nu pot respinge H0 - mediile sunt egale
pairwise.t.test(data_fo$pre_freq,as.factor(data_fo$pre_intens))
pairwise.t.test(data_fo$pre_freq,as.factor(data_fo$pre_intens),p.adj="fdr")
#False Discovery Rate (FDR)
#a less restricting adjustment of the a-level that can be applied to a larger number of simultaneous tests
#plot the spread of data values for each group
par(mfrow=c(1, 2));
stripchart(pre_freq~as.factor(pre_intens))
stripchart(data_fo$pre_freq~as.factor(data_fo$pre_intens))

par(mfrow=c(2,2))
plot(aov(data_fo$pre_freq~as.factor(data_fo$pre_intens)))
plot(aov(pre_freq~as.factor(pre_intens)))


#regresie
par(mfrow=c(1, 1));
plot(post_freq,post_intens)
cor(post_freq,post_intens)

#m1=lm(post_intens~post_freq+gender+age+aura);summary(m1)
m1=lm(post_intens~post_freq+age);summary(m1)
m1=lm(post_intens~post_freq+gender+age);summary(m1)
m1=lm(post_freq~age);summary(m1)

m1=lm(post_intens~post_freq+gender);summary(m1)
#m1=lm(post_freq~gender);summary(m1) # test - ok

par(mfrow=c(1, 1));
boxplot(post_intens~gender, main="post_intens")
plot(lm(post_intens~post_freq+gender))

res=residuals(m1)
plot(fitted(m1), res)
abline(0,0)


qqnorm(res)
qqline(res) 
plot(density(res))

#testul de normalitate pentru reziduuri
plot(m1$residuals)
hist(m1$residuals)
plot(m1$fitted)
shapiro.test(m1$residuals)

#predict
datenoi=data.frame(post_freq=control_patient$post_freq,gender=control_patient$gender)
predicted=predict(m1, newdata = datenoi, interval="predict",level=0.95)
actual=control_patient$post_intens
c(predicted[1],actual)
predicted
actual>=predicted[1]&actual<=predicted[3]




#logit - prezenta aura
# date_a -> fara a treia stare AN
# comparatii in functie de aura
data_a=data[-which(data$aura=="AN") ,]
class(data_a$aura)
data_a$aura=as.character(data_a$aura)
data_a$aura[as.character(data_a$aura)=="A"]=1
data_a$aura[as.character(data_a$aura)=="N"]=0
data_a$aura=as.numeric(data_a$aura)

par(mfrow=c(2, 2)); par(mar=c(5, 2, 2, 2))
boxplot(data_a$pre_freq~data_a$aura, main="pre_freq")
boxplot(data_a$post_freq~data_a$aura, main="post_freq")
boxplot(data_a$pre_intens~data_a$aura, main="pre_intens")
boxplot(data_a$post_intens~data_a$aura, main="post_intens")

# distributie non-normala si esantioane mici - teste neparametrice 
# esantioane independente 
# H0: miu1=miu2, H1: miu1>miu2
pre_freq_a=data_a$pre_freq[data_a$aura==1] #14
pre_freq_na=data_a$pre_freq[data_a$aura==0] #13

par(mfrow=c(1, 1));par(mar=c(4, 4, 4, 4))
plot(density(pre_freq_a));lines(density(pre_freq_na),col="red")
shapiro.test(pre_freq_a); shapiro.test(pre_freq_na)
#P val<0.05 - resping H0 - nu provin din distributie normala

wilcox.test(pre_freq_a,pre_freq_na,alternative="less",paired=FALSE)
#pval<0.05 -> resping H0 -> media pre freq cu aura e mai mica 


# https://stackoverflow.com/questions/27464893/getting-warning-newdata-had-1-row-but-variables-found-have-32-rows-on-pred
#data_a$aura=as.factor(data_a$aura)
#plot(data_a$pre_freq, data_a$aura)
m2=glm(data_a$aura~data_a$pre_freq +data_a$age+data_a$gender,family=binomial(logit),data=data_a); summary(m2)
b=data_a$pre_freq #b=data_a$post_intens
m2=glm(data_a$aura~b,family=binomial(logit),data=data_a); summary(m2)
# odds-ul pentru prezenta aurei scade cu 0.5 unitati cu fiecare unitate in plus a frecventei pre-tratament

#  P = e(beta0 + beta1 X+ epsilor i)/ (e(beta0 + beta1 X+ epsilor i) +1)
odds=exp(cbind(OR = coef(m2), confint(m2)))
#p/1-p=0.6 ->p=0.6-0.6*p ->
p=0.6/1.6 
# probabilitatea de a avea aura scade cu 37,5% cu cresterea unei unitati de frecventa

for (i in 1:10)
{predicted[i]=predict(m2,data.frame(b=i),type="response")}
predicted
plot(predicted)

#chi2
m2.nul=glm(data_a$aura~1,family=binomial(logit),data=data_a); summary(m2)
#install.packages("lmtest")
library(lmtest)
lrtest(m2,m2.nul)
#pval<0.05 -> resping H0 - modelul e ok
hoslem.test(as.integer(data_a$aura[1:10]),fitted(m2)[1:10])


library(boot)
d <- data.frame(x=data_a$pre_freq, y=data_a$aura)
m <- glm(y ~ x, data=d)
m.cv <- cv.glm(d, m, K=2) # Sometimes succeeds
acuratete=1-m.cv$delta[2]
acuratete




#####################
#altele

# sort vs order
sort(data$age) #valorile sortate
order(data$age) #pozitiile vectorului aferente sortarii

#install.packages("carData")
# library(carData)
# summary(KosteckiDillon)
# data(KosteckiDillon)


# piechart in functie de varsta
slices <- c(length(gender[gender=="F"])/length(gender), length(gender[gender=="M"])/length(gender))
lbls <- c("F","M")
pie(slices, labels = lbls)

#kruscal
kruskal.test(pre_freq~as.factor(pre_intens))
#p-value = 0.05521>0.05 ->nu putem respinge H0 conform careia mediile de la nivelul celor 4 grupuri de intensitate sunt egale 
kruskal.test(data_fo$pre_freq~as.factor(data_fo$pre_intens))


#grupul cu 7 - prea mic -il putem grupa cu 8
data_k=data_fo
data_k$pre_intens[data_k$pre_intens==7]=8
table(data_k$pre_intens)
x=numeric(3)
for (i in 8:10){
  x[i-7]=sd(data_k$pre_freq[data_k$pre_intens==i])}

pairwise.t.test(data_k$pre_freq,as.factor(data_k$pre_intens))
#pval>0.05 - nu pot respinge H0 - mediile sunt egale

data_a=data[-which(data$aura=="AN") ,]
class(data_a$aura)
data_a$aura=as.character(data_a$aura)

boxplot(data_a$pre_freq~data_a$aura, main="pre_freq")
boxplot(data_a$post_freq~data_a$aura, main="post_freq")

boxplot(data_a$pre_intens~data_a$aura, main="pre_intens")
boxplot(data_a$post_intens~data_a$aura, main="post_intens")

prop.table(table(data_a$pre_intens,data_a$aura),1)
prop.table(table(data_a$post_intens,data_a$aura),1)