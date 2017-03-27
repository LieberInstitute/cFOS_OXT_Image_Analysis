##########################
# analyze cFos immuno data
library(ggplot2)
library(reshape2)
library(mclust)
library(beeswarm)
library(jaffelab)
library(lme4)
library(lmerTest)

getValues = function(N) c(mean = mean(N), SEM = sd(N)/sqrt(length(N)), N = length(N))

load('rdas/cFos_pheno.rda')
rownames(pd) = pd$FileID
pd$Animal = factor(as.character(pd$Animal))

datList = lapply(pd$datRda,function(x) {
  load(x);
  list(nuc = nuc,c1 = c1,vol = vol)
  })
names(datList) = pd$FileID

####################
# nuclei intensities & cfos
nucList = lapply(datList,'[[','nuc')
oxtList = lapply(datList,'[[','c1')

load(pd$datRda[1], envir =dat <- new.env())
channels = c('Oxytocin','Cfos') # in order "Cy5" "EGFP"
pd_vars = c('Genotype','Animal','FileID')

###############################
# give nuclei expression labels
nuc = cbind(img = rep(seq(length(nucList)),sapply(nucList,nrow)),
            do.call('rbind',nucList))
nuc = cbind(nuc,apply(nuc[,c(8:9)],2,function(x){
  mod = Mclust(sqrt(x),G=2)
  labs = unique(mod$classification)
  ind = which.max(tapply(x, mod$classification,mean))
  factor(mod$classification,labels = c('Lo','Hi'),
         levels = c(labs[labs!=ind],labs[labs==ind]))
}),pd[ss(rownames(nuc),'\\.'),pd_vars])
names(nuc) = c(names(nuc)[1:7],channels,paste0(channels,'_lab'),pd_vars)

####################
# channel intensities
nucLong = melt(nuc,id.var = c('img','lab','Oxytocin_lab','Cfos_lab',pd_vars))
ggplot(data = subset(nucLong, variable %in% channels),aes(x = value)) + 
  geom_histogram(bins = 100) + facet_wrap(~variable) + xlab('Log10 Intensity')+
  ylab('Count') + scale_x_log10()

##############################
# test nuclei size difference
boxplot(rad~Genotype, main = 'Nuclei Radius',ylab = 'Radius (um)',data = nuc)
summary(lmer(rad~Genotype + (1|Animal) + (1|FileID), data = nuc))

####################
# test number nuclei
datCount = cbind(pd[names(nucList),pd_vars],
                 vol = sapply(datList,'[[','vol'),
                 nNuc = sapply(nucList,nrow))
datCount$nOxt = sapply(datCount$FileID,function(n) with(nuc,sum(FileID==n & Oxytocin_lab =='Hi')))
datCount$nCfos = sapply(datCount$FileID,function(n) with(nuc,sum(FileID==n & Cfos_lab =='Hi')))
datCount$pvnSize = (datCount$vol)^(1/3)
datCount$pvnArea =  (datCount$vol)/h[3]
datCount$pOxt = datCount$nOxt/datCount$nNuc
datCount$pCfos = datCount$nCfos/datCount$nNuc
datCount$nOxtpVol = datCount$nOxt/datCount$pvnSize
datCount$nCfospVol = datCount$nCfos/datCount$pvnSize
datCount$nOxt2 = sapply(datCount$FileID,function(n) nrow(oxtList[[n]]))
datCount$iOxt = sapply(datCount$FileID,function(n) with(oxtList[[n]],
                         mean(intensity*vol)/sum(vol))) #weighted average
datCount$iOxt2 = sapply(datCount$FileID,function(n) with(oxtList[[n]],
                           mean(intensity*vol)))/1000000

# number nuclei
boxplot(nNuc~Genotype, main = '# Nuclei',ylab = 'Count',data = datCount)
beeswarm(nNuc~Genotype,data = datCount,pch = 21, pwbg = datCount$Animal,add = T)
summary(glmer(nNuc~Genotype + pvnSize+ (1|Animal) , data = datCount,family = poisson()))
tapply(datCount$nNuc,datCount$Genotype,getValues)


# pvn size
boxplot(pvnArea~Genotype, main = 'PVN size',ylab = 'Length (um)',data = datCount)
beeswarm(pvnArea~Genotype,data = datCount,pch = 21, pwbg = datCount$Animal,add = T)
summary(lmer(pvnArea/1e6~Genotype + (1|Animal) , data = datCount))
tapply(datCount$pvnArea/1e6,datCount$Genotype,getValues)

########################
# number Oxytocin nuclei
boxplot(nOxt~Genotype, main = '# Oxt(+) Nuclei',ylab = 'Count',data = datCount)
beeswarm(nOxt~Genotype,data = datCount,pch = 21, pwbg = datCount$Animal,add = T)
summary(glmer(nOxt~Genotype+(1|Animal) , data = datCount,family = poisson()))

# number cFos nuclei
boxplot(nCfos~Genotype, main = '# cFos(+) Nuclei',ylab = 'Count',data = datCount)
beeswarm(nCfos~Genotype,data = datCount,pch = 21, pwbg = datCount$Animal,add = T)
summary(glmer(nCfos~Genotype+(1|Animal) , data = datCount,family = poisson()))
tapply(datCount$nCfos/datCount$nNuc,datCount$Genotype,getValues)

#####################################
# % of nuclei Oxt or cFos (+)
boxplot(pOxt~Genotype, main = '% Oxt(+) nuclei',ylab = 'Fraction',data = datCount)
beeswarm(pOxt~Genotype,data = datCount,pch = 21, pwbg = datCount$Animal,add = T)
summary(lmer(pOxt~Genotype+ (1|Animal) , data = datCount))
tapply(datCount$pOxt,datCount$Genotype,getValues)

# number cFos nuclei per PVN size
boxplot(pCfos~Genotype, main = '% cFos(+) nuclei',ylab = 'Fraction',data = datCount)
beeswarm(pCfos~Genotype,data = datCount,pch = 21, pwbg = datCount$Animal,add = T)
summary(lmer(pCfos~Genotype+ (1|Animal) , data = datCount))
tapply(datCount$pCfos,datCount$Genotype,getValues)

#####################################
# number Oxytocin nuclei per PVN size
boxplot(nOxtpVol~Genotype, main = '# Oxt(+) per PVN size',ylab = 'Density (#/um)',data = datCount)
summary(lmer(nOxtpVol~Genotype+ (1|Animal) , data = datCount))

# number cFos nuclei per PVN size
boxplot(nCfospVol~Genotype, main = '# cFos(+) per PVN size',ylab = 'Density (#/um)',data = datCount)
summary(lmer(nCfospVol~Genotype+ (1|Animal) , data = datCount))

##########################
# oxt sizes, oxytocin cells
oxt = cbind(img = rep(seq(length(oxtList)),sapply(oxtList,nrow)),
            do.call('rbind',oxtList))
oxt = cbind(oxt,pd[ss(rownames(oxt),'\\.'),pd_vars])

########################
# number Oxytocin cells
boxplot(nOxt2~Genotype,main='# Oxt cells',ylab='Count',data = datCount)
beeswarm(nOxt2~Genotype,add =T,data = datCount,pch =21,pwbg = datCount$Animal)
summary(glmer(nOxt2~Genotype + (1|Animal),data = datCount,family = poisson))

########################
# mean intensity of Oxt
boxplot(iOxt~Genotype,main='Mean Oxt intensity',ylab='Weighted Avg. Intensity',data = datCount)
beeswarm(iOxt~Genotype,add =T,data = datCount,pch =21,pwbg = datCount$Animal)
summary(lmer(iOxt~Genotype + (1|Animal),data = datCount))

boxplot(iOxt2~Genotype,main='Sum Oxt intensity',ylab='Total Intensity (1e6)',data = datCount)
beeswarm(iOxt2~Genotype,add =T,data = datCount,pch =21,pwbg = datCount$Animal)
summary(lmer(iOxt2~Genotype+ (1|Animal),data = datCount))

########################
# Oxytocin cell diameter
boxplot(diam~Genotype,main='Oxytocin Cell Diameter',ylab='Diameter (um)',data = oxt)
summary(lmer(diam~Genotype +(1|FileID)+ (1|Animal),data = oxt))

boxplot(vol~Genotype,main='Oxytocin Cell Volume',ylab='Volume (um^3)',data = oxt,log ='y')
summary(lmer(vol~Genotype +(1|FileID)+ (1|Animal),data = oxt))

boxplot(log(intensity+1)~Genotype,main='Oxytocin Intensity',
        ylab='Volume (um^3)',data = oxt)
summary(lmer(log(intensity+1)~Genotype +(1|FileID)+ (1|Animal),data = oxt))

ggplot(data = nuc,aes(x = Oxytocin,fill = Genotype))+
  geom_density(alpha = .5)+scale_x_log10()+xlab('Oxytocin Intensity')
summary(lmer(log(Oxytocin+1)~Genotype +(1|FileID) + (1|Animal),data = nuc))

ggplot(data = nuc,aes(x = Cfos,fill = Genotype))+
  geom_density(alpha = .5)+scale_x_log10()+xlab('cFOS Intensity')
summary(lmer(log(Cfos+1)~Genotype +(1|FileID) + (1|Animal),data = nuc))

##########################
# percent of cellA+ that are cellB+
datCount$nCfosPnOxt = sapply(datCount$FileID,function(n) 
  with(nuc,sum(FileID==n & Cfos_lab =='Hi' & Oxytocin_lab =='Hi')))/datCount$nOxt
datCount$nOxtPncFos = sapply(datCount$FileID,function(n) 
  with(nuc,sum(FileID==n & Cfos_lab =='Hi' & Oxytocin_lab =='Hi')))/datCount$nCfos

################################
# % of cFOS+ cells that are OXT+ 
boxplot(nCfosPnOxt~Genotype, main = '% Oxt(+) that are cFos(+)',ylab = 'Fraction',data = datCount)
beeswarm(nCfosPnOxt~Genotype,data = datCount,pch = 21, 
         pwbg = datCount$Animal,add = T,corral = 'wrap')
summary(glmer(nCfosPnOxt~Genotype+(1|Animal) , data = datCount,family = binomial()))
tapply(datCount$nCfosPnOxt,datCount$Genotype,getValues)

# % of Oxt+ cells that are cFOS+ 
boxplot(nOxtPncFos~Genotype, main = '% cFos(+) that are Oxt(+)',ylab = 'Fraction',data = datCount)
beeswarm(nOxtPncFos~Genotype,data = datCount,pch = 21, pwbg = datCount$Animal,add = T)
summary(glmer(nOxtPncFos~Genotype+ (1|Animal) , data = datCount,family = binomial()))
tapply(datCount$nOxtPncFos,datCount$Genotype,getValues)

summary((lm2 = lmer(Cfos~Oxytocin+(1|Animal) +(1|FileID), data = nuc)))
plot(Cfos~Oxytocin, data = nuc, col = c('#00000020','#FF000020')[as.numeric(Genotype)],
     pch=20, log = 'xy',xlab ='Oxytocin Intensity',ylab = 'cFos Intensity')




