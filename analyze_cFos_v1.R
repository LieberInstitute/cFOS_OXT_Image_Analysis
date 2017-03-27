##########################
# analyze cFos immuno data
library(ggplot2)
library(reshape2)
library(mclust)
library(beeswarm)
load('rdas/cFos_pheno.rda')

datList = lapply(pd$datRda,function(x) {
  load(x);
  list(nuc = nuc,c1 = c1,c2 = c2)
  })

load(pd$datRda[1], envir =dat <- new.env())
channels = c('Oxytocin','tDtomato','Cfos') # in order "Cy5" "DsRed"  "EGFP"
####################
# nuclei intensities
nucList = lapply(datList,'[[','nuc')
nuc = cbind(img = rep(seq(length(nucList)),sapply(nucList,nrow)),
            do.call('rbind',nucList))

###############################
# give nuclei expression labels
nuc = cbind(nuc,apply(nuc[,c(8:10)],2,function(x){
  mod = Mclust(log10(x+1),G=2)
  labs = unique(mod$classification)
  ind = which.max(tapply(x, mod$classification,mean))
  factor(mod$classification,labels = c('Lo','Hi'),
         levels = c(labs[labs!=ind],labs[labs==ind]))
  }))

names(nuc) = c(names(nuc)[1:7],channels,paste0(channels,'_lab')) 


#########################
# nuclei bodies
boxplot(nuc$rad, main = 'Nuclei Radius',ylab = 'Radius (um)')

####################
# channel intensities
nucLong = melt(nuc,id.var = c('img','lab'))
ggplot(data = subset(nucLong, variable %in% channels),aes(x = value)) + 
  geom_histogram() + facet_wrap(~variable) + xlab('Log10 Intensity')+
  ylab('Count') + scale_x_log10()

################################
# % cells with cFos and oxytocin
apply(nuc[,paste0(channels,'_lab')],2,function(x) sum(x =='Hi'))/nrow(nuc) # each
sum(apply(nuc[,paste0(channels,'_lab')],1,function(x) all(x =='Hi')))/nrow(nuc) # all
sum(apply(nuc[,c('Oxytocin_lab','Cfos_lab')],1,function(x) all(x =='Hi')))/nrow(nuc)
sum(apply(nuc[,c('Oxytocin_lab','tDtomato_lab')],1,function(x) all(x =='Hi')))/nrow(nuc)
sum(apply(nuc[,c('Cfos_lab','tDtomato_lab')],1,function(x) all(x =='Hi')))/nrow(nuc)

##########################
# scatterplot matrix
nuc2 = nuc[apply(nuc[,paste0(channels,'_lab')],1,function(x) all(x =='Hi')),]
nucLong2 = melt(nuc2,id.var = c('img','lab'))
pairs(nuc2[, channels],log = 'xy')


####################
# c1 sizes, oxytocin cells
c1List = lapply(datList,'[[','c1')
c1 = cbind(img = rep(seq(length(c1List)),sapply(c1List,nrow)),
            do.call('rbind',c1List))
boxplot(c1$rad,main = 'Oxytocin Cell Radius',ylab = 'Radius (um)')
nrow(c1)

####################
# c2 sizes, tDtomato cells
c2List = lapply(datList,'[[','c2')
c2 = cbind(img = rep(seq(length(c2List)),sapply(c2List,nrow)),
           do.call('rbind',c2List))
boxplot(c2$rad,main = 'tDtomato Cell Radius',ylab = 'Radius (um)')
nrow(c2)

