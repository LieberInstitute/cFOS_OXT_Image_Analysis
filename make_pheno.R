# make cFos phenotype
library(jaffelab)
options(stringsAsFactors = F)

#########################
# find all the files
dir = paste0('/dcl01/lieber/ajaffe/Keri/Imaging/cFos/Stitched_for_Badoi')
fns = ss(list.files(dir,pattern = '.czi',recursive = T),'/',2)

##########################
# make the phenotype table
pd = data.frame(FileID = gsub('e1 ','',ss(fns,'.czi')))
pd$CZI = list.files(dir,pattern = '.czi',full.names = T, recursive = T)
pd$Genotype = factor(ss(pd$FileID,'_'),levels = c('WT','KI'))
pd$Animal = ss(pd$FileID,'_',2)
pd$imgMat = paste0('/dcl01/lieber/ajaffe/Keri/Imaging/cFos/rdas/',pd$FileID,'_img.mat')
pd$datMat = paste0('/dcl01/lieber/ajaffe/Keri/Imaging/cFos/rdas/',pd$FileID,'_dat.mat')
pd$datRda = paste0('/dcl01/lieber/ajaffe/Keri/Imaging/cFos/rdas/',pd$FileID,'_dat.rda')

############################
# save the data
save(pd,file = 'rdas/cFos_pheno.rda')