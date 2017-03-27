##########################################
# script to read and extract RNAscope data
library(matlabr)
library(jaffelab)
library(EBImage)
library(R.matlab)

options(matlab.path='~/bin/MATLAB')

disp = function(x,all=T,...) EBImage::display(x,method = 'raster',all=all,...)
#########################
# load the phenotype data
load('rdas/cFos_pheno.rda')

###################################################
# extract RNAscope data with custom matlab function
for (i in seq(nrow(pd))){
  cat(paste0('Preprocessing ',i,' out of ',nrow(pd),' images.\n'))
  code = c(paste0("filename='",pd$CZI[i],"';"),
           '[tvmask, hypomask, nmask, img, dat] = cellSegCFos(filename);',
           paste0("saveImg = '",pd$imgMat[i],"';"),
           paste0("saveDat = '",pd$datMat[i],"';"),
           "save(saveImg,'img','tvmask','hypomask','nmask');",
           "save(saveDat,'dat');")
  res = run_matlab_code(code)
  
}

###############################
# reformat RNAscope data into R
for (i in seq(nrow(pd))){
  ######################
  # read in matlab file
  dat = readMat(pd$datMat[i])$dat
  
  ################
  # extract tables
  cName = unlist(dat[1,,]);
  h = unlist(dat[2,,]);
  nuc = data.frame(dat[3,,])
  names(nuc) = c('x','y','z','lab','vol','rad','c1','c2','c3')
  nuc = nuc[nuc$vol > 100 & nuc$vol <1000,]
   c1 = data.frame(dat[4,,])
  names(c1) = c('x','y','z','lab','vol','rad')
  
  c2 = data.frame(dat[5,,])
  names(c2) = c('x','y','z','lab','vol','rad')
  
  save(cName,h,nuc,c1,c2,file= pd$datRda[i])
}


table(file.exists(pd$datRda))