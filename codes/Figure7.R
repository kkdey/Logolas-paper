library(pmsignature)
library(Logolas)
cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category ==
                                       'qual',]
col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))
total_chars = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O",
                "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "zero", "one", "two",
                "three", "four", "five", "six", "seven", "eight", "nine", "dot", "comma",
                "dash", "colon", "semicolon", "leftarrow", "rightarrow")

set.seed(20)
color_profile <- list("type" = "per_symbol",
                      "col" = sample(col_vector, length(total_chars), replace=FALSE))
load('../data/Lymphoma-B-cell.4.Rdata')

##############################
visPMSignature(resultForSave[[1]],2,isScale=T)
################################
mat=resultForSave[[1]]@signatureFeatureDistribution[2,,]
mat1=cbind(t(mat[2:3,1:4]),rep(NA,4),t(mat[4:5,1:4]))
rownames(mat1)=c('A','C','G','T')
colnames(mat1) = c("-2", "-1", "0", "1", "2")
mat2=cbind(rep(NA,6),rep(NA,6),(mat[1,]),rep(NA,6),rep(NA,6))
colnames(mat2) = c("-2", "-1", "0", "1", "2")
rownames(mat2) = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
table = rbind(mat1, mat2)

#logomaker(table,color_profile = color_profile,frame_width = 1,pop_name = paste('logo plot',name))

nlogomaker(table,
           logoheight = 'log',
           color_profile = color_profile,
           frame_width = 1,
           xlab = "Position",
           pop_name = '',
           #control = list(epsilon=0.25,gap_ylab=3.5)
)
