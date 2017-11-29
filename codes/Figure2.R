
########################
#EBF1 family
library(Logolas)
#read all the motif data from http://compbio.mit.edu/encode-motifs/motifs.txt
library(atSNP)
pfm=LoadMotifLibrary('http://compbio.mit.edu/encode-motifs/motifs.txt',tag = ">",
                     transpose = F, field = 1, sep = c("\t", " ", ">"), skipcols = 1, skiprows = 1,
                     pseudocount = 0)
color_profile = list("type" = "per_row",
                     "col" = RColorBrewer::brewer.pal(4,name ="Spectral"))



EBF1_disc1=t(pfm[[52]]);rownames(EBF1_disc1)=c('A','C','G','T');colnames(EBF1_disc1)=1:ncol(EBF1_disc1)
######################
#logo plot of EBF1
logomaker(EBF1_disc1,color_profile = color_profile,frame_width = 1,pop_name = 'EBF1_disc1')
#neg logo plot
#The log type negative logo plot does not look good


######################
#neg logo plot for the other EBF1 family
#the logoheight is ratio
index=grep('EBF1',names(pfm))


library(grid)
get_viewport_logo(3,2)

seekViewport(paste0("plotlogo", 1))
nlogomaker(EBF1_disc1,logoheight = 'log',color_profile = color_profile,
           frame_width = 1,pop_name = 'EBF1 (disc1) TF', newpage = FALSE,
           control = list(gap_ylab=3.5, epsilon = 0.01, round_off = 0),
           xaxis = FALSE, col_line_split = "white")

#disc2
disc2=t(pfm[[249]]);rownames(disc2)=c('A','C','G','T');colnames(disc2)=1:ncol(disc2)
seekViewport(paste0("plotlogo", 2))
nlogomaker(disc2,logoheight = 'log',color_profile = color_profile,
           frame_width = 1,pop_name = 'EBF1 (disc2) TF', newpage = FALSE,
           control = list(gap_ylab=3.5, epsilon = 0.01, round_off = 0),
           xaxis = FALSE, col_line_split = "white")
#known1
known1=t(pfm[[464]]);rownames(known1)=c('A','C','G','T');colnames(known1)=1:ncol(known1)
seekViewport(paste0("plotlogo", 3))
nlogomaker(known1,logoheight = 'log',color_profile = color_profile,
           frame_width = 1,pop_name = 'EBF1 (known1) TF', newpage = FALSE,
           control = list(gap_ylab=3.5, epsilon = 0.01, round_off = 0),
           xaxis = FALSE, col_line_split = "white")
#known2
known2=t(pfm[[747]]);rownames(known2)=c('A','C','G','T');colnames(known2)=1:ncol(known2)
seekViewport(paste0("plotlogo", 4))
nlogomaker(known2,logoheight = 'log',color_profile = color_profile,
           frame_width = 1,pop_name = 'EBF1 (known2) TF', newpage = FALSE,
           control = list(gap_ylab=3.5, epsilon = 0.01, round_off = 0),
           xaxis = FALSE, col_line_split = "white")
#knwon3
knwon3=t(pfm[[944]]);rownames(knwon3)=c('A','C','G','T');colnames(knwon3)=1:ncol(knwon3)
seekViewport(paste0("plotlogo", 5))
nlogomaker(knwon3,logoheight = 'log',color_profile = color_profile,
           frame_width = 1,pop_name = 'EBF1 (known3) TF', newpage = FALSE,
           control = list(gap_ylab=3.5, epsilon = 0.01, round_off = 0),
           xaxis = FALSE, col_line_split = "white")
#known4
knwon4=t(pfm[[1325]]);rownames(knwon4)=c('A','C','G','T');colnames(knwon4)=1:ncol(knwon4)
seekViewport(paste0("plotlogo", 6))
nlogomaker(knwon4,logoheight = 'log',color_profile = color_profile,
           frame_width = 1,pop_name = 'EBF1 (known4) TF', newpage = FALSE,
           control = list(gap_ylab=3.5, epsilon = 0.01, round_off = 0),
           xaxis = FALSE, col_line_split = "white")

##############  save this image as 15 by 14 in pdf format in R   ########################



