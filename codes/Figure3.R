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

library(grid)
grid.newpage()
layout.rows <- 1
layout.cols <- 2
top.vp <- viewport(layout=grid.layout(layout.rows, layout.cols,
                                      widths=unit(rep(6,layout.cols), rep("null", 2)),
                                      heights=unit(c(20, 20), rep("lines", 2))))

plot_reg <- vpList()
l <- 1
for(i in 1:layout.rows){
  for(j in 1:layout.cols){
    plot_reg[[l]] <- viewport(layout.pos.col = j, layout.pos.row = i, name = paste0("plotlogo", l))
    l <- l+1
  }
}


plot_tree <- vpTree(top.vp, plot_reg)

pushViewport(plot_tree)

seekViewport(paste0("plotlogo", 1))
#logo plot of EBF1
logomaker(EBF1_disc1,color_profile = color_profile,
          frame_width = 1,pop_name = 'EBF1_disc1', newpage = FALSE)
#neg logo plot
#The log type negative logo plot does not look good

seekViewport(paste0("plotlogo", 2))
nlogomaker(EBF1_disc1,logoheight = 'log',color_profile = color_profile,
           frame_width = 1,pop_name = 'EBF1 (disc1) TF', newpage = FALSE,
           control = list(gap_ylab=3.5, epsilon = 0.01))
