

##########################    PANEL   (A)    #####################################

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

#neg logo plot
#The log type negative logo plot does not look good

bgmat <- matrix(0.25, dim(EBF1_disc1)[1], dim(EBF1_disc1)[2])
rownames(bgmat) <- rownames(EBF1_disc1)
colnames(bgmat) <- colnames(EBF1_disc1)

library(grid)
grid.newpage()
layout.rows <- 4
layout.cols <- 3
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
Logolas::logomaker(EBF1_disc1,color_profile = color_profile,frame_width = 1,
          pop_name = 'EBF1_disc1 (Logo)', newpage = FALSE)

seekViewport(paste0("plotlogo", 2))
Logolas::nlogomaker(EBF1_disc1,logoheight = 'log',color_profile = color_profile,
           frame_width = 1,pop_name = 'EBF1_disc1 (EDLogo)',
           bg = bgmat, newpage = FALSE,
           control = list(gap_ylab=3.5, epsilon = 0.01))

seekViewport(paste0("plotlogo", 3))
Logolas::nlogomaker(EBF1_disc1,logoheight = 'wKL', color_profile = color_profile,
           frame_width = 1,pop_name = 'EBF1_disc1 (wKL-Logo)',
           bg = bgmat, newpage = FALSE,
           control = list(quant = 0, gap_ylab=3.5, epsilon = 0.1))


######################    PANEL  (B)   ############################################

proteins_pwm <- get(load("../data/protein_pwm.RData"))
m <- proteins_pwm$pwm00389


cols1 <- c(rev(RColorBrewer::brewer.pal(12, "Paired"))[c(3,4,7,8,11,12,5,6,9,10)],
           RColorBrewer::brewer.pal(12, "Set3")[c(1,2,5,8,9)],
           RColorBrewer::brewer.pal(9, "Set1")[c(9,7)],
           RColorBrewer::brewer.pal(8, "Dark2")[c(3,4,8)])

color_profile <- list("type" = "per_row",
                      "col" = cols1)

seekViewport(paste0("plotlogo", 4))
Logolas::logomaker(m,color_profile = color_profile,frame_width = 1,
                   pop_name = 'PF00389 motif 2 (Logo)', newpage = FALSE)

seekViewport(paste0("plotlogo", 5))
Logolas::nlogomaker(m, logoheight = "log",
           color_profile = color_profile,
           frame_width = 1, pop_name = "PF00389 motif 2 (EDLogo)",
           newpage = FALSE,
           control = list(quant = 0.5,
                          depletion_weight = 0, epsilon = 0.1))

seekViewport(paste0("plotlogo", 6))
Logolas::nlogomaker(m, logoheight = "wKL",
                    color_profile = color_profile,
                    frame_width = 1, pop_name = "PF00389 motif 2 (wKL-Logo)",
                    newpage = FALSE,
                    control = list(quant = 0,
                                   depletion_weight = 0, epsilon = 0.1))



######################    PANEL  (C)   ############################################

load("../data/Lymphoma-B-cell.4.RData")
mat=resultForSave[[1]]@signatureFeatureDistribution[2,,]
mat1=cbind(t(mat[2:3,1:4]),rep(NA,4),t(mat[4:5,1:4]))
rownames(mat1)=c('A','C','G','T')
colnames(mat1) = c("-2", "-1", "0", "1", "2")
mat2=cbind(rep(NA,6),rep(NA,6),(mat[1,]),rep(NA,6),rep(NA,6))
colnames(mat2) = c("-2", "-1", "0", "1", "2")
rownames(mat2) = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
table = rbind(mat1, mat2)

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

seekViewport(paste0("plotlogo", 7))
Logolas::logomaker(table,
          color_profile = color_profile,
          frame_width = 1, newpage = FALSE,
          xlab = "Position",
          pop_name='B cell mutations (Logo)')

seekViewport(paste0("plotlogo", 8))
Logolas::nlogomaker(table,
           logoheight = 'log',
           color_profile = color_profile,
           frame_width = 1, newpage = FALSE,
           xlab = "Position",
           pop_name = 'B cell mutations (EDLogo)',
           control = list(epsilon=0.1,
                          gap_ylab=3.5))

seekViewport(paste0("plotlogo", 9))
Logolas::nlogomaker(table,
                    logoheight = 'wKL',
                    color_profile = color_profile,
                    frame_width = 1, newpage = FALSE,
                    xlab = "Position",
                    pop_name = 'B cell mutations (wKL-Logo)',
                    control = list(quant = 0, epsilon=0.1,
                                   gap_ylab=3.5))



########################   PANEL (D)   ###################################


mat <- rbind(c(326, 296, 81, 245, 71),
             c(258, 228, 55, 273, 90),
             c(145, 121, 29, 253, 85),
             c(60, 52, 23, 180, 53),
             c(150, 191, 63, 178, 63))

bgmat <- rbind(c(542, 218, 118, 108, 33),
               c(480, 204, 107, 87, 26),
               c(346, 131, 71, 66, 20),
               c(199, 72, 41, 43, 13),
               c(368, 101, 67, 80, 30))

rownames(mat) <- c("H3K4ME1", "H3K4ME2", "H3K4ME3", "H3AC", "H4AC")
colnames(mat) <- c("Intergenic","Intron","Exon",
                   "Gene start ","Gene end")

set.seed(201)
color_profile <- list("type" = "per_row",
                      "col" = sample(RColorBrewer::brewer.pal(10,name = "Spectral"),
                                     dim(mat)[1]))

seekViewport(paste0("plotlogo", 10))
Logolas::logomaker(mat,
          color_profile = color_profile,
          frame_width = 1,
          ic.scale = TRUE,
          pop_name = "Histone marks (Logo)",
          xlab = "", newpage = FALSE,
          yscale_change = TRUE,
          col_line_split = "black",
          control = list(gap_ylab = 3.5))


seekViewport(paste0("plotlogo", 11))
Logolas::nlogomaker(mat,
                    logoheight = 'unscaled_log',
                    bg = bgmat,
                    color_profile = color_profile,
                    frame_width = 1, newpage = FALSE,
                    xlab = "Position",
                    pop_name = 'Histone marks (EDLogo, unscaled)',
                    control = list(epsilon=0.1,
                                   gap_ylab=3.5))


seekViewport(paste0("plotlogo", 12))
Logolas::nlogomaker(mat,
                    logoheight = 'unscaled_log',
                    bg = bgmat,
                    color_profile = color_profile,
                    frame_width = 1, newpage = FALSE,
                    xlab = "Position",
                    pop_name = 'Histone marks (wKL-Logo, unscaled)',
                    control = list(quant = 0, epsilon=0.1,
                                   gap_ylab=3.5))


##########  save this file in R pdf as 22 by 17   ########################


