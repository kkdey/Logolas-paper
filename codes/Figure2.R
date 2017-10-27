
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
Logolas::nlogomaker(mat,
                    logoheight = 'log',
                    bg = bgmat,
                    color_profile = color_profile,
                    frame_width = 1, newpage = FALSE,
                    xlab = "Position",
                    pop_name = 'Histone marks (EDLogo : log)',
                    control = list(epsilon=0.1,
                                   gap_ylab=3.5))

seekViewport(paste0("plotlogo", 2))
Logolas::nlogomaker(mat,
                    logoheight = 'unscaled_log',
                    bg = bgmat,
                    color_profile = color_profile,
                    frame_width = 1, newpage = FALSE,
                    xlab = "Position",
                    pop_name = 'Histone marks (EDLogo :  unscaled log)',
                    control = list(epsilon=0.1,
                                   gap_ylab=3.5))

