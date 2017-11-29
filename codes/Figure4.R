
#################           Figure 4           ########################


library(Logolas)
color_profile = list("type" = "per_row",
                     "col" = RColorBrewer::brewer.pal(4,name ="Spectral"))
pwm=rbind(0.328,0.332,0.33,0.01)
rownames(pwm)=c('A','C','G','T')
colnames(pwm)=1:ncol(pwm)

library(grid)
grid.newpage()
layout.rows <- 1
layout.cols <- 3
top.vp <- viewport(layout=grid.layout(layout.rows, layout.cols,
                                      widths=unit(rep(8,layout.cols), rep("null", 2)),
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
logomaker(pwm,color_profile=color_profile,frame_width = 1,
          newpage = FALSE,
          col_line_split="white",
          xlab = "",
          pop_name = '',control = list(gap_ylab=3.5, totbins = 3))

seekViewport(paste0("plotlogo", 2))
nlogomaker(pwm,logoheight = 'log', color_profile = color_profile,
           frame_width = 1,pop_name = '', newpage = FALSE,
           col_line_split="white",
           xlab = "",
           control = list(gap_ylab=3.5, epsilon = 0.1, negbins = 3,
                          posbins = 1, round_off = 1))

seekViewport(paste0("plotlogo", 3))
nlogomaker(pwm,logoheight = 'wKL', color_profile = color_profile,
           col_line_split="white",
           xlab = "",
          frame_width = 1,pop_name = '', newpage = FALSE,
                    control = list(quant = 0, gap_ylab=3.5, epsilon = 0.1,
                                   negbins = 3, posbins = 2, round_off = 1))

