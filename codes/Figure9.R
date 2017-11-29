

###############   PSSM EDLogo plot (Supplementary figure 7) ##################


pssm <- get(load("../data/pssm_00389_motif2.RData"))

cols1 <- c(rev(RColorBrewer::brewer.pal(12, "Paired"))[c(3,4,7,8,11,12,5,6,9,10)],
           RColorBrewer::brewer.pal(12, "Set3")[c(1,2,5,8,9)],
           RColorBrewer::brewer.pal(9, "Set1")[c(9,7)],
           RColorBrewer::brewer.pal(8, "Dark2")[c(3,4,8)])

color_profile <- list("type" = "per_row",
                      "col" = cols1)

library(grid)
grid.newpage()
layout.rows <- 2
layout.cols <- 1
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
logo_pssm(pssm,xlab = 'position',
          color_profile = color_profile,
          frame_width = 1,
          pop_name = "",
          newpage = FALSE,
          control = list(quant = 0, round_off = 0,
                         gap_ylab = 3, posbins = 2, negbins =3,
                         viewport.margin.bottom = 3,
                         viewport.margin.left = 5,
                         viewport.margin.top = 2,
                         viewport.margin.right = 2),
          xaxis = FALSE,
          col_line_split = "white")


seekViewport(paste0("plotlogo", 2))
logo_pssm(pssm,xlab = 'position',
          color_profile = color_profile,
          frame_width = 1,
          pop_name = "",
          newpage = FALSE,
          control = list(quant = 0.5, round_off = 0,
                         gap_ylab = 3, posbins = 3, negbins = 3,
                         viewport.margin.bottom = 3,
                         viewport.margin.left = 5,
                         viewport.margin.top = 2,
                         viewport.margin.right = 2),
          xaxis = FALSE,
          col_line_split = "white")


##############  save the figure as 7 by 10 in in R  #####################


