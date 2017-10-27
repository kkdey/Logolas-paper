library(Logolas)
library(grid)
library(dplyr)
###############################
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

grid.newpage()
layout.rows <- 4
layout.cols <- 6
top.vp <- viewport(layout=grid.layout(layout.rows, layout.cols,
                                      widths=unit(rep(5,layout.cols), rep("null", 2)),
                                      heights=unit(rep(5,layout.rows), rep("null", 1))))
plot_reg <- vpList()
l <- 1
for(i in 1:layout.rows){
  for(j in 1:layout.cols){
    plot_reg[[l]] <- viewport(layout.pos.col = j, layout.pos.row = i, name = paste0("plotlogo", l))
    l <- l+1
  }
}




files=c('ALL.2','Bladder.2','Breast.2','Cervix.2','Head-and-Neck.3',
        'Kidney-Papillary.5','Lung-Adeno.5',
        'Lung-Squamous.3','Lymphoma-B-cell.4',
        'Myeloma.4','Pancreas.4','Thyroid.4',
        'Head-and-Neck.4','Lung-Adeno.3','Lung-Squamous.4',
        'Lung-Small-Cell.2','Colorectum.6','Uterus.4','Colorectum.5','Uterus.6',
        'Head-and-Neck.6','Melanoma.2','Lung-Small-Cell.6',
        'Stomach.5'
)
names=c('ALL_2','Bladder_2','Breast_1','Cervix_1','Head-and-Neck_2',
        'Kidney-Papillary_2','Lung-Adeno_2',
        'Lung-Squamous_1','Lymphoma-B-cell_2',
        'Myeloma_2','Pancreas_2','Thyroid_1',
        'Head-and-Neck_1','Lung-Adeno_3','Lung-Squamous_2',
        'Lung-Small-Cell_1','Colorectum_4','Uterus_1','Colorectum_3','Uterus_2',
        'Head-and-Neck_3','Melanoma_3','Lung-Small-Cell_2',
        'Stomach_1')
idx=c(1,1,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,1,3,1,1,1,1,1)


plot_tree <- vpTree(top.vp, plot_reg)

pushViewport(plot_tree)

for(j in 1:length(files)){
  file=files[j]
  file=paste(file,'.Rdata',sep = '')
  name=names[j]
  i=idx[j]
  load(paste0("../data/Param_ind5_dir/",file))

  mat=resultForSave[[1]]@signatureFeatureDistribution[i,,]
  mat1=cbind(t(mat[2:3,1:4]),rep(NA,4),t(mat[4:5,1:4]))
  rownames(mat1)=c('A','C','G','T')
  colnames(mat1) = c("-2", "-1", "0", "1", "2")
  mat2=cbind(rep(NA,6),rep(NA,6),(mat[1,]),rep(NA,6),rep(NA,6))
  colnames(mat2) = c("-2", "-1", "0", "1", "2")
  rownames(mat2) = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  table = rbind(mat1, mat2)


  seekViewport(paste0("plotlogo", j))
  Logolas::nlogomaker(table,
             logoheight = 'log',
             color_profile = color_profile,
             frame_width = 1,
             xlab = "Position",
             newpage = F,
             pop_name = paste(name),
             control = list(epsilon=0.2,gap_ylab=3.5)
  )
}
