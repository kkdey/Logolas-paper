#setup
scale1=function(x){return(x/sum(x))}
getpwm=function(pfm,bg=NULL,type='dash'){
  if(is.null(bg)){bg=rep(1/nrow(pfm),nrow(pfm))}
  pfm=as.matrix(pfm)
  if(type=='pse'){
    pseudo=sqrt(dim(pfm)[2])*bg
    pwm=apply((pfm+pseudo),2,scale1)
  }else if(type=='prop'){
    pwm=apply(pfm,2,scale1)
  }else if(type=='dash'){
    pwm=(dash((pfm),optmethod = 'mixEM',mode = bg)$posmean)
  }else{pwm=NULL}
  if(nrow(pwm)==4){rownames(pwm)=c('A','C','G','T')}
  if(nrow(pwm)==20){rownames(pwm)=c('A' ,'R','N','D',   'C' ,  'Q',   'E' ,  'G',   'H' ,  'I',   'L' , 'K'  , 'M' ,  'F',   'P' ,  'S'  , 'T' ,  'W',   'Y',   'V')}
  colnames(pwm)=1:dim(pfm)[2]
  return(pwm)
}

#read data
library(Logolas)
readprotein=function(file,skip=3,nsites,bg=bg,adash,mat){
  nrows=length(readLines(file))-9
  rawdata=read.table(file = file,skip = skip,nrows = nrows)
  if(mat=='pfm'){
    pwm=as.matrix(rawdata[,23:42]/100)
    pfm=round(pwm*nsites)
    colnames(pfm)=c('A' ,'R','N','D',   'C' ,  'Q',   'E' ,  'G',   'H' ,  'I',   'L' , 'K'  , 'M' ,  'F',   'P' ,  'S'  , 'T' ,  'W',   'Y',   'V')
    rownames(pfm)=1:nrow(pfm)
    return(t(pfm))
  }
  if(mat=='pwm'){
    pwm=as.matrix(rawdata[,23:42]/100)
    pfm=round(pwm*nsites)
    if(adash==T){
      pwm=dash(as.matrix(pfm),mode = bg,optmethod = 'mixEM')$posmean
    }
    colnames(pwm)=c('A' ,'R','N','D',   'C' ,  'Q',   'E' ,  'G',   'H' ,  'I',   'L' , 'K'  , 'M' ,  'F',   'P' ,  'S'  , 'T' ,  'W',   'Y',   'V')
    rownames(pwm)=1:nrow(pwm)
    return(t(pwm))
  }

  if(mat=='pssm'){
    pssm=as.matrix(rawdata[,3:22])
    colnames(pssm)=c('A' ,'R','N','D',   'C' ,  'Q',   'E' ,  'G',   'H' ,  'I',   'L' , 'K'  , 'M' ,  'F',   'P' ,  'S'  , 'T' ,  'W',   'Y',   'V')
    rownames(pssm)=1:nrow(pssm)
    return(t(pssm))
  }
  if(mat=='pssm_median'){
    pssm=as.matrix(rawdata[,3:22])
    pssm=apply(pssm,2,function(x){x-median(x)})
    colnames(pssm)=c('A' ,'R','N','D',   'C' ,  'Q',   'E' ,  'G',   'H' ,  'I',   'L' , 'K'  , 'M' ,  'F',   'P' ,  'S'  , 'T' ,  'W',   'Y',   'V')
    rownames(pssm)=1:nrow(pssm)
    return(t(pssm))
  }
}

cols1 <- c(rev(RColorBrewer::brewer.pal(12, "Paired"))[c(3,4,7,8,11,12,5,6,9,10)],
           RColorBrewer::brewer.pal(12, "Set3")[c(1,2,5,8,9)],
           RColorBrewer::brewer.pal(9, "Set1")[c(9,7)],
           RColorBrewer::brewer.pal(8, "Dark2")[c(3,4,8)])
color_profile <- list("type" = "per_row",
                      "col" = cols1)

file='http://caps.ncbs.res.in/cgi-bin/mini/databases/3pfdb/3pfdb_pssm_download.cgi?id=PF00389&data_dir=SDB_folder'

pwm389=readprotein(file,mat = 'pwm',nsites = 319,adash = F)
pwm_mt2 = pwm389[,257:267]

library(grid)
grid.newpage()
layout.rows <- 3
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


library(Logolas)
#log odds
seekViewport(paste0("plotlogo", 1))
nlogomaker(pwm_mt2,logoheight = 'log_odds',color_profile = color_profile,
           newpage = FALSE, pop_name = "PF00389 motif 2 (log odds)")
#ratio
seekViewport(paste0("plotlogo", 2))
nlogomaker(pwm_mt2,logoheight = 'ratio',color_profile = color_profile,
           newpage = FALSE, pop_name = "PF00389 motif 2 (ratio)")
#ic log
seekViewport(paste0("plotlogo", 3))
nlogomaker(pwm_mt2,logoheight = 'ic_log',color_profile = color_profile,
           newpage = FALSE, pop_name = "PF00389 motif 2 (ic log)")
#ic log
seekViewport(paste0("plotlogo", 4))
nlogomaker(pwm_mt2,logoheight = 'ic_ratio',color_profile = color_profile,
           newpage = FALSE, pop_name = "PF00389 motif 2 (ic ratio)")
# probKL
seekViewport(paste0("plotlogo", 5))
nlogomaker(pwm_mt2,logoheight = 'ic_log_odds',color_profile = color_profile,
           newpage = FALSE, pop_name = "PF00389 motif 2 (ic log odds)")
# weighted KL
seekViewport(paste0("plotlogo", 6))
nlogomaker(pwm_mt2,logoheight = 'probKL',color_profile = color_profile,
           control = list(quant=0.5), pop_name = "PF00389 motif 2 (probKL)",
           newpage = FALSE)

#########  save the plot as pdf 20 by 15   ######################



