####### pcaplot with the same layout and themes for all future plots
### Alida Kindt-Dunjko
### 14 June 2018

library(ggplot2)
library(ggfortify)

pcaplot <- function(data1, group, title = "", numeric = F, legend.title = ""){
  
  if(numeric == F){
    group = as.character(group)
  }else{
    group = as.numeric(group)
  }
  df1 = cbind(data1, group)
  autoplot(prcomp(scale(data1)), data  = df1, colour = 'group') + ggtitle(title) +
    guides(colour=guide_legend(title=legend.title)) +
    theme_bw()
  
}
