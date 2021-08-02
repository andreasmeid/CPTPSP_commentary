require(ggplot2)
require(Rmisc)

pdp_plot <- function(v.obj=rf.obj, 
                       grid=1:8,  ## Grid, auf dem der Plot fuer eine bestimmte Variable ausgewertet werden soll
                       variable="METSITES", work_data=dat2,
                     x_label=NULL, y_label=NULL, box_labels=NULL) {
  

  pdp <- data.frame(pred=numeric(), grid=numeric())
  
  ## Jetzt einmal ueber das gesamte festegelegte grid gehen
  for(i in 1:length(grid)){
    
    ## Die Zielvariable hier einsetzen und mit den Werten in jeder Iteration ueberschreiben
    work_data[,variable] <- grid[i]
    
    if ("rfsrc" %in% class(v.obj)) {
      predictions <- predict.rfsrc(v.obj, newdata=work_data)$predicted
    } else {
      predictions <- predict(v.obj, newdata=work_data)
    }

      pdp <- rbind(pdp, data.frame(pred=predictions, grid=rep(grid[i], length(predictions))))
  }
  
  if (range(grid)[2] == 1) { ## 0/1-Grid --> binary indicator --> boxplot
    
    pd_plot <- pdp %>% #filter(grid <= 2) %>%
      ggplot(aes(x=factor(grid), y=pred)) + geom_boxplot() +
      scale_x_discrete(name=x_label, labels=box_labels) +
      scale_y_continuous(name=y_label, limits=c(0,400), expand=c(0,0)) + 
      theme_minimal() +
      theme(legend.position = "none",
            axis.line = element_line(size = 1, colour = "black", linetype = "solid"),
            axis.ticks.x = element_line(size = 0.8),
            axis.ticks.y = element_line(size = 0.8),
            axis.ticks.length = unit(.1, "cm"),
            axis.text.x = element_text(colour="black", size=14,angle=0, vjust=0, hjust=0.5),
            axis.text.y = element_text(colour="black", size=12),
            axis.title.y =  element_text(colour="black", size=16, face="bold"),#
            axis.title.x =  element_text(colour="black", size=16, face="bold"),#
            strip.text= element_text(colour="black", size=14),
            plot.title = element_text(hjust = 0.25, size=18),
            panel.spacing = unit(0.2, "lines"),
            plot.margin = unit(c(5,10,5,5), "pt"), #t,r,b,l
            strip.background = element_blank(),
            strip.placement = "outside")
    
  } else { ## continuous Grid --> scatter-plot

    pd_plot <- pdp %>% summarySE(data=., measurevar="pred", groupvars="grid", na.rm=FALSE, conf.interval=.95) %>%
      ggplot(aes(x=grid, y=pred, group=1)) + geom_line() +
      geom_errorbar(width=.1, aes(ymin=pred-ci, ymax=pred+ci)) +
      geom_point(shape=21, size=3, fill="white") +
      scale_x_continuous(name=x_label, expand=c(0,0)) +
      scale_y_continuous(name=y_label, limits=c(0,400), expand=c(0,0)) + 
      theme_minimal() +
      theme(legend.position = "none",
            axis.line = element_line(size = 1, colour = "black", linetype = "solid"),
            axis.ticks.x = element_line(size = 0.8),
            axis.ticks.y = element_line(size = 0.8),
            axis.ticks.length = unit(.1, "cm"),
            axis.text.x = element_text(colour="black", size=14,angle=0, vjust=0, hjust=0.5),
            axis.text.y = element_text(colour="black", size=12),
            axis.title.y =  element_text(colour="black", size=16, face="bold"),#
            axis.title.x =  element_text(colour="black", size=16, face="bold"),#
            strip.text= element_text(colour="black", size=14),
            plot.title = element_text(hjust = 0.25, size=18),
            panel.spacing = unit(0.2, "lines"),
            plot.margin = unit(c(5,10,5,5), "pt"), #t,r,b,l
            strip.background = element_blank(),
            strip.placement = "outside")
    
  }
  return(pd_plot)
}
