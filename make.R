# environment, packages & functions ---------------------------------------

  # set path
    #rm(list=ls())
    #dateipfad <- rstudioapi::getSourceEditorContext()$path
    library(tidyverse, quietly = TRUE)
    getCurrentFileLocation <-  function()
    {
      this_file <- commandArgs() %>% 
        tibble::enframe(name = NULL) %>%
        tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
        dplyr::filter(key == "--file") %>%
        dplyr::pull(value)
      if (length(this_file)==0)
      {
        this_file <- rstudioapi::getSourceEditorContext()$path
      }
      return(dirname(this_file))
    }
    hauptpfad <- hauptpfad_new <- paste0(getCurrentFileLocation(), "/")
    setwd(hauptpfad)
  
  #library(plyr); 
    library(dplyr); 
    library(stringr); library(magrittr)
    library(data.table); library(tables);
    library(purrr); library(furrr)
    library(survival); library(glmnet); library(plotmo);
    library(randomForestSRC); library(pec); library(CoxBoost)
    library(scales)
    library(ggplot2)
    suppressPackageStartupMessages(library(cowplot))
    library(rattle)
    library(SubgrPlots)
    library(ggdendro)
    require(rpart)
    
# recapture objects -------------------------------------------------------

  # data
    yvar <- rf.obj$yvar
    xvar <- rf.obj$xvar
    xvar <- xvar %>% mutate(SEX = SEX - 1, ECOGGR = ECOGGR - 1)
    dat2 <- cbind(yvar,xvar)
    cov.vec <- names(xvar)
  
                   
  # simple treatment effect
    global_effect_model <- coxph(Surv(time, status) ~ ., data=cbind(yvar,xvar))
  
  # re-fit forest
    set.seed(1)
    rf.obj <- rfsrc(as.formula(paste("Surv(time,status) ~ ", paste(cov.vec,collapse = "+"),sep = "")),
                  dat2, nsplit = 10, mtry=6,ntree = 1000,tree.err=TRUE)
    #var.select(rf.obj, importance="permute")
    
  # re-fit lasso
    
    set.seed(123)
    cox.lasso.cv.fit <- cv.glmnet(as.matrix(dat2[,match(cov.vec, colnames(dat2))]), as.matrix(dat2[,1:2]), family="cox",
                                  standardize=T, type.measure = "deviance")
    #coef(cox.lasso.cv.fit, s=c(cox.lasso.cv.fit$lambda.1se))
    cox.lasso.fit <- glmnet(as.matrix(dat2[,match(cov.vec, colnames(dat2))]), as.matrix(dat2[,1:2]), family="cox",
                                  standardize=T)
    
    par(mfrow=c(1,1)); par(mar=c(4.2,3.1,3.1,2))
    plotmo::plot_glmnet(cox.lasso.fit, label=T)
    abline(v=log(cox.lasso.cv.fit$lambda.1se), lty=4)
    
    cox.lasso.fit <- glmnet(as.matrix(dat2[,match(cov.vec, colnames(dat2))]), as.matrix(dat2[,1:2]), family="cox",
                            standardize=T, lambda = cox.lasso.cv.fit$lambda.1se)
    
    lasso_coefs <- as.numeric(coef(cox.lasso.cv.fit, s=c(cox.lasso.cv.fit$lambda.1se)))[which(as.numeric(coef(cox.lasso.cv.fit, s=c(cox.lasso.cv.fit$lambda.1se))) != 0)]
    lasso_coefnames <- cov.vec[which(as.numeric(coef(cox.lasso.cv.fit, s=c(cox.lasso.cv.fit$lambda.1se))) != 0)]
    lasso_formula <- as.formula(paste0("Surv(time, status) ~ ", paste(lasso_coefnames, collapse = " + ")))
    
    lasso_model_cox <- coxph(lasso_formula, data=dat2)
    lasso_model_cox$coefficients <- lasso_coefs
    
    
  # re-fit boosting
    set.seed(1)
    coxBoost.obj <- CoxBoost(time=dat2$time, status=dat2$status, x=as.matrix(dat2[,match(cov.vec, colnames(dat2))]))
    
    boosting_coefs <- coefficients(coxBoost.obj)[which(coefficients(coxBoost.obj) != 0)]
    boosting_coefnames <- cov.vec[which(coefficients(coxBoost.obj) != 0)]
    boosting_formula <- as.formula(paste0("Surv(time, status) ~ ", paste(boosting_coefnames, collapse = " + ")))
    
    boosting_model_cox <- coxph(boosting_formula, data=dat2)
    boosting_model_cox$coefficients <- boosting_coefs


# variable importance -----------------------------------------------------

    # muesste auf predict-Fkt basieren, damit vergleichbar zwischen Methoden
      
    source(as.character(paste0(hauptpfad, "lib/breiman.R")))
    
    set.seed(123)
    vi_global <- breiman(model=global_effect_model, time=dat2$time, event=dat2$status, covdata=xvar,nrep=250, plot=FALSE)
    vi_boosting <- breiman(model=coxBoost.obj, time=dat2$time, event=dat2$status, covdata=xvar,nrep=250, plot=FALSE)
    vi_lasso <- breiman(model=lasso_model_cox, time=dat2$time, event=dat2$status, covdata=xvar,nrep=250, plot=FALSE)
    vi_rfsrc <- vimp(rf.obj, xvar.names=names(xvar), importance="permute", seed=123)
    #vi_rfsrc <- predict(rf.obj, get.tree=1:250, importance = TRUE)$importance
    
    # --> skalieren auf jeweils 100 % der wichtigsten Variable im jeweiligen Modell?
    
    vi_lasso.df <- data.frame(vi_lasso=vi_lasso, modeltype="lasso", predictor=names(vi_lasso)) %>% 
      mutate(vimp=rescale(vi_lasso, to = c(0, 100))) %>% select(predictor, vimp, modeltype)
    vi_boosting.df <- data.frame(vi_boosting=vi_boosting, modeltype="boosting", predictor=names(vi_boosting)) %>% 
      mutate(vimp=rescale(vi_boosting, to = c(0, 100))) %>% select(predictor, vimp, modeltype)
    vi_rfsrc.df <- data.frame(vi_rfsrc=as.vector(vi_rfsrc$importance), modeltype="rfsrc", predictor=names(vi_rfsrc$importance),
                              vimp = rescale(as.vector(vi_rfsrc$importance), to = c(0, 100))) %>% select(predictor, vimp, modeltype)
    
    vi.df <- rbind(vi_lasso.df, vi_boosting.df, vi_rfsrc.df)
    
    vi.df$indikator <- rep(rev(letters[1:21]),3)
    var_labels <- c("Age", "Sex", "Body weight", "ECOG performance status", "Number of prior\nchemotherapy regimens",  "Total protein",
                    "Albumin", "Alkaline phosphatase", "Aspartate aminotransferase", "Lactate dehydrogenase",
                    "Serum creatinine", "Number of metastatitc sites", "Estimated\nglomerular filtration rate",
                    "PD-L1 expression TC123", "PD-L1 expression TC23", "PD-L1 expression TC3",
                    "Smoking", "Histology", "More than two\nmetastatic sites", "More than four\nmetastatic sites", "Treatment")
    
    
    vi_plot <- vi.df %>% ggplot(aes(x=vimp, y=indikator, group=factor(modeltype), shape=factor(modeltype))) + 
      geom_point(size=4, position = position_dodge((0.8))) +
      scale_x_continuous(name="Variable importance\n(rescaled to 100 percent)") +
      scale_y_discrete(name="", labels = rev(var_labels)) + #Predictor variable
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
    

# pdp ---------------------------------------------------------------------

  source(as.character(paste0(hauptpfad, "lib/pdp_plot.R")))
  
  # top 6 predictors
    names(sort(vi_rfsrc$importance))[(length(vi_rfsrc$importance)-7):length(vi_rfsrc$importance)]
    
  # --> Boxplot fuer binaere Variable?!  
    MET_g2_pdp <- pdp_plot(v.obj=rf.obj, grid=0:1, variable="MET_g2", work_data=dat2,
                           x_label="Number of metastatic sites", y_label=NULL, box_labels=c("<= 2", "> 2"))
    METSITES_pdp <- pdp_plot(v.obj=rf.obj, grid=1:8, variable="METSITES", work_data=dat2,
                           x_label="Number of metastatic sites", y_label="predicted risk")
    ALB_pdp <- pdp_plot(v.obj=rf.obj, grid=seq(0,55, by=5), variable="ALB", work_data=dat2,
                             x_label="Albumin", y_label="predicted risk")
    HISTOLOGY_pdp <- pdp_plot(v.obj=rf.obj, grid=0:1, variable="HISTOLOGY", work_data=dat2,
                           x_label="Histology", y_label=NULL, box_labels=c("non-squamous", "squamous"))
    LDH_pdp <- pdp_plot(v.obj=rf.obj, grid=seq(0,3250, by=250), variable="LDH", work_data=dat2,
                        x_label="Lactate dehydrogenase\n", y_label="predicted risk")
    ECOGGR_pdp <- pdp_plot(v.obj=rf.obj, grid=0:1, variable="ECOGGR", work_data=dat2,
                           x_label="ECOG group\n", y_label=NULL, box_labels=c("1", "2"))
    AST_pdp <- pdp_plot(v.obj=rf.obj, grid=seq(0,100, by=10), variable="AST", work_data=dat2,
                           x_label="Aspartate aminostransferase", y_label="predicted risk")
    ARM_pdp <- pdp_plot(v.obj=rf.obj, grid=0:1, variable="ARM", work_data=dat2,
                        x_label="Treatment", y_label=NULL, box_labels=c("Docetaxel", "Atezolizumab"))
    
    pd_plots <- plot_grid(AST_pdp, ARM_pdp,  METSITES_pdp, MET_g2_pdp, ALB_pdp,HISTOLOGY_pdp,LDH_pdp,ECOGGR_pdp,
                          #labels = c('A', 'B','C', 'D','E', 'F'),
                          ncol = 2)
    
    Fig1 <- plot_grid(vi_plot, pd_plots, ncol=2, labels=c("A", "B"), rel_widths = c(2.2,3))
    
    
# split model -------------------------------------------------------------

  cox_treated <- dat2 %>% filter(ARM == 1) %>% select(-ARM) %>% coxph(Surv(time, status) ~ ., data=.)
  cox_control <- dat2 %>% filter(ARM == 0) %>% select(-ARM) %>% coxph(Surv(time, status) ~ ., data=.)
    
  cox_pred1 <- predict(cox_treated, newdata = dat2)
  cox_pred0 <- predict(cox_control, newdata = dat2)
  
  cate <- cox_pred1 - cox_pred0

  rf_treated <- dat2 %>% filter(ARM == 1) %>% select(-ARM) %>% rfsrc(Surv(time,status) ~ ., data=., nsplit = 10, mtry=6,ntree = 1000,tree.err=TRUE)
  rf_control <- dat2 %>% filter(ARM == 0) %>% select(-ARM) %>% rfsrc(Surv(time,status) ~ ., data=., nsplit = 10, mtry=6,ntree = 1000,tree.err=TRUE)
  
  rf_pred1 <- predict(rf_treated, newdata = dat2)$predicted
  rf_pred0 <- predict(rf_control, newdata = dat2)$predicted
  
  cate <- rf_cate <- rf_pred1 - rf_pred0
  
# waterfall plot ----------------------------------------------------------

  waterfall_plot <- data.frame(y=sort(cate), #lwr=post_contrasts[[1]][,2], upr=post_contrasts[[1]][,3], 
                   x=seq(0,100,length.out=length(cate))) %>% 
    ggplot(aes(x=x, y=y)) + geom_hline(yintercept=0, colour="gray") +
    geom_point() +  #geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.3) +
    scale_x_continuous(name='Ordered Patients', breaks=seq(0,100,20), labels=paste0(seq(0,100,20),'%'), limits=c(0,102), expand=c(0,0)) +
    scale_y_continuous(name='Predicted risk difference', expand=c(0,0)) +
    theme_minimal() +
    theme(axis.line = element_line(size = 1.2, colour = "black", linetype = "solid"),
          axis.ticks = element_line(size = 1.2),
          axis.ticks.length = unit(.5, "cm"),
          axis.text.x = element_text(colour="black", size=14, angle=0, vjust=0.45),
          axis.text.y = element_text(colour="black", size=16),
          axis.title.y = element_text(hjust = 0.95, size=18),
          axis.title.x = element_text(hjust = 0.5, size=18),
          plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),
          #legend.title=element_blank(),
          legend.text = element_text(colour=colors, size = 14, face = "bold"),
          legend.position = "bottom",#,
          legend.key.height=unit(2, "cm"))
  

# fit-the-fit -------------------------------------------------------------

  cart_data_cox <-  dat2 %>% bind_cols(data.frame(y=cate)) 
 
  cart_formula <- as.formula(paste("y ~ ", paste0(names(xvar), collapse=" +"), sep=""))
  cart_formula <- as.formula(paste("y ~ AGE + factor(SEX) + BWT + factor(ECOGGR) + PRIORTXC + PROT + ALB + ALP + 
    AST + LDH + CREAT + METSITES + eGFR + factor(TC123) + factor(TC23) + factor(TC3) + 
    factor(SMOKING) + factor(HISTOLOGY) + MET_g2 + MET_g4"))
  
  # grow tree
    fit_cart_cox <- rpart(cart_formula, method="anova", data=cart_data_cox)
    #printcp(fit_cart_cox)
    #plotcp(fit_cart_cox)
    #summary(fit_cart_cox)
  
  # create additional plots
    #par(mfrow=c(1,2)) # two plots on one page
    #rsq.rpart(fit_cart_cox) # visualize cross-validation results  
  
  # plot (full) tree
    #plot(fit_cart_cox, uniform=TRUE,main="Regression Tree for Mileage ")
    #text(fit_cart_cox, use.n=TRUE, all=TRUE, cex=.8)
  
  # prune the tree
    pfit <- prune(fit_cart_cox, cp=0.03)
    #plot(pfit, uniform=TRUE, main="Pruned Regression Tree (Fit-the-Fit)")
    #text(pfit, use.n=TRUE, all=TRUE, cex=.8)
 
    fitr <- dendro_data(pfit)
    tree <- ggplot() +
      geom_segment(data = fitr$segments, 
                   aes(x = x, y = y+0.025, xend = xend, yend = yend+0.025)
      ) +
      geom_text(data = fitr$labels, aes(x = x, y = y, label = c("ECOG group > 1\n(yes / no)",
                                                                "TC3 positive\n(yes / no)",
                                                                "LDH > 459\n(yes / no)",
                                                                "Metastatic sites > 4\n(yes / no)",
                                                                "TC3 positive\n(yes / no)"
                                                                )), size=5) +
      geom_text(data = fitr$leaf_labels, aes(x = x, y = y, label = label), colour=c("black", "black",
                                                                                   "darkgray", "black",
                                                                                   "black", "darkgray"), size=5) +
      theme_dendro() + theme(legend.position = "none")
    
    
   
    top <- plot_grid(waterfall_plot, tree, labels = c('A', 'B'),ncol = 2)


# Upset plot --------------------------------------------------------------
    
    dat2$LDH_group <- ifelse(dat2$LDH<459, 0, 1)
    
    Upset_plot <- dat2 %>% mutate(trt=factor(ifelse(ARM == 1, "Exp_Atezolizumab", "Cont_Docetaxel"))) %>% 
      subgroupset(.,
                  order.by = "freq", #degree
                  empty.intersections = "on",
                  sets = c("MET_g4", "TC3","LDH_group", "ECOGGR"),#, "LDH_group", "ECOGGR"
                  text.scale = 1.5,
                  mb.ratio = c( 0.25, 0.45, 0.30),
                  treatment.var = "trt",
                  outcome.type = "survival",
                  effects.summary = c("time", "status"),
                  query.legend = "bottom", icon = "pm", transpose = T,
                  cutoff=6, color.pal="Greys", 
                  mainbar.y.label="Intersection size",
                  #sets.bar.color= c("gray", "black"),
                  nsets = 12,
                  min.n = 80
                  )  
   
   
