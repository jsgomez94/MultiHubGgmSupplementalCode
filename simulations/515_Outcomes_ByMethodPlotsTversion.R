############################################
############################################
## SETTING PARAMETERS:
wd <- getwd()
.libPaths(paste0(wd,"/req_lib"))

source("001_requirements.R")
library(readr)
library(plot.matrix)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(stringr)
library(lazyeval)
library(pROC)
library(ggh4x)
library(tidyverse)
plot_version <- "1"

source("061_Simulation_CreatingParameters.R")
args <- CreateParameters(1)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


###################### Parameter table:
runtype       <- 2 # FOR EXPERIMENT RUNS
#runtype       <- 3 # FOR FULL RUNS
index_old     <- 1 # run index to use
sim_par_table <- expand.grid(
  K              = 3,
  hmin1          = 4, hmax1         = 5,
  hmin2          = 4, hmax2         = 5,
  nhmin          = 4, nhmax         = 5,
  neffmin        = 4, neffmax       = 5,
  shuffle       = FALSE,
  type          = "unif",
  running_days  = ifelse(runtype <= 2, 1, 5),
  threshold     = 2,
    
  r2              = c(3),
  r1              = c(5),
  pneff           = c(0.01),
  pnh             = c(0.05),
  ph2             = c(0.05, 0.25, 0.5),
  ph1             = c(0.25, 0.5),
    
  nsim            = ifelse(runtype <= 2, 2, 5),
  diagonal_shift  = c(2,5),
  #n_prop          = c(0.25, 0.5, 0.75, 1),
  n_prop          = c(0.5, 0.75, 1, 1.25),
  T0_prop         = c(0.5, 0.75, 1),
  p               = c(100, 200, 500))
attach(sim_par_table)

###################### Creating folders:
subfolder_new        <- "500_AggregatedDataExperiments/"
subfolder_data_new   <- paste0(subfolder_new, "data_all/")
subfolder_plots_new  <- paste0(subfolder_new, "plots_bymethod/")

if (!dir.exists(subfolder_new)) {
  dir.create(subfolder_new)
}
if (!dir.exists(subfolder_data_new)) {
  dir.create(subfolder_data_new)
}
if (!dir.exists(subfolder_plots_new)) {
  dir.create(subfolder_plots_new)
}


##################################################################
##################################################################
###################### Generating plots: WLINE
##################################################################
##################################################################

## Choosing data to load
method_names <- c(
    "GL.CORR.d",          ## GLASSO-methods
    "HWGL.CORR.d",        ## HWGL-methdos.
    "ST.OVER.CORR.IM", "ST.ORAC.CORR.IM", 
    "ST.OVER.THR.IM", "ST.ORAC.THR.IM")  ## IPC-HD
method_names_clean <- c(
    "GLASSO",          ## GLASSO-methods
    "HWGL",        ## HWGL-methdos.
    "SampStiefelOver", "SampStiefelOrac.",
    "ThrStiefelOver", "ThrStiefelOrac.")  ## IPC-HD

rm(list = grep("output", ls(), value = TRUE))
for (p_val in c(100, 200)) {
  for (T0_prop_val in c(0.5, 0.75, 1)) {
    for (d_shift in c(2,5)) {
      
      ## Set simulation parameters:
      args <- CreateParameters(1)
      T0_val          <- p_val * T0_prop_val
      sim_ind_load    <- which(diagonal_shift == d_shift & p == p_val & T0_prop == T0_prop_val)
      type            <- "bymethod"
      results_dir     <- paste0(subfolder_new, "plots_", type, "/")


      ## Load data corresponding to specified parameters:
      for (sim_ind in sim_ind_load) {
        load(paste0(
          subfolder_new, "data_all/",
          "output", sim_ind, ".RData"))
      }
      
      ## Merging all datasets:
      output_merged   <- distinct(bind_rows(mget(ls(pattern = '^output\\d+')))) %>%
        filter(METHOD %in% method_names) %>%
        mutate(METHOD = str_replace_all(METHOD, setNames(method_names_clean, method_names))) %>%
        mutate(METHOD = str_replace_all(METHOD, setNames("HWGL","HWGLASSO"))) %>%
        mutate(METHOD = factor(METHOD, levels = method_names_clean))

      ## Processing data:
      output_summarized <- output_merged %>% 
        pivot_longer(
          cols = starts_with("var"),
          names_to = "var_ind",
          names_prefix = "var",
          values_to = "var_deg") %>%

        group_by(p, ph1, ph2, T0, n, var_ind, METHOD, K_MAT_NUM) %>%
        
        summarize(mean = mean(var_deg), sd = sd(var_deg)) %>%
        
        mutate(
          var_ind = as.numeric(var_ind),
          GLASSO_Method = ifelse(METHOD %in% method_names_clean[1:2], "GLASSO-Based", "Non-GL-Based"),
          hub = ifelse(var_ind %in% args$Hjoint, 1, 0)) %>%
        
        ungroup()

      ## Plot 1: TPR
      file_name <- paste0(
        results_dir, p_val,"_T", T0_val, "_d", d_shift, "_DegByMethod.pdf")
      pdf(file_name, width = 12, height = 8)      

      p1 <- output_summarized %>%
        filter(METHOD %in% method_names_clean[1:2]) %>%
        mutate(
          ph1_name = ifelse(ph1 == 0.5, "p1[h] == 0.5", "p1[h] == 0.25"),
          ph2_name = ifelse(ph2 == 0.5, "p2[h] == 0.5", ifelse(ph2 == 0.25, "p2[h] == 0.25", "p2[h] == 0.05"))) %>%

        ggplot(aes(x = n, y = mean, group = factor(var_ind):factor(K_MAT_NUM))) +
          geom_line(aes(col = factor(hub)), linewidth = 0.75) + 
          geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = factor(hub)), alpha =  0.1) +
          facet_grid(ph2_name ~ METHOD + ph1_name, labeller = label_parsed, scales='free_y')
      print("It happened 1")
      print(p1)

      p2 <- output_summarized %>%
        filter(METHOD %in% method_names_clean[3:6]) %>%
        mutate(
          ph1_name = ifelse(ph1 == 0.5, "p1[h] == 0.5", "p1[h] == 0.25"),
          ph2_name = ifelse(ph2 == 0.5, "p2[h] == 0.5", ifelse(ph2 == 0.25, "p2[h] == 0.25", "p2[h] == 0.05")),
          METHOD   = factor(METHOD)) %>%
  
        ggplot(aes(x = n, y = mean, group = var_ind)) +
          geom_line(aes(col = factor(hub), group = var_ind), linewidth = 0.75) + 
          geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = factor(hub)), alpha =  0.1) +
          facet_grid(ph2_name ~  METHOD + ph1_name, labeller = label_parsed, scales='free_y')
      print("It happened 2")
      print(p2)
      dev.off()

      
      rm(list = grep("output", ls(), value = TRUE))
      rm(list = grep("args", ls(), value = TRUE))
      print(ls())
    }
  }
}
