############################################
############################################
## SETTING PARAMETERS:
.libPaths("/nas/longleaf/home/jsgomez/github/Project-2/ExperimentsCluster4/req_lib")
library(readr)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(stringr)
library(lazyeval)
library(pROC)

plot_version <- "1"

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
subfolder_new        <- paste0("500_AggregatedDataExperiments/")
subfolder_data_new   <- paste0(subfolder_new, "data_all/")
subfolder_plots_new  <- paste0(subfolder_new, "plots_all/")

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
###################### Generating plots: DEGREE
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


         

for(d_shift in c(2, 5)) {
  for(p_val in c(100, 200, 500)) {
    ## Set simulation parameters:
    sim_ind_load    <- which(diagonal_shift == d_shift & p == p_val)
    type            <- "all"
    results_dir     <- paste0(subfolder_new, "plots_", type, "/")
    
    ## Load data corresponding to specified parameters:
    for (sim_ind in sim_ind_load) {
      load(paste0(
          subfolder_new, "data_all/",
          "output", sim_ind, ".RData"))
    }

    ## Turn into long dataframe to facilitate plotting
    ## from   var1,var2,var3,...,varP ->  var  im
    ##        im1, im2, im3,..., imP      1    im1
    ##        im1, im2, im3,..., imP      2    im2
    ##                ...                   ...
    ## I also eliminate all exactly zero values of im.
    var_names <- paste0("var", 1:p_val)
    
    ##############################
    ##############################
    ## TRUE POSITIVE RATE OF JIC-HD METHODS:

    ## Merge dataset of JIC-HD-derived data.
    output_merged_jic  <- distinct(bind_rows(mget(ls(pattern = '^output\\d+')))) %>%
      filter(METHOD %in% method_names[-c(1,2)]) %>%
      mutate(METHOD = str_replace_all(METHOD, setNames(method_names_clean, method_names))) %>%
      arrange(TASK_ID, SIM_NUM, K_MAT_NUM, METHOD) %>%
      select(-K_MAT_NUM, -TIME)
    output_merged_jic$microrun <- rep(1:10, nrow(output_merged_jic) / 10) 
    output_merged_jic <- output_merged_jic %>% 
      relocate(microrun, 1)
    dim(output_merged_jic)
    head(output_merged_jic[,1:15], 15)

    ## For each row, we calculate the TPR/FPR:
    mat_jic <- t(apply(
      output_merged_jic, MARGIN = 1, 
      function(x) {
        nhubs     <- 5
        p_val     <- length(x) - 9
        id_task   <- x[2]
        args_temp <- get(gsub(" ", "", paste0("args", id_task, sep = "")))
        trueHubs  <- (1:p_val) %in% c(args_temp$Hjoint)
        nhubs     <- length(args_temp$Hjoint)

        vals      <- as.numeric(x[-(1:9)])
        vals_pos  <- vals
        
        tr_mean   <- mean(vals_pos)
        tr_sd     <- sd(vals_pos)

        hubshat   <- vals > tr_mean + 2 * tr_sd # 2.32 * tr_sd

        tp <- sum(hubshat & trueHubs) / (nhubs)
        fp <- sum(hubshat & !trueHubs) / (p_val - nhubs)
        
        return(c(tp, fp))
      }
    ))
    output_merged_jic$tp <- mat_jic[,1]
    output_merged_jic$fp <- mat_jic[,2]
    output_merged_jic <- output_merged_jic %>%
      dplyr::select(!starts_with("var"))

    dim(output_merged_jic)
    colnames(output_merged_jic)
    head(output_merged_jic)

    ##############################
    ##############################

    output_merged_gl   <- distinct(bind_rows(mget(ls(pattern = '^output\\d+')))) %>%
      filter(METHOD %in% method_names[c(1,2)]) %>%
      filter(K_MAT_NUM != 0) %>%
      mutate(METHOD = str_replace_all(METHOD, setNames(method_names_clean, method_names))) %>%
      arrange(TASK_ID, SIM_NUM, K_MAT_NUM, METHOD)
    output_merged_gl$microrun <- rep(1:10, nrow(output_merged_gl) / 10) ## adding micro-run identifier...
    output_merged_gl <- output_merged_gl %>% 
      relocate(microrun, 1)
    dim(output_merged_gl)
    head(output_merged_gl[, 1:15], 10)

    ## Calculate hubs for each of the simulations ran
    hubsdata_gl <- t(apply(
      output_merged_gl, MARGIN = 1, 
      function(x) {
        nhubs     <- 5
        p_val     <- length(x) - 11
        id_task   <- x[2]
        args_temp <- get(gsub(" ", "", paste0("args", id_task, sep = "")))
        trueHubs  <- (1:p_val) %in% c(args_temp$Hjoint)
        nhubs     <- length(args_temp$Hjoint)

        vals      <- as.numeric(x[-c(1:11)])
        vals_pos  <- vals
        
        tr_mean   <- mean(vals_pos)
        tr_sd     <- sd(vals_pos)

        hubshat   <- (vals > tr_mean + 2 * tr_sd) # 2.32 * tr_sd
        return(hubshat)

      }
    ))
    
    ## Aggregate hubs of K_MAT_NUM = 1,2,3 to obtain common hub estimation rate.
    colnames(hubsdata_gl) <- paste0("ishub", 1:p_val)
    output_merged_gl <- cbind(output_merged_gl, hubsdata_gl) %>%
      dplyr::select(!starts_with("var")) %>%
      group_by(microrun, TASK_ID, SIM_NUM, p, T0, n, ph1, ph2, METHOD) %>%
      summarise_at(vars(starts_with("ishub")),
        function(x) {1 * (sum(x) == 3)}) %>%
      ungroup()
    dim(output_merged_gl)
    head(output_merged_gl, 20)

    mat_gl <- t(apply(
      output_merged_gl, MARGIN = 1, 
      function(x) {
        nhubs     <- 5
        p_val     <- length(x) - 9
        id_task   <- x[2]
        args_temp <- get(gsub(" ", "", paste0("args", id_task, sep = "")))
        trueHubs  <- (1:p_val) %in% c(args_temp$Hjoint)
        nhubs     <- length(args_temp$Hjoint)

        hubshat   <- as.numeric(x[-(1:9)])

        tp <- sum(hubshat & trueHubs) / (nhubs)
        fp <- sum(hubshat & !trueHubs) / (p_val - nhubs)
        
        return(c(tp, fp))
      }
    ))
    output_merged_gl$tp <- mat_gl[,1]
    output_merged_gl$fp <- mat_gl[,2]
    output_merged_gl <- output_merged_gl %>%
      dplyr::select(!starts_with("ishub"))


    ##############################
    ##############################
    output_merged <- rbind(output_merged_gl, output_merged_jic)    

    output_summarised <- output_merged %>%
      dplyr::select(!starts_with("var")) %>%
      pivot_longer(
        cols          = c("tp", "fp"),
        names_to      = "eval_par",
        values_to     = "eval") %>%
      
      group_by(METHOD, p, ph1 , ph2, T0, n, eval_par) %>%
      summarise(mean = mean(eval), sd = sd(eval))
    
    file_name <- paste0(
      results_dir, "v", plot_version,
      "_517_TPR_trimmed_deg_p", p_val,"_d", d_shift,".pdf")

    pdf(file_name, width = 7.5, height = 5)
    ## Plot 1: TPR for ph = 0.25
    p1 <-  output_summarised %>% filter(eval_par == "tp", ph1 == 0.25) %>%
      ggplot(aes(x = n, y = mean)) + 
      geom_line(aes(col = METHOD, linetype = METHOD), linewidth = 1) + 
      #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = METHOD), alpha = 0.1) +
      geom_hline(yintercept = c(0,1), linetype = 2) +
      facet_grid(rows = vars(ph2), cols = vars(T0))
    print(p1)
  
    ## Plot 3: TPR for ph = 0.5
    p3 <-  output_summarised %>% filter(eval_par == "tp", ph1 == 0.5) %>%
      ggplot(aes(x = n, y = mean)) + 
      geom_line(aes(col = METHOD, linetype = METHOD), linewidth = 1) + 
      #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = METHOD), alpha = 0.1) +
      geom_hline(yintercept = c(0,1), linetype = 2) +
      facet_grid(rows = vars(ph2), cols = vars(T0))
    print(p3)
    dev.off()


    file_name <- paste0(
      results_dir, "v", plot_version,
      "_517_FPR_trimmed_deg_p", p_val,"_d", d_shift,".pdf")
    pdf(file_name, width = 7.5, height = 5)
    ## Plot 2: FPR for ph = 0.25
    p2 <- output_summarised %>% filter(eval_par == "fp", ph1 == 0.25) %>% 
      ggplot(aes(x = n, y = mean)) +  
      geom_line(aes(col = METHOD, linetype = METHOD), linewidth = 1) + 
      #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = method), alpha = 0.1) +
      geom_hline(yintercept = c(0,1), linetype = 2) +
      facet_grid(rows = vars(ph2), cols = vars(T0))
    print(p2)
    
    ## Plot 4: FPR for ph = 0.5
    p4 <- output_summarised %>% filter(eval_par == "fp", ph1 == 0.5) %>% 
      ggplot(aes(x = n, y = mean)) +  
      geom_line(aes(col = METHOD, linetype = METHOD), linewidth = 1) + 
      #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = method), alpha = 0.1) +
      geom_hline(yintercept = c(0,1), linetype = 2) +
      facet_grid(rows = vars(ph2), cols = vars(T0))
    print(p4)
    dev.off()

    rm(list = grep("output", ls(), value = TRUE))
    rm(list = grep("args", ls(), value = TRUE))
    rm(p1, p2, p3, p4)
  }
}

