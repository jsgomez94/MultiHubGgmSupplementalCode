############################################
############################################
## SETTING PARAMETERS:
.libPaths("/nas/longleaf/home/jsgomez/github/Project-2/ExperimentsCluster4/req_lib")
library(readr)
library(plot.matrix)
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
###################### TPR / FPR CALCULATION FUNCTIONS
##################################################################
##################################################################

## The following functions are aimed at estimating common hubs
## across multiple GGMs, using the individual estimated degree.

## IDEA: estimate hubs forr each population, then intersect. 

## Function that, given a reduced dataframe,
## It calculates the set of common hubs in 
## wide TRUE/FALSE format:
hub_intsersect <- function(df, sd_threshold = 2) {
    p_val  <- ncol(df) - 10
    n_row  <- nrow(df)
    if (nrow(df) > 1) {
        thresholds <- apply(df, MARGIN = 1, 
            function(x){
                p_val     <- sum(str_count(colnames(df), "var") > 0)
                n_pars    <- ncol(df) - p_val
                K         <- x[3]
                id_task   <- x[1]
                args_temp <- get(gsub(" ", "", paste0("args", id_task, sep = "")))
                
                vals      <- as.numeric(x[str_count(colnames(df), "var") > 0])
                vals_pos  <- vals
        
                tr_mean   <- mean(vals_pos, trim = 0)
                tr_sd     <- sd(vals_pos)#sd(vals_pos[vals_pos >= quantile(vals_pos, 0.05) & vals_pos <= quantile(vals_pos, 0.95)])

                return(tr_mean + sd_threshold * tr_sd)
            }
        )
        mat <- matrix(F, ncol = p_val, nrow = nrow(df))
        for (.i in 1:nrow(df)) {
            mat[.i, ] <- df[.i, str_count(colnames(df), "var") > 0] < thresholds[.i]
        }
        common_hubs <- apply(mat, MARGIN = 2, function(x) {sum(x) == 0})
        names(common_hubs) <- paste0("var", 1:(p_val))
        
        output <- as.data.frame(t(common_hubs))
        output <- cbind(data.frame(nrow = n_row), output)
        return(output)
    } else if (nrow(df) == 1) {
        p_val     <- sum(str_count(colnames(df), "var") > 0)
        n_pars    <- ncol(df) - p_val
        K         <- df[1,3]
        id_task   <- df[1,1]
        args_temp <- get(gsub(" ", "", paste0("args", id_task, sep = "")))
            
        vals      <- as.numeric(df[1, str_count(colnames(df), "var") > 0])
        vals_pos  <- vals
        tr_mean   <- mean(vals_pos, trim = 0.05)
        tr_sd     <- sd(vals_pos[vals_pos >= quantile(vals_pos, 0.05) & vals_pos <= quantile(vals_pos, 0.95)])

        common_hubs <-  (df[1, str_count(colnames(df), "var") > 0] > tr_mean + sd_threshold * tr_sd)
        names(common_hubs) <- paste0("var", 1:(p_val))

        output <- as.data.frame(t(common_hubs))
        output <- cbind(data.frame(nrow = n_row), output)
        return(output)
    } else return(NULL)
}
## Given a wide table of TRUE/FALSE format that
## captures common hub estimation, this function
## provides the TPR/FPR summary.
hub_stat_summary <- function(df) {
        p_val     <- sum(str_count(colnames(df), "var") > 0)
        n_pars    <- ncol(df) - p_val

        id_task   <- df[1,1]
        args_temp <- get(gsub(" ", "", paste0("args", id_task, sep = "")))
        trueHubs  <- (1:p_val) %in% c(args_temp$Hjoint)
        nhubs     <- length(args_temp$Hjoint)

        tp <- sum(as.logical(df[1, str_count(colnames(df), "var") > 0]) & trueHubs) / nhubs
        fp <- sum(as.logical(df[1, str_count(colnames(df), "var") > 0]) & !trueHubs) / (p_val - nhubs)
        
        return(data.frame(tp = tp, fp = fp, nrow = nrow(df), ncol = ncol(df)))
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

for (d_shift in c(2, 5)) {
  for (p_val in c(100, 200)) {
    ## Set simulation parameters:
    sd_threshold    <- 2
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
    output_merged   <- distinct(bind_rows(mget(ls(pattern = '^output\\d+')))) %>%
      filter(METHOD %in% method_names) %>%
      mutate(METHOD = str_replace_all(METHOD, setNames(method_names_clean, method_names)))
    var_names <- paste0("var", 1:p_val) #output_merged %>% select() %>% names()

    hub_measures_intsersec <- output_merged %>%
        group_by(TASK_ID, SIM_NUM, p, T0, n, ph1, ph2, METHOD) %>%
        do({
            hub_intsersect(.)
        }) %>%
        group_by(TASK_ID, SIM_NUM, p, T0, n, ph1, ph2, METHOD) %>%
        do({
            hub_stat_summary(.)
        })


    output_summarised <- hub_measures_intsersec %>%
      pivot_longer(
        cols          = c("tp", "fp"),
        names_to      = "eval_par",
        values_to     = "eval") %>%
      group_by(p, T0, n, ph1, ph2, METHOD, eval_par) %>%
      summarise(mean = mean(eval), sd = sd(eval)) %>% 
      mutate(
        ph2_name = ifelse(ph2 == 0.5, "p2[h] == 0.5", 
                                      ifelse(ph2 == 0.25, "p2[h] == 0.25", "p2[h] == 0.05")),
        T0_name = ifelse(T0 == p, paste0("T0 == ", p), 
                                  ifelse(T0 == 0.75 * p, paste0("T0 == ", 0.75 * p), paste0("T0 == ", 0.5 * p))))


    
    file_name <- paste0(
      results_dir, "516_v", plot_version,
      "_trimmed_deg_p", p_val,"_d", d_shift,".pdf")

    pdf(file_name, width = 7.5, height = 5)
    ## Plot 1: TPR for ph = 0.25
    p1 <-  output_summarised %>% filter(eval_par == "tp", ph1 == 0.25) %>%
      ggplot(aes(x = n, y = mean)) + 
      geom_line(aes(col = METHOD, linetype = METHOD), linewidth = 1) + 
      #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = METHOD), alpha = 0.1) +
      geom_hline(yintercept = c(0,1), linetype = 2) +
      facet_grid(rows = vars(ph2_name), cols = vars(T0_name)) +
      ggtitle(paste("TPR, ph1 = 0.25"))
    print(p1)
    
    ## Plot 2: FPR for ph = 0.25
    p2 <- output_summarised %>% filter(eval_par == "fp", ph1 == 0.25) %>% 
      ggplot(aes(x = n, y = mean)) +  
      geom_line(aes(col = METHOD, linetype = METHOD), linewidth = 1) + 
      #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = method), alpha = 0.1) +
      geom_hline(yintercept = c(0,1), linetype = 2) +
      facet_grid(rows = vars(ph2_name), cols = vars(T0_name)) +
      ggtitle(paste("FPR, ph1 = 0.25"))
    print(p2)
  
    ## Plot 3: TPR for ph = 0.5
    p3 <-  output_summarised %>% filter(eval_par == "tp", ph1 == 0.5) %>%
      ggplot(aes(x = n, y = mean)) + 
      geom_line(aes(col = METHOD, linetype = METHOD), linewidth = 1) + 
      #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = METHOD), alpha = 0.1) +
      geom_hline(yintercept = c(0,1), linetype = 2) +
      facet_grid(rows = vars(ph2_name), cols = vars(T0_name)) +
      ggtitle(paste("TPR, ph1 = 0.5"))
    print(p3)

    ## Plot 4: FPR for ph = 0.5
    p4 <- output_summarised %>% filter(eval_par == "fp", ph1 == 0.5) %>% 
      ggplot(aes(x = n, y = mean)) +  
      geom_line(aes(col = METHOD, linetype = METHOD), linewidth = 1) + 
      #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = method), alpha = 0.1) +
      geom_hline(yintercept = c(0,1), linetype = 2) +
      facet_grid(rows = vars(ph2_name), cols = vars(T0_name)) +
      ggtitle(paste("FPR, ph1 = 0.5"))
    print(p4)
  

    dev.off()

    rm(list = grep("output", ls(), value = TRUE))
    rm(list = grep("args", ls(), value = TRUE))
    rm(p1, p2, p3, p4)
  }
}
