library(tidyverse)
library(BGLR)
library(parallel)
library(coda)


# Load the data -----------------------------------------------------------
# Phenotype and weather data
data <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$"))
ped_site <- data %>%
  dplyr::select(Site, PedigreeNew) %>%
  split(., .$Site)

# Get the selected environmental variables
ga <- read_rds("data/weather/ga_select_R.rds")[["TXH1_2016"]]$ga$g[-1]

tibble(V = ga) %>%
  separate(V, c("Variable", "Start", "End"), sep = "_", remove = TRUE) %>%
  mutate_at(c("Start", "End"), str_replace, pattern = "X", replacement = "") %>%
  mutate_at(c("Start", "End"), as.numeric) %>%
  mutate(Ypos = factor(Variable) %>% as.integer(), 
         Ypos = Ypos + c(0.125, -0.125, 0, 0, 0)) %>%
  ggplot(.) + theme_classic() +
    geom_segment(aes(x = Start, xend = End, y = Ypos, yend = Ypos, 
                     colour = Variable), size = 2) +
    geom_vline(xintercept = 1, linetype = 2) +
    labs(x = "% CHU to anthesis", y = "") + guides(colour = "none") +
    scale_colour_manual(values = c("TMIN" = "blue", "TMAX" = "red", "PPT" = 
                                     "skyblue", "NET" = "orange")) +
    scale_y_continuous(breaks = 1:4, labels = c("NET", "PPT", "TMAX", "TMIN")) +
    scale_x_continuous(limits = c(0, 1.5), labels = scales::percent)
ggsave("figures/select/final_variables.pdf", width = 6, height = 4, units = "in", dpi = 300)

data <- data[, c("Site", "PedigreeNew", "BLUE", ga)]
data <- data %>%
  separate(Site, c("Location", "Year"), sep = "_", remove = FALSE)

# Save this data table for future use
write_rds(data, "data/phenotype/yield_blue_final.rds")

# A site-wise table for post hoc analysis of single-site results
data %>%
  group_by(Site) %>%
  summarise_at(c("BLUE", "TMIN_X0.225_X0.45", "TMIN_X0.45_X1.1", 
                 "TMAX_X0.75_X0.95", "PPT_X0.1_X0.225", "NET_X0.425_X1.125"), 
               median) %>%
  ungroup() %>%
  write_rds(., "data/phenotype/yield_site_final.rds")


# Center and scale the predictors -----------------------------------------
means <- apply(data[, 6:10], 2, mean)
sds <- apply(data[, 6:10], 2, sd)

for (i in 6:10) {
  data[[i]] <- (data[[i]] - means[i - 5])/sds[i - 5]
}


# Compute joint reaction norms --------------------------------------------
# Construct a model object
ETA <- list(fixed = list(~ factor(Site) + TMIN_X0.225_X0.45 + TMIN_X0.45_X1.1 + 
                           TMAX_X0.75_X0.95 + PPT_X0.1_X0.225 + NET_X0.425_X1.125, 
                         data = data, model = "FIXED", saveEffects = TRUE), 
            hybrid = list(~ 0 + factor(PedigreeNew), data = data, model = "BRR", 
                          saveEffects = TRUE), 
            TMIN_X0.225_X0.45 = list(~ 0 + factor(PedigreeNew):TMIN_X0.225_X0.45, 
                                     data = data, model = "BRR", saveEffects = TRUE), 
            TMIN_X0.45_X1.1 = list(~ 0 + factor(PedigreeNew):TMIN_X0.45_X1.1, 
                                   data = data, model = "BRR", saveEffects = TRUE), 
            TMAX_X0.75_X0.95 = list(~ 0 + factor(PedigreeNew):TMAX_X0.75_X0.95, 
                                    data = data, model = "BRR", saveEffects = TRUE), 
            PPT_X0.1_X0.225 = list(~ 0 + factor(PedigreeNew):PPT_X0.1_X0.225, 
                                   data = data, model = "BRR", saveEffects = TRUE), 
            NET_X0.425_X1.125 = list(~ 0 + factor(PedigreeNew):NET_X0.425_X1.125, 
                                     data = data, model = "BRR", saveEffects = TRUE))

# Construct groups for heterogeneous error variances
groups <- as.integer(factor(data$Site))

# Run 4 chains
cl <- makeCluster(4)
clusterEvalQ(cl, library(BGLR))
clusterExport(cl, list("ETA", "groups", "data"))

rxn <- parLapply(cl, 1:4, function(i) {
  BGLR(y = data$BLUE, ETA = ETA, burnIn = 1e5, nIter = 1.1e6, thin = 1e3, 
       saveAt = paste0("data/bglr/rxn_chain", i, "_"), groups = groups)
})

stopCluster(cl)

write_rds(rxn, "data/phenotype/joint_rxn_norm_model.rds")


# Check chain convergence -------------------------------------------------
site_levels <- sort(unique(data$Site))
files <- list.files("data/bglr", "*", full.names = TRUE) %>%
  split(., rep(1:4, each = 15))
mcmc_res_list <- lapply(files, function(f) {
  temp <- lapply(f, function(f2) {
    if (str_detect(f2, "dat$") & str_detect(f2, "fixed")) {
      temp <- read_delim(f2, delim = " ")
    } else if (str_detect(f2, "dat$") & !str_detect(f2, "fixed")) {
      temp <- read_delim(f2, delim = " ", col_names = FALSE)
    } else {
      temp <- readBinMat(f2)
    }
    
    if (!str_detect(f2, "fixed")) {
      if (str_detect(f2, "varB")) {
        colnames(temp) <- str_replace(f2, "data/bglr/rxn_chain[1-4]_ETA_", "") %>%
          str_replace(., "\\.dat", "")
      } else if (str_detect(f2, "mu")) {
        colnames(temp) <- "mu"
      } else if (str_detect(f2, "varE")) {
        colnames(temp) <- site_levels[as.integer(names(rxn[[1]]$varE))]
      } else {
        comp <- str_replace(f2, "data/bglr/rxn_chain[1-4]_ETA_", "") %>%
          str_replace(., "_b\\.bin", "")
        colnames(temp) <- rxn[[1]]$ETA[[comp]]$colNames %>%
          str_replace(., "factor\\(PedigreeNew\\)", "")
      }
    }
    
    if (nrow(temp) > 1000) temp <- temp[101:1100, ]
    
    return(as.matrix(temp))
  })
  
  mcmc(do.call("cbind", temp), start = 1, end = 1000)
}) %>%
  as.mcmc.list()

gdiag <- geweke.diag(mcmc_res_list)
idx <- lapply(gdiag, function(g) which(abs(g$z) >= 1.96))
length(reduce(idx, intersect)) 
# Each parameter is properly mixed in >= 3/4 chains


# Collect the Gibbs samples -----------------------------------------------
gibbs <- lapply(mcmc_res_list, as.matrix)
gibbs <- do.call("rbind", gibbs)
params <- apply(gibbs, 2, mean) %>%
  enframe()

mus <- filter(params, name %in% c("mu", "TMIN_X0.225_X0.45", "TMIN_X0.45_X1.1", 
                                  "TMAX_X0.75_X0.95", "PPT_X0.1_X0.225", 
                                  "NET_X0.425_X1.125"))
ped <- filter(params, str_detect(name, "/")) %>%
  separate("name", c("PedigreeNew", "Variable"), sep = ":", remove = TRUE) %>%
  mutate(Variable = if_else(is.na(Variable), "mu", Variable)) %>%
  full_join(., mus, by = c("Variable" = "name")) %>%
  mutate(Value = value.x + value.y) %>%
  select(-value.x, -value.y)

ggplot(ped, aes(x = Value)) + theme_classic() +
  geom_density() + facet_wrap(~ Variable, scales = "free") + 
  labs(x = "", y = "Density")

write_rds(ped, "data/phenotype/rxn_norm_parameters.rds")
