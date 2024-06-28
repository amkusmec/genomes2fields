library(tidyverse)
library(parallel)
library(grid)
library(gridExtra)


# Prepare the data --------------------------------------------------------
data <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$")) %>%
  filter(Site != "NEH3_2015") %>%
  separate(Site, c("Environment", "Year"), sep = "_", remove = FALSE) %>%
  select(BLUE, PedigreeNew, Site, Environment, Year, everything())
net <- names(data)[which(str_detect(names(data), "NET"))]
data <- mutate_at(data, net, function(x) -1*x)
ped_site <- data %>%
  select(Site, PedigreeNew) %>%
  split(., .$Site)

# Response
y <- data$BLUE

# Predictors
X <- cbind(model.matrix(~ 1 + Site, data = data), 
           as.matrix(data[, -c(1:5)]))
X[, -c(1:45)] <- scale(X[, -c(1:45)], center = TRUE, scale = TRUE)


# Variance-covariance matrix of the BLUEs
d <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
d <- d[names(d) %in% names(ped_site)]
d <- mapply(function(x, y) {
  idx <- which(rownames(x$vcov) %in% y$PedigreeNew)
  temp <- x$vcov[idx, idx]
  c(temp[1, 1], temp[1, 1] + diag(temp)[-1] + 2*temp[-1, 1])
}, x = d, y = ped_site)
wts <- 1/unlist(d, use.names = FALSE)



g <- read_rds("data/weather/ga_het_resid_mean.rds")

# Reduce the number of variables ------------------------------------------
# Get all the variables that were selected in any model
nvars <- sapply(g, function(x) length(x$g))
vars <- lapply(g, function(x) x$g) %>% unlist(use.names = FALSE)
minima <- sapply(g, function(x) min(x$minima))
niters <- sapply(g, function(x) length(x$minima))


tt <- table(vars)
var_table <- tibble(Variable = names(tt)) %>%
  separate(Variable, c("EVar", "Start", "End"), sep = "_", remove = FALSE) %>%
  mutate(Start = str_remove(Start, "X") %>% as.numeric(),
         End = str_remove(End, "X") %>% as.numeric()) %>%
  group_by(EVar) %>%
  arrange(Start, End) %>%
  mutate(Index = 1:n()) %>%
  ungroup()

X2 <- X[, c(1:45, which(colnames(X) %in% var_table$Variable))]
g2 <- read_rds("data/weather/ga_final_model.rds")

# Behavior of GA over replicates and final model
minima <- g %>%
  map(function(x) {
    c(x$minima, rep(NA, max(niters) - length(x$minima)))
  }) %>%
  unlist() %>%
  matrix(., nrow = max(niters), ncol = 100, byrow = FALSE)

perf <- as_tibble(minima) %>%
  gather(Replicate, Minimum, everything()) %>%
  filter(!is.na(Minimum)) %>%
  group_by(Replicate) %>%
  mutate(Iteration = 1:n()) %>%
  ungroup()


# Evaluate all possible models --------------------------------------------
models <- lapply(1:5, function(i) {
  temp <- combn(g2$g, i)
  if (nrow(temp) < 5) {
    temp <- rbind(temp, matrix(NA, ncol = ncol(temp), nrow = 5 - i))
  }
  return(temp)
})
models <- do.call("cbind", models)

bic <- apply(models, 2, function(v) {
  v <- v[!is.na(v)]
  temp <- X[, c(1:45, which(colnames(X) %in% v))]
  BIC(lm(y ~ 0 + temp, weights = wts)) + 2*1*log(choose(ncol(X) - 45, length(v)))
})

# Add the model without environmental variables
# N.B.: The ith column of `models` corresponds to the (i+1)th element of `bic`.
bic <- c(BIC(lm(y ~ 0 + X[, 1:45], weights = wts)) + 
           2*1*log(choose(ncol(X) - 45, 0)), bic)


selected <- models[, which.min(bic) - 1]

pA <- ggplot(perf, aes(x = Iteration, y = Minimum)) + theme_classic() +
  geom_line(aes(group = Replicate), colour = "grey80", alpha = 0.5) +
  geom_line(aes(y = M), data = perf %>%
              group_by(Iteration) %>%
              summarise(M = mean(Minimum))) +
  geom_line(data = tibble(Iteration = 1:length(g2$minima),
                          Minimum = g2$minima),
            colour = "red", linewidth = 1) +
  geom_hline(yintercept = min(bic), linetype = 2) + 
  labs(y = "eBIC", tags = "(a)") + 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10))

evar_lbls <- function(y) {
  sapply(y, function(x) {
    if (x == "NET") {
      expression(bold("NET (mm)"))
    } else if (x == "PPT") {
      expression(bold("PPT (mm)")) 
    } else if (x == "SR") {
      expression(bold("SR (W "*m^-2*")"))
    } else if (x == "TMAX") {
      expression(bold("TMAX ("*degree*"C)"))
    } else if (x == "TMIN") {
      expression(bold("TMIN ("*degree*"C)"))
    }
  })
}

pB <- ggplot(var_table, aes(y = Index)) + theme_classic() +
  geom_segment(aes(yend = Index, x = Start, xend = End), linewidth = 0.25) +
  facet_wrap(~ EVar, ncol = 2, 
             labeller = as_labeller(evar_lbls, 
                                    default = label_parsed)) + 
  labs(x = "% CHU to Anthesis", y = "", tag = "(b)") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.5)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size = 8), 
        axis.title = element_text(size = 10), 
        strip.text = element_text(size = 8))

# Plot of selected variables
pC <- tibble(V = selected[!is.na(selected)]) %>%
  separate(V, c("Variable", "Start", "End"), sep = "_", remove = TRUE) %>%
  mutate_at(c("Start", "End"), str_remove, pattern = "X") %>%
  mutate_at(c("Start", "End"), as.numeric) %>%
  mutate(Ypos = factor(Variable) %>% as.integer(),
         Ypos = Ypos + c(0.125, -0.125, 0, 0)) %>%
  ggplot(.) + theme_classic() +
  geom_segment(aes(x = Start, xend = End, y = Ypos, yend = Ypos,
                   colour = Variable), size = 2) +
  geom_vline(xintercept = 1, linetype = 2) +
  labs(x = "% CHU to anthesis", y = "", tag = "(c)") + guides(colour = "none") +
  scale_colour_manual(values = c("TMAX" = "red", "SR" = "goldenrod", "NET" = "brown")) +
  scale_y_continuous(breaks = 1:3, labels = c("NET", "SR" ,"TMAX")) +
  scale_x_continuous(limits = c(0, 1.5), labels = scales::percent) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 12))



lay <- matrix(c(1, 2, 3, 2), ncol = 2, byrow = TRUE)
gp <- arrangeGrob(ggplotGrob(pA), ggplotGrob(pB), ggplotGrob(pC), 
                  layout_matrix = lay)
ggsave("figures/final/Figure_2.pdf", gp, width = 180, height = 90, units = "mm", dpi = 600)
