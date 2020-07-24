### Implements the deltaK method.

library(tidyverse)


# Parse fastStructure output ----------------------------------------------
ll <- list.files("data/structure", "g2f_hyb_rep[0-9]{1,2}\\.[0-9]{1,2}\\.log", 
                 full.names = TRUE) %>%
  map_dbl(function(f) {
    temp <- read_lines(f)
    str_split(temp[length(temp) - 2], " ") %>%
      sapply(., function(x) x[4]) %>%
      as.double()
  })
ll <- tibble(K = list.files("data/structure", "g2f_hyb_rep[0-9]{1,2}\\.[0-9]{1,2}\\.log") %>%
               str_remove(., "g2f_hyb_rep[0-9]{1,2}\\.") %>%
               str_remove(., "\\.log") %>% as.integer(), 
             Rep = list.files("data/structure", "g2f_hyb_rep[0-9]{1,2}\\.[0-9]{1,2}\\.log") %>%
               str_remove(., "g2f_hyb_rep") %>%
               str_remove(., "\\.[0-9]{1,2}\\.log") %>% as.integer(), 
             LL = ll) %>%
  arrange(K, Rep)


# Summarize changes in the log-likelihood with respect to K ---------------
lk <- ll %>%
  group_by(K) %>%
  summarise(LK = mean(LL), 
            SD = sd(LL)) %>%
  ungroup()

ggplot(lk, aes(x = K, y = LK)) + theme_bw() + 
  geom_pointrange(aes(ymin = LK - 2*SD, ymax = LK + 2*SD), size = 0.5) +
  labs(y = expression(L(K)))
ggsave("figures/munge/fastSTRUCTURE_ll.pdf", width = 6, height = 4, units = "in", dpi = 300)


# Calculate the delta-K statistic -----------------------------------------
# Evanno et al. (2011)
lpk <- ll %>%
  group_by(K) %>%
  mutate(LpK = c(0, LL[-1] - LL[-n()]),
         LppK = abs(c(LpK[-1] - LpK[-n()], 0))) %>%
  ungroup()
deltaK <- lpk %>%
  group_by(K) %>%
  summarise(m = mean(LppK), 
            s = sd(LL), 
            DK = m/s) %>%
  ungroup()
ggplot(deltaK, aes(x = K, y = DK)) + theme_classic() + geom_line() +
  scale_x_continuous(breaks = 2:15) +
  labs(y = expression(Delta*K))
ggsave("figures/munge/fastSTRUCTURE_deltaK.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)

deltaK$K[which.max(deltaK$DK)]
