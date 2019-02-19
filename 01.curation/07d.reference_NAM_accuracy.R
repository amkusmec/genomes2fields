library(tidyverse)


# Chromosome 1 ------------------------------------------------------------
chrom1 <- list.files("data/gbs", "NAM_merged_chrom1_acc_[0-9]{1,2}\\.csv", 
                     full.names = TRUE) %>%
  map(read_csv)
acc <- lapply(chrom1, function(x) x$acc) %>%
  bind_cols() %>% as.matrix() %>%
  rowMeans()
mean(acc) # ~85% imputation accuracy
summary(lm(acc ~ chrom1[[1]]$maf))

tibble(MAF = chrom1[[1]]$maf, 
       Accuracy = acc) %>%
  ggplot(., aes(x = MAF, y = Accuracy)) + theme_bw() +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
    ggtitle("Chromosome 1") + scale_fill_viridis_c() +
    labs(fill = "Density")
ggsave("figures/munge/NAM_reference_chrom1_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
rm(chrom1)


# Chromosome 2 ------------------------------------------------------------
chrom2 <- list.files("data/gbs", "NAM_merged_chrom2_acc_[0-9]{1,2}\\.csv", 
                     full.names = TRUE) %>%
  map(read_csv)
acc <- lapply(chrom2, function(x) x$acc) %>%
  bind_cols() %>% as.matrix() %>%
  rowMeans()
mean(acc) # ~85% imputation accuracy
summary(lm(acc ~ chrom2[[1]]$maf))

tibble(MAF = chrom2[[1]]$maf, 
       Accuracy = acc) %>%
  ggplot(., aes(x = MAF, y = Accuracy)) + theme_bw() +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
    ggtitle("Chromosome 2") + scale_fill_viridis_c() +
    labs(fill = "Density")
ggsave("figures/munge/NAM_reference_chrom2_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
rm(chrom2)


# Chromosome 3 ------------------------------------------------------------
chrom3 <- list.files("data/gbs", "NAM_merged_chrom3_acc_[0-9]{1,2}\\.csv", 
                     full.names = TRUE) %>%
  map(read_csv)
acc <- lapply(chrom3, function(x) x$acc) %>%
  bind_cols() %>% as.matrix() %>%
  rowMeans()
mean(acc) # ~85% imputation accuracy
summary(lm(acc ~ chrom3[[1]]$maf))

tibble(MAF = chrom3[[1]]$maf, 
       Accuracy = acc) %>%
  ggplot(., aes(x = MAF, y = Accuracy)) + theme_bw() +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
    ggtitle("Chromosome 3") + scale_fill_viridis_c() +
    labs(fill = "Density")
ggsave("figures/munge/NAM_reference_chrom3_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
rm(chrom3)


# Chromosome 4 ------------------------------------------------------------
chrom4 <- list.files("data/gbs", "NAM_merged_chrom4_acc_[0-9]{1,2}\\.csv", 
                     full.names = TRUE) %>%
  map(read_csv)
acc <- lapply(chrom4, function(x) x$acc) %>%
  bind_cols() %>% as.matrix() %>%
  rowMeans()
mean(acc) # ~85% imputation accuracy
summary(lm(acc ~ chrom4[[1]]$maf))

tibble(MAF = chrom4[[1]]$maf, 
       Accuracy = acc) %>%
  ggplot(., aes(x = MAF, y = Accuracy)) + theme_bw() +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
    ggtitle("Chromosome 4") + scale_fill_viridis_c() +
    labs(fill = "Density")
ggsave("figures/munge/NAM_reference_chrom4_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
rm(chrom4)


# Chromosome 5 ------------------------------------------------------------
chrom5 <- list.files("data/gbs", "NAM_merged_chrom5_acc_[0-9]{1,2}\\.csv", 
                     full.names = TRUE) %>%
  map(read_csv)
acc <- lapply(chrom5, function(x) x$acc) %>%
  bind_cols() %>% as.matrix() %>%
  rowMeans()
mean(acc) # ~85% imputation accuracy
summary(lm(acc ~ chrom5[[1]]$maf))

tibble(MAF = chrom5[[1]]$maf, 
       Accuracy = acc) %>%
  ggplot(., aes(x = MAF, y = Accuracy)) + theme_bw() +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
    ggtitle("Chromosome 5") + scale_fill_viridis_c() +
    labs(fill = "Density")
ggsave("figures/munge/NAM_reference_chrom5_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
rm(chrom5)


# Chromosome 6 ------------------------------------------------------------
chrom6 <- list.files("data/gbs", "NAM_merged_chrom6_acc_[0-9]{1,2}\\.csv", 
                     full.names = TRUE) %>%
  map(read_csv)
acc <- lapply(chrom6, function(x) x$acc) %>%
  bind_cols() %>% as.matrix() %>%
  rowMeans()
mean(acc) # ~85% imputation accuracy
summary(lm(acc ~ chrom6[[1]]$maf))

tibble(MAF = chrom6[[1]]$maf, 
       Accuracy = acc) %>%
  ggplot(., aes(x = MAF, y = Accuracy)) + theme_bw() +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
    ggtitle("Chromosome 6") + scale_fill_viridis_c() +
    labs(fill = "Density")
ggsave("figures/munge/NAM_reference_chrom6_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
rm(chrom6)


# Chromosome 7 ------------------------------------------------------------
chrom7 <- list.files("data/gbs", "NAM_merged_chrom7_acc_[0-9]{1,2}\\.csv", 
                     full.names = TRUE) %>%
  map(read_csv)
acc <- lapply(chrom7, function(x) x$acc) %>%
  bind_cols() %>% as.matrix() %>%
  rowMeans()
mean(acc) # ~85% imputation accuracy
summary(lm(acc ~ chrom7[[1]]$maf))

tibble(MAF = chrom7[[1]]$maf, 
       Accuracy = acc) %>%
  ggplot(., aes(x = MAF, y = Accuracy)) + theme_bw() +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
    ggtitle("Chromosome 7") + scale_fill_viridis_c() +
    labs(fill = "Density")
ggsave("figures/munge/NAM_reference_chrom7_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
rm(chrom7)


# Chromosome 8 ------------------------------------------------------------
chrom8 <- list.files("data/gbs", "NAM_merged_chrom8_acc_[0-9]{1,2}\\.csv", 
                     full.names = TRUE) %>%
  map(read_csv)
acc <- lapply(chrom8, function(x) x$acc) %>%
  bind_cols() %>% as.matrix() %>%
  rowMeans()
mean(acc) # ~85% imputation accuracy
summary(lm(acc ~ chrom8[[1]]$maf))

tibble(MAF = chrom8[[1]]$maf, 
       Accuracy = acc) %>%
  ggplot(., aes(x = MAF, y = Accuracy)) + theme_bw() +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
    ggtitle("Chromosome 8") + scale_fill_viridis_c() +
    labs(fill = "Density")
ggsave("figures/munge/NAM_reference_chrom8_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
rm(chrom8)


# Chromosome 9 ------------------------------------------------------------
chrom9 <- list.files("data/gbs", "NAM_merged_chrom9_acc_[0-9]{1,2}\\.csv", 
                     full.names = TRUE) %>%
  map(read_csv)
acc <- lapply(chrom9, function(x) x$acc) %>%
  bind_cols() %>% as.matrix() %>%
  rowMeans()
mean(acc) # ~85% imputation accuracy
summary(lm(acc ~ chrom9[[1]]$maf))

tibble(MAF = chrom9[[1]]$maf, 
       Accuracy = acc) %>%
  ggplot(., aes(x = MAF, y = Accuracy)) + theme_bw() +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
    ggtitle("Chromosome 9") + scale_fill_viridis_c() +
    labs(fill = "Density")
ggsave("figures/munge/NAM_reference_chrom9_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
rm(chrom9)


# Chromosome 10 -----------------------------------------------------------
chrom10 <- list.files("data/gbs", "NAM_merged_chrom10_acc_[0-9]{1,2}\\.csv", 
                     full.names = TRUE) %>%
  map(read_csv)
acc <- lapply(chrom10, function(x) x$acc) %>%
  bind_cols() %>% as.matrix() %>%
  rowMeans()
mean(acc) # ~85% imputation accuracy
summary(lm(acc ~ chrom10[[1]]$maf))

tibble(MAF = chrom10[[1]]$maf, 
       Accuracy = acc) %>%
  ggplot(., aes(x = MAF, y = Accuracy)) + theme_bw() +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") + 
    ggtitle("Chromosome 10") + scale_fill_viridis_c() +
    labs(fill = "Density")
ggsave("figures/munge/NAM_reference_chrom10_accuracy.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)
rm(chrom10)
