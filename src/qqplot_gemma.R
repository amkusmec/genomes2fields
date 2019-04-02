require(tidyverse)

qqplot_gemma <- function(p_values, title = NULL) {
  p_values <- p_values[!is.na(p_values) & p_values > 0 & p_values <= 1]
  if (length(p_values) < 1) return(NULL)
  
  # Compute quantiles
  p_values <- p_values[order(p_values)]
  quantiles <- (1:length(p_values))/(length(p_values) + 1)
  log_p_values <- -log10(p_values)
  log_quantiles <- -log10(quantiles)
  
  res <- tibble(p_values, quantiles, log_p_values, log_quantiles)
  
  # Compute 90% confidence interval
  res <- res %>%
    mutate(j = ceiling((10^-log_quantiles)*length(log_quantiles)),
           j = if_else(j == 0, 1, j),
           ci95 = qbeta(0.95, j, n() - j + 1),
           ci05 = qbeta(0.05, j, n() - j + 1),
           ci95 = -log10(ci95),
           ci05 = -log10(ci05))
  
  p <- ggplot(res) + theme_bw() +
    scale_x_continuous(breaks = 0:ceiling(max(log_quantiles)),
                       labels = 0:ceiling(max(log_quantiles))) +
    scale_y_continuous(breaks = 0:ceiling(max(log_p_values)),
                       labels = 0:ceiling(max(log_p_values))) +
    geom_ribbon(aes(x = log_quantiles, ymin = ci05, ymax = ci95), fill = "grey80") +
    geom_point(aes(x = log_quantiles, y = log_p_values), size = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "red") +
    labs(x = expression(Expected~-log[10](p-value)),
         y = expression(Observed~-log[10](p-value))) +
    theme(panel.grid.minor = element_blank())
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}
