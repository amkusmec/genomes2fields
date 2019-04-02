require(tidyverse)

manhattan_gemma <- function(assoc, cutoff = NULL, title = NULL,
                            ylab = expression(-log[10](p-value))) {
  # Re-format the association table
  assoc <- assoc %>%
    mutate(ps = 1:n(),
           chr_color = factor(chr %% 2),
           pval = -log10(pval))
  
  # Define chromosome boundaries
  bounds <- assoc %>%
    group_by(chr) %>%
    summarise(Bound = max(ps)) %>%
    ungroup() %>%
    bind_rows(tibble(chr = 0, Bound = 0), .)
  
  # Define x-axis ticks
  ticks <- assoc %>%
    group_by(chr) %>%
    summarise(Tick = (max(ps) - min(ps))/2 + min(ps)) %>%
    ungroup()
  
  # Get y-axis limits
  ymax <- ceiling(max(assoc$pval, na.rm = TRUE))
  
  # Make the basic manhattan plot
  p <- ggplot(assoc, aes(x = ps, y = pval, colour = chr_color)) + theme_bw() +
    geom_vline(data = bounds, mapping = aes(xintercept = Bound), colour = "grey80") +
    geom_point(size = 0.5) +
    scale_colour_manual(values = c("grey60", "black")) +
    scale_x_continuous(breaks = ticks$Tick, labels = ticks$chr) +
    # scale_y_continuous(breaks = 0:ymax, labels = 0:ymax) +
    labs(colour = "", x = "Chromosome", y = ylab) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "none")
  
  if (!is.null(cutoff)) {
    p <- p + geom_hline(yintercept = cutoff, linetype = 2, colour = "red")
  }
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}
