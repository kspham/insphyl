library(ggplot2)
library(reshape2)
library(ggdendro)
library(cowplot)
library(jsonlite)


plot_boolean <- function(meta_json, path_cor, path_out, refname) {
    # path_cor <- "mrsa_process/result/cor.csv"

    data <- read.csv(path_cor)
    meta <- fromJSON(meta_json)
    # Collect status
    status <- as.character(unlist(lapply(meta, function(x) x$status)))
    status <- data.frame("sample" = names(meta), status, stringsAsFactors = F)
    status <- rbind(status, list("sample" = refname, "status" = "reference"))
    # Replace digit names
    digit_names <- !is.na(lapply(substr(status$sample, 1, 1), as.numeric))
    status$sample[digit_names] <- paste0("X", status$sample[digit_names])
    # Dendro
    hc <- hclust(dist(x = t(data[, -1])), "average")
    Order <- c(1, hc$order + 1)
    group <- dendro_data(hc)$labels
    group$status <- status$status[match(group$label, status$sample)]
    dd <- ggdendrogram(hc)
    dd <- dd + geom_point(data = group, aes(x = x, y = y, color = status))
    dd <- dd + guides(color = F) + theme(
        axis.ticks.y = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0, 0, -1, 1), "cm")
    )
    # Heatmap
    hcr <- hclust(dist(x = data[, -1]), "average")
    hm_data <- melt(data[hcr$order, Order], id.vars = "X")
    hm_data$variable <- factor(hm_data$variable, levels = colnames(data)[Order])
    hm <- ggplot(hm_data, aes(x = variable, y = X, fill = as.character(value)))
    hm <- hm + geom_raster() + guides(fill = F)
    hm <- hm + theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
    hm <- hm + scale_fill_brewer(palette = "Paired")
    p <- plot_grid(dd, hm, nrow = 2, align = "v", rel_heights = c(1, 1.2))
    ggsave(p, device = "pdf", filename = path_out)
    return ("Done")
}


if (!interactive()) {
    args <- commandArgs(trailing = T)
    do.call(args[1], as.list(args[-1]))
}