library(readxl)
library(vegan)
library(car)
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Load data
data <- read_excel("C:/Users/USER/Desktop/Antibiotics/Eran Elinav/Data/spontaneous_combined.xlsx", col_names = FALSE)

data <- apply(data, 2, as.numeric)

baseline <- data[1:130,]

PCoA <- function(data, method = "bray") {
  
  dissimilarity_matrix <- vegdist(data, method = method)
  coordinates <- cmdscale(dissimilarity_matrix, k = 2)
  
  return(coordinates)
}

coordinates = PCoA(baseline)

colors <- rep("grey", dim(coordinates)[1])

colors[119:123] <- "green"
colors[124:130] <- "red"
colors[99:105] <- "brown"
colors[1:7] <- "black"
colors[20:25] <- "blue"

create_custom_plot <- function(coordinates, colors, m, title = NULL) {
  data_df <- data.frame(PCoA1 = coordinates[, 1],
                        PCoA2 = coordinates[, 2], 
                        CustomColor = colors)
  
  p <- ggplot(data_df, aes(x = PCoA1, y = PCoA2)) +
    geom_point(aes(color = CustomColor), size = 3, alpha = 0.6, shape = 16) +
    scale_color_manual(values = c("grey" = "grey", "green" = "black"
                                  , "red" = "black", "blue" = "black",
                                  "brown" = "black", "black" = "black")) +
    xlab(TeX("$PCo1$")) +
    ylab(TeX("$PCo2$")) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = "black"),
          panel.grid = element_blank(),
          axis.title.x = element_text(size = 20),  
          axis.title.y = element_text(size = 20),  
          axis.text.x = element_text(size = 14),   
          axis.text.y = element_text(size = 14)) +
    guides(color = FALSE)  
  
  groups <- unique(colors)
  for (group in groups) {
    if (group != "grey") {  
      p <- p + stat_ellipse(data = subset(data_df, CustomColor == group),
                            aes(x = PCoA1, y = PCoA2, group = CustomColor), 
                            color = "black", type = "norm", level = 0.95)
    }
  }
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

p = create_custom_plot(coordinates, colors, 18)

print(p)