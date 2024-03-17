library(readxl)
library(vegan)
library(car)
library(ggplot2)
library(gridExtra)

# Load data
data <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics_project/spontaneous_combined.xlsx", col_names = FALSE)

data_abx_follow_up <- read.csv("C:/Users/shaya/OneDrive/Desktop/Antibiotics_project/total_abx_follow_up.csv", header = FALSE)

data <- apply(data, 2, as.numeric)
data_abx_follow_up <- apply(data_abx_follow_up, 2, as.numeric)

baseline <- data[1:130,]

spo_abx_1 = data_abx_follow_up[1:5, ]
spo_abx_2 = data_abx_follow_up[19:25, ]
spo_abx_3 = data_abx_follow_up[42:48, ]
spo_abx_4 = data_abx_follow_up[64:70, ]
spo_abx_5 = data_abx_follow_up[83:89, ]
spo_abx_6 = data_abx_follow_up[105:110, ]
spo_abx_7 = data_abx_follow_up[122:128, ]
  
spo_follow_up_1 = data_abx_follow_up[6:18, ]
spo_follow_up_2 = data_abx_follow_up[26:41, ]
spo_follow_up_3 = data_abx_follow_up[49:63, ]
spo_follow_up_4 = data_abx_follow_up[71:82, ]
spo_follow_up_5 = data_abx_follow_up[90:104, ]
spo_follow_up_6 = data_abx_follow_up[111:121, ]
spo_follow_up_7 = data_abx_follow_up[129:139, ]

PCoA <- function(data, method = "bray") {
 
  dissimilarity_matrix <- vegdist(data, method = method)
  coordinates <- cmdscale(dissimilarity_matrix, k = 2)
  
  return(coordinates)
}

data_spo_1 <- rbind(baseline, spo_abx_1, spo_follow_up_1)
data_spo_2 <- rbind(baseline, spo_abx_2, spo_follow_up_2)
data_spo_3 <- rbind(baseline, spo_abx_3, spo_follow_up_3)
data_spo_4 <- rbind(baseline, spo_abx_4, spo_follow_up_4)
data_spo_5 <- rbind(baseline, spo_abx_5, spo_follow_up_5)
data_spo_6 <- rbind(baseline, spo_abx_6, spo_follow_up_6)
data_spo_7 <- rbind(baseline, spo_abx_7, spo_follow_up_7)

coordinates_spo_1 = PCoA(data_spo_1)
coordinates_spo_2 = PCoA(data_spo_2)
coordinates_spo_3 = PCoA(data_spo_3)
coordinates_spo_4 = PCoA(data_spo_4)
coordinates_spo_5 = PCoA(data_spo_5)
coordinates_spo_6 = PCoA(data_spo_6)
coordinates_spo_7 = PCoA(data_spo_7)

colors_spo_1 <- rep("grey", dim(coordinates_spo_1)[1])

colors_spo_1[87:91] <- "green"
colors_spo_1[131:135] <- "red"        
colors_spo_1[136:148] <- "blue"

colors_spo_2 <- rep("grey", dim(coordinates_spo_2)[1])

colors_spo_2[92:98] <- "green"
colors_spo_2[131:137] <- "red"        
colors_spo_2[138:148] <- "blue"
colors_spo_2[149:153] <- "blue"

colors_spo_3 <- rep("grey", dim(coordinates_spo_3)[1])

colors_spo_3[99:105] <- "green"
colors_spo_3[131:137] <- "red"        
colors_spo_3[138:147] <- "blue"
colors_spo_3[148:152] <- "blue"

colors_spo_4 <- rep("grey", dim(coordinates_spo_4)[1])

colors_spo_4[106:111] <- "green"
colors_spo_4[131:137] <- "red"        
colors_spo_4[138:144] <- "blue"
colors_spo_4[145:149] <- "blue"

colors_spo_5 <- rep("grey", dim(coordinates_spo_5)[1])

colors_spo_5[112:118] <- "green"
colors_spo_5[131:137] <- "red"        
colors_spo_5[138:147] <- "blue"
colors_spo_5[148:152] <- "blue"

colors_spo_6 <- rep("grey", dim(coordinates_spo_6)[1])

colors_spo_6[119:123] <- "green"
colors_spo_6[131:136] <- "red"        
colors_spo_6[137:142] <- "blue"
colors_spo_6[143:147] <- "blue"

colors_spo_7 <- rep("grey", dim(coordinates_spo_7)[1])

colors_spo_7[124:130] <- "green"
colors_spo_7[131:137] <- "red"        
colors_spo_7[138:143] <- "blue"
colors_spo_7[144:148] <- "blue"

create_custom_plot <- function(coordinates, colors, m, title = NULL) {
  
  data_df <- data.frame(PCoA1 = coordinates[, 1],
                        PCoA2 = coordinates[, 2], CustomColor = colors)
  
  p <- ggplot(data_df, aes(x = PCoA1, y = PCoA2, color = CustomColor)) +
    geom_point(size = 3, alpha = 0.6, shape = 16) +  # Set alpha and shape for transparency and filled points
    geom_point(data = subset(data_df, CustomColor == "green"),
               aes(color = "green"), size = 3, , alpha = 0.6, shape = 16, fill = NA) +
    xlab("PCo1") +
    ylab("PCo2") +
    scale_color_identity() +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = "black"),
          panel.grid = element_blank())
  
  p <- p + stat_ellipse(data = subset(data_df, CustomColor == "green"),
                        aes(group = CustomColor), type = "norm", level = 0.95)
  
  if (m > 0 && m < nrow(coordinates)) {
    # Extract the last m rows of coordinates
    last_m_points <- coordinates[(nrow(coordinates) - m + 1):nrow(coordinates), ]
    
    # Create a data frame to store line segments
    line_df <- data.frame()
    
    # Create line segments connecting the last m points in chronological order
    for (i in 1:(m - 1)) {
      segment <- data.frame(PCoA1 = last_m_points[i, 1],
                            PCoA2 = last_m_points[i, 2],
                            PCoA1_end = last_m_points[i + 1, 1],
                            PCoA2_end = last_m_points[i + 1, 2])
      line_df <- rbind(line_df, segment)
    }
    
    # Add line segments to the plot with a thinner line size (e.g., size = 0.5)
    p <- p + geom_segment(data = line_df, aes(xend = PCoA1_end, yend = PCoA2_end), color = "black", size = 0.2)
  }
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

#p1 = create_custom_plot(coordinates_spo_1, colors_spo_1, 18, 'Spontaneous 1')
p2 = create_custom_plot(coordinates_spo_2, colors_spo_2, 23, 'Subject #2')
p3 = create_custom_plot(coordinates_spo_3, colors_spo_3, 22, 'Subject #3')
p4 = create_custom_plot(coordinates_spo_4, colors_spo_4, 19, 'Subject #4')
p5 = create_custom_plot(coordinates_spo_5, colors_spo_5, 22, 'Subject #5')
#p6 = create_custom_plot(coordinates_spo_6, colors_spo_6, 17, 'Spontaneous 6')
#p7 = create_custom_plot(coordinates_spo_7, colors_spo_7, 18, 'Spontaneous 7')

#grid.arrange(p1, p2, p3, p4, p5, p6, p7, ncol = 4, nrow = 2)

grid.arrange(p2, p3, p4, p5, ncol = 2, nrow = 2)




























