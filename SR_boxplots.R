library(ggplot2)
library(readxl)

## Eran Elinav

# Load data
baseline_data <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/Eran Elinav/Data/num_species_base.xlsx", col_names = TRUE)
abx_data <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/Eran Elinav/Data/num_species_ABX.xlsx", col_names = TRUE)
post_abx_data <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/Eran Elinav/Data/num_species_follow.xlsx", col_names = TRUE)

# Transform to vectors
baseline_data <- baseline_data$`Species richness`
abx_data <- abx_data$`Species richness`
post_abx_data <- post_abx_data$`Species richness`

# Create a data frame with the data
df <- data.frame(
  TimePoint = factor(rep(c('Baseline', 'ABX', 'Post ABX'), each = 7)),
  Value = c(baseline_data, abx_data, post_abx_data)
)

# Specify the order of factor levels for 'TimePoint'
df$TimePoint <- factor(df$TimePoint, levels = c('Baseline', 'ABX', 'Post ABX'))

# Create a new variable 'PointOrder' representing the order of points within each category
df$PointOrder <- 1:21  

custom_colors <- c(
  'red', 'blue', 'green', 'purple', 'brown', 'pink', 'yellow',  # Baseline
  'red', 'blue', 'green', 'purple', 'brown', 'pink', 'yellow',  # ABX
  'red', 'blue', 'green', 'purple', 'brown', 'pink', 'yellow'   # Post ABX
)

# Create a scatter plot using ggplot2 with custom colors
scatter_plot <- ggplot(df, aes(x = TimePoint, y = Value, color = as.factor(PointOrder))) +
  geom_point(size = 3) +  # Add scatter points
  labs(
    x = '',
    y = 'Species Richness'
  ) +
  scale_color_manual(values = custom_colors) +  
  theme_minimal() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 14),  # Adjust axis title size
        axis.text = element_text(size = 12),   # Adjust axis text size
        axis.text.y = element_text(size = 12)) +  # Adjust y-axis text size
  ylab("Species Richness")  # Adjust the y-axis label

# Display the scatter plot
print(scatter_plot)

## Recovery

# Load data
baseline_data_recovery <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/Recovery/Data/num_species_base_recovery.xlsx", col_names = TRUE)
abx_data_recovery <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/Recovery/Data/num_species_ABX_recovery.xlsx", col_names = TRUE)
post_abx_data_recovery <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/Recovery/Data/num_species_follow_recovery.xlsx", col_names = TRUE)

# Transform to vectors
baseline_data_recovery <- baseline_data_recovery$`Species richness`
abx_data_recovery <- abx_data_recovery$`Species richness`
post_abx_data_recovery <- post_abx_data_recovery$`Species richness`

# Create a data frame with the data
df_recovery <- data.frame(
  TimePoint = factor(rep(c('Baseline', 'ABX', 'Post ABX'), each = 9)),
  Value = c(baseline_data_recovery, abx_data_recovery, post_abx_data_recovery)
)

# Specify the order of factor levels for 'TimePoint'
df_recovery$TimePoint <- factor(df_recovery$TimePoint, levels = c('Baseline', 'ABX', 'Post ABX'))

# Create a new variable 'PointOrder' representing the order of points within each category
df_recovery$PointOrder <- 1:27  

custom_colors <- c(
  'red', 'blue', 'green', 'purple', 'brown', 'pink', 'yellow','black','cyan',  # Baseline
  'red', 'blue', 'green', 'purple', 'brown', 'pink', 'yellow','black','cyan',  # ABX
  'red', 'blue', 'green', 'purple', 'brown', 'pink', 'yellow','black','cyan'   # Post ABX
)

# Create a scatter plot using ggplot2 with custom colors
scatter_plot <- ggplot(df_recovery, aes(x = TimePoint, y = Value, color = as.factor(PointOrder))) +
  geom_point(size = 3) +  # Add scatter points
  labs(
    x = '',
    y = 'Species Richness'
  ) +
  scale_color_manual(values = custom_colors) +  
  theme_minimal() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 14),  # Adjust axis title size
        axis.text = element_text(size = 12),   # Adjust axis text size
        axis.text.y = element_text(size = 12)) +  # Adjust y-axis text size
  ylab("Species Richness")  # Adjust the y-axis label

# Display the scatter plot
print(scatter_plot)

## Gut Bacterial Microbiota

# Load data
baseline_data_gut <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/gut/Data/num_species_base_gut.xlsx", col_names = TRUE)
abx_data_gut <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/gut/Data/num_species_ABX_gut.xlsx", col_names = TRUE)
post_abx_data_gut <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/gut/Data/num_species_follow_gut.xlsx", col_names = TRUE)

# Transform to vectors
baseline_data_gut <- baseline_data_gut$`Species richness`
abx_data_gut <- abx_data_gut$`Species richness`
post_abx_data_gut <- post_abx_data_gut$`Species richness`

# Create a data frame with the data
df_gut <- data.frame(
  TimePoint = factor(rep(c('Baseline', 'ABX', 'Post ABX'), each = 35)),
  Value = c(baseline_data_gut, abx_data_gut, post_abx_data_gut)
)

# Specify the order of factor levels for 'TimePoint'
df_gut$TimePoint <- factor(df_gut$TimePoint, levels = c('Baseline', 'ABX', 'Post ABX'))

# Create a new variable 'PointOrder' representing the order of points within each category
df_gut$PointOrder <- 1:105  

custom_colors <- c(rep('grey', each=105)
)

# Create a scatter plot using ggplot2 with custom colors
scatter_plot <- ggplot(df_gut, aes(x = TimePoint, y = Value, color = as.factor(PointOrder))) +
  geom_point(size = 3) +  # Add scatter points
  labs(
    x = '',
    y = 'Species Richness'
  ) +
  scale_color_manual(values = custom_colors) +  
  theme_minimal() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 14),  # Adjust axis title size
        axis.text = element_text(size = 12),   # Adjust axis text size
        axis.text.y = element_text(size = 12)) +  # Adjust y-axis text size
  ylab("Species Richness")  # Adjust the y-axis label

# Display the scatter plot
print(scatter_plot)

## Effects

# Load data
baseline_data_effects <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/Effects of different amoxicillin/Data/num_species_base_effects.xlsx", col_names = TRUE)
abx_data_effects <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/Effects of different amoxicillin/Data/num_species_ABX_effects.xlsx", col_names = TRUE)
post_abx_data_effects <- read_excel("C:/Users/shaya/OneDrive/Desktop/Antibiotics project/Effects of different amoxicillin/Data/num_species_follow_effects.xlsx", col_names = TRUE)

# Transform to vectors
baseline_data_effects <- baseline_data_effects$`Species richness`
abx_data_effects <- abx_data_effects$`Species richness`
post_abx_data_effects <- post_abx_data_effects$`Species richness`

# Create a data frame with the data
df_effects <- data.frame(
  TimePoint = factor(rep(c('Baseline', 'ABX', 'Post ABX'), each = 15)),
  Value = c(baseline_data_effects, abx_data_effects, post_abx_data_effects)
)

# Specify the order of factor levels for 'TimePoint'
df_effects$TimePoint <- factor(df_effects$TimePoint, levels = c('Baseline', 'ABX', 'Post ABX'))

# Create a new variable 'PointOrder' representing the order of points within each category
df_effects$PointOrder <- 1:45  

custom_colors <- c(
  'red', 'blue', 'green', 'purple', 'brown', 'pink', 'yellow','black','cyan','grey', 'salmon' ,'coral','deepskyblue','orange', 'darkgreen',  # Baseline
  'red', 'blue', 'green', 'purple', 'brown', 'pink', 'yellow','black','cyan','grey', 'salmon' ,'coral','deepskyblue','orange', 'darkgreen',  # ABX
  'red', 'blue', 'green', 'purple', 'brown', 'pink', 'yellow','black','cyan','grey', 'salmon' ,'coral','deepskyblue','orange', 'darkgreen'   # Post ABX
)

# Create a scatter plot using ggplot2 with custom colors
scatter_plot <- ggplot(df_effects, aes(x = TimePoint, y = Value, color = as.factor(PointOrder))) +
  geom_point(size = 3) +  # Add scatter points
  labs(
    x = '',
    y = 'Species Richness'
  ) +
  scale_color_manual(values = custom_colors) +  
  theme_minimal() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 14),  # Adjust axis title size
        axis.text = element_text(size = 12),   # Adjust axis text size
        axis.text.y = element_text(size = 12)) +  # Adjust y-axis text size
  ylab("Species Richness")  # Adjust the y-axis label

# Display the scatter plot
print(scatter_plot)