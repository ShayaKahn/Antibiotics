library(ggplot2)
library(readxl)
library(vegan)
library(gridExtra)

########## Functions ##########

filter_low_abundance <- function(cohort, threshold) {
  ind <- which(cohort > threshold)
  cohort[-ind] <- 0.0
  cohort <- cohort / rowSums(cohort)
  return(cohort)
}

pairwise_jaccard_similarity <- function(ABX_cohort, post_ABX_cohort) {
  # Number of rows in the matrix
  n <- nrow(ABX_cohort)
  
  # Initialize a vector to store the Jaccard similarities
  jaccard_sim_ABX <- numeric()
  jaccard_sim_post_ABX <- numeric()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Find indices of nonzero elements in each row in ABX_cohort
      nonzero_i_ABX <- which(ABX_cohort[i, ] != 0)
      nonzero_j_ABX <- which(ABX_cohort[j, ] != 0)
      
      # Calculate intersection and union for presence/absence in ABX_cohort
      inersection_ind_ABX <- intersect(nonzero_i_ABX, nonzero_j_ABX)
      intersection_ABX <- length(inersection_ind_ABX)
      union_ABX_ind <- union(nonzero_i_ABX, nonzero_j_ABX)
      union_ABX <- length(union_ABX_ind)
      
      # Calculate Jaccard similarity
      similarity_ABX <- intersection_ABX / union_ABX
      
      # Append the similarity to the vector
      jaccard_sim_ABX <- c(jaccard_sim_ABX, similarity_ABX)
      
      # Find indices of nonzero elements in each row in post_ABX_cohort
      nonzero_i_post_ABX <- which(post_ABX_cohort[i, ][-inersection_ind_ABX] != 0)
      nonzero_j_post_ABX <- which(post_ABX_cohort[j, ][-inersection_ind_ABX] != 0)
      #nonzero_i_post_ABX <- which(post_ABX_cohort[i, ][-union_ABX] != 0)
      #nonzero_j_post_ABX <- which(post_ABX_cohort[j, ][-union_ABX] != 0)
      
      # Calculate intersection and union for presence/absence in post_ABX_cohort
      intersection_post_ABX <- length(intersect(nonzero_i_post_ABX, nonzero_j_post_ABX))
      union_post_ABX <- length(union(nonzero_i_post_ABX, nonzero_j_post_ABX))
      
      # Calculate Jaccard similarity
      similarity_post_ABX <- intersection_post_ABX / union_post_ABX
      
      # Append the similarity to the vector
      jaccard_sim_post_ABX <- c(jaccard_sim_post_ABX, similarity_post_ABX)
    }
  }
  return(list("jaccard_ABX" = jaccard_sim_ABX,
              "jaccard_post_ABX" = jaccard_sim_post_ABX))
}

create_correlation_plot <-function(x,y,width,hight,x_title,y_title){
  # Function to generate the plots
  p = ggplot(data = data.frame(x, y),
             aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "blue") + 
    geom_jitter(width=width, height=hight) +
    theme_minimal() +
    labs(x = x_title,
         y = y_title) +
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),  
          axis.text.x = element_text(size = 14),  
          axis.text.y = element_text(size = 14),  
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16),
          axis.line = element_line(color = "black"))
  return(p)
}

########## Load data ##########

# Eran Elinav

spo_baseline <- read_excel("C:/Users/USER/Desktop/Antibiotics/Eran Elinav/Data/Spo_baseline_cohort.xlsx",
                           col_names = FALSE)
spo_baseline <- as.matrix(spo_baseline)
spo_baseline <- spo_baseline[, -1]
spo_baseline <- filter_low_abundance(spo_baseline, 0.5*1e-2)
spo_ABX <- read_excel("C:/Users/USER/Desktop/Antibiotics/Eran Elinav/Data/Spo_ABX_cohort.xlsx",
                      col_names = FALSE)
spo_ABX <- as.matrix(spo_ABX)
spo_ABX <- spo_ABX[, -1]
spo_ABX <- filter_low_abundance(spo_ABX, 0.5*1e-2)
spo_post_ABX <- read_excel("C:/Users/USER/Desktop/Antibiotics/Eran Elinav/Data/Spo_post_ABX_cohort.xlsx",
                           col_names = FALSE)
spo_post_ABX <- as.matrix(spo_post_ABX)
spo_post_ABX <- spo_post_ABX[, -1]
spo_post_ABX <- filter_low_abundance(spo_post_ABX, 0.5*1e-2)

# Recovery

recovery_baseline <- read_excel("C:/Users/USER/Desktop/Antibiotics/Recovery/Data/baseline_rel_abund_rarefied_appear_4.xlsx",
                                col_names = FALSE)
recovery_baseline <- as.matrix(recovery_baseline)
recovery_baseline <- filter_low_abundance(recovery_baseline, 1e-3)
recovery_ABX <- read_excel("C:/Users/USER/Desktop/Antibiotics/Recovery/Data/rel_abund_rarefied_4.xlsx",
                           col_names = FALSE)
recovery_ABX <- as.matrix(recovery_ABX)
recovery_ABX <- filter_low_abundance(recovery_ABX, 1e-3)
recovery_post_ABX <- read_excel("C:/Users/USER/Desktop/Antibiotics/Recovery/Data/rel_abund_rarefied_180_appear_4.xlsx",
                                col_names = FALSE)
recovery_post_ABX <- as.matrix(recovery_post_ABX)
recovery_post_ABX <- filter_low_abundance(recovery_post_ABX, 1e-3)

# Gut

gut_baseline <- read_excel("C:/Users/USER/Desktop/Antibiotics/Gut/Data/placebo_data_v2.xlsx",
                           col_names = FALSE)
gut_baseline <- as.matrix(gut_baseline)
gut_ABX <- read_excel("C:/Users/USER/Desktop/Antibiotics/Gut/Data/placebo_data_v3.xlsx",
                      col_names = FALSE)
gut_ABX <- as.matrix(gut_ABX)
gut_post_ABX <- read_excel("C:/Users/USER/Desktop/Antibiotics/Gut/Data/placebo_data_v5.xlsx",
                           col_names = FALSE)
gut_post_ABX <- as.matrix(gut_post_ABX)

# Transpose and normalize the Gut data

column_sums_gut_baseline <- colSums(gut_baseline)
normalized_gut_baseline <- sweep(gut_baseline, 2, column_sums_gut_baseline, "/")
gut_baseline = t(normalized_gut_baseline)
gut_baseline <- filter_low_abundance(gut_baseline, 0.5*1e-2)
column_sums_gut_ABX <- colSums(gut_ABX)
normalized_gut_ABX <- sweep(gut_ABX, 2, column_sums_gut_ABX, "/")
gut_ABX = t(normalized_gut_ABX)
gut_ABX <- filter_low_abundance(gut_ABX, 0.5*1e-2)
column_sums_gut_post_ABX <- colSums(gut_post_ABX)
normalized_gut_post_ABX <- sweep(gut_post_ABX, 2, column_sums_gut_post_ABX, "/")
gut_post_ABX = t(normalized_gut_post_ABX)
gut_post_ABX <- filter_low_abundance(gut_post_ABX, 0.5*1e-2)

# simulations

sim_ABX <- read.csv("C:/Users/USER/Desktop/Antibiotics/sim_results/perturbed_state.csv", header = FALSE)
sim_ABX <- as.matrix(sim_ABX)
sim_post_ABX <- read.csv("C:/Users/USER/Desktop/Antibiotics/sim_results/filtered_post_perturbed_state.csv", header = FALSE)
sim_post_ABX <- as.matrix(sim_post_ABX)
sim_post_ABX_shuffled <- read.csv("C:/Users/USER/Desktop/Antibiotics/sim_results/shuffled_state.csv", header = FALSE)
sim_post_ABX_shuffled <- as.matrix(sim_post_ABX_shuffled)

# Results recovery

recovery_results <- list()

results_jaccard <- pairwise_jaccard_similarity(recovery_ABX, recovery_post_ABX)

x_recovery_jaccard <- results_jaccard$"jaccard_ABX"
y_recovery_jaccard <- results_jaccard$"jaccard_post_ABX"

recovery_results$jaccard_ABX <- x_recovery_jaccard
recovery_results$jaccard_post_ABX <- y_recovery_jaccard

p_recovery_jaccard <- create_correlation_plot(recovery_results$jaccard_ABX,
                                              recovery_results$jaccard_post_ABX,
                                              0,0,"Jaccard ABX",
                                              "Jaccard post-ABX")

recovery_results$cor_recovery_jaccard <- cor(recovery_results$jaccard_ABX,
                                             recovery_results$jaccard_post_ABX,
                                             method = "pearson")
test_recovery_jaccard <- cor.test(recovery_results$jaccard_ABX,
                                  recovery_results$jaccard_post_ABX,
                                  method = "pearson")
recovery_results$p_value_jaccard <- test_recovery_jaccard$p.value

## Results gut

gut_results <- list()

results_jaccard <- pairwise_jaccard_similarity(gut_ABX, gut_post_ABX)

x_gut_jaccard <- results_jaccard$"jaccard_ABX"
y_gut_jaccard <- results_jaccard$"jaccard_post_ABX"

indices_jaccard <- y_gut_jaccard != 0 & y_gut_jaccard != 1 &
  !is.na(y_gut_jaccard) & !is.na(x_gut_jaccard) & x_gut_jaccard != 0 &
  x_gut_jaccard != 1

x_gut_jaccard <- x_gut_jaccard[indices_jaccard]
y_gut_jaccard <- y_gut_jaccard[indices_jaccard]

gut_results$jaccard_ABX <- x_gut_jaccard
gut_results$jaccard_post_ABX <- y_gut_jaccard

#ind = gut_results$jaccard_ABX>0.1

#p_gut_jaccard <- create_correlation_plot(gut_results$jaccard_ABX[ind],
#                                         gut_results$jaccard_post_ABX[ind],
#                                         0,0,"Jaccard ABX",
#                                         "Jaccard post-ABX")

p_gut_jaccard <- create_correlation_plot(gut_results$jaccard_ABX,
                                         gut_results$jaccard_post_ABX,
                                         0,0,"Jaccard ABX",
                                         "Jaccard post-ABX")

gut_results$cor_gut_jaccard <- cor(gut_results$jaccard_ABX,
                                   gut_results$jaccard_post_ABX,
                                   method = "pearson")
test_gut_jaccard <- cor.test(gut_results$jaccard_ABX,
                             gut_results$jaccard_post_ABX,
                             method = "pearson")
gut_results$p_value_jaccard <- test_gut_jaccard$p.value

# Results spo

spo_results <- list()

results_jaccard <- pairwise_jaccard_similarity(spo_ABX, spo_post_ABX)

x_spo_jaccard <- results_jaccard$"jaccard_ABX"
y_spo_jaccard <- results_jaccard$"jaccard_post_ABX"

spo_results$jaccard_ABX <- x_spo_jaccard
spo_results$jaccard_post_ABX <- y_spo_jaccard

p_spo_jaccard <- create_correlation_plot(spo_results$jaccard_ABX,
                                         spo_results$jaccard_post_ABX,
                                         0,0,"Jaccard ABX",
                                         "Jaccard post-ABX")

spo_results$cor_spo_jaccard <- cor(spo_results$jaccard_ABX,
                                   spo_results$jaccard_post_ABX,
                                   method = "pearson")
test_spo_jaccard <- cor.test(spo_results$jaccard_ABX,
                             spo_results$jaccard_post_ABX,
                             method = "pearson")
spo_results$p_value_jaccard <- test_spo_jaccard$p.value

grid.arrange(p_spo_jaccard, p_recovery_jaccard, p_gut_jaccard,
             ncol = 1, nrow = 3)

###### simulations #########

## real results

sim_results <- list()

results_jaccard <- pairwise_jaccard_similarity(sim_ABX, sim_post_ABX)

x_sim_jaccard <- results_jaccard$"jaccard_ABX"
y_sim_jaccard <- results_jaccard$"jaccard_post_ABX"

sim_results$jaccard_ABX <- x_sim_jaccard
sim_results$jaccard_post_ABX <- y_sim_jaccard

p_sim_jaccard <- create_correlation_plot(sim_results$jaccard_ABX,
                                         sim_results$jaccard_post_ABX,
                                         0.03,0.03,"Jaccard ABX",
                                         "Jaccard post-ABX")

sim_results$cor_sim_jaccard <- cor(sim_results$jaccard_ABX,
                                   sim_results$jaccard_post_ABX,
                                   method = "pearson")
test_sim_jaccard <- cor.test(sim_results$jaccard_ABX,
                             sim_results$jaccard_post_ABX,
                             method = "pearson")
sim_results$p_value_jaccard <- test_sim_jaccard$p.value

## shuffled results

sim_shuffled_results <- list()

results_jaccard <- pairwise_jaccard_similarity(sim_ABX, sim_post_ABX_shuffled)

x_sim_shuffled_jaccard <- results_jaccard$"jaccard_ABX"
y_sim_shuffled_jaccard <- results_jaccard$"jaccard_post_ABX"

sim_shuffled_results$jaccard_ABX <- x_sim_shuffled_jaccard
sim_shuffled_results$jaccard_post_ABX <- y_sim_shuffled_jaccard

p_sim_shuffled_jaccard <- create_correlation_plot(sim_shuffled_results$jaccard_ABX,
                                                  sim_shuffled_results$jaccard_post_ABX,
                                                  0.03,0.03,"Jaccard ABX",
                                                  "Jaccard post-ABX")

sim_shuffled_results$cor_sim_shuffled_jaccard <- cor(sim_shuffled_results$jaccard_ABX,
                                                     sim_shuffled_results$jaccard_post_ABX,
                                                     method = "spearman")
test_sim_shuffled_jaccard <- cor.test(sim_shuffled_results$jaccard_ABX,
                                      sim_shuffled_results$jaccard_post_ABX,
                                      method = "spearman")
sim_shuffled_results$p_value_jaccard <- test_sim_shuffled_jaccard$p.value

grid.arrange(p_sim_jaccard, p_sim_shuffled_jaccard, ncol = 1, nrow = 2)



#pairwise_nonzero_common_elements <- function(mat) {
# Number of rows in the matrix
#  n <- nrow(mat)

# Initialize a vector to store the pairwise counts
#  pairwise_counts <- numeric()

# Iterate over all unique row pairs
#  for (i in 1:(n-1)) {
#    for (j in (i+1):n) {
# Find indices of nonzero elements in each row
#      nonzero_i <- which(mat[i, ] != 0)
#      nonzero_j <- which(mat[j, ] != 0)

# Count common nonzero elements
#      common_nonzero <- length(intersect(nonzero_i, nonzero_j))

# Append the count to the vector
#      pairwise_counts <- c(pairwise_counts, common_nonzero)
#    }
#  }

#  return(pairwise_counts)
#}

#pairwise_jaccard_similarity <- function(mat) {
# Number of rows in the matrix
#  n <- nrow(mat)

# Initialize a vector to store the Jaccard similarities
#  jaccard_similarities <- numeric()
#  intersections <- list()

#  for (i in 1:(n-1)) {
#    for (j in (i+1):n) {
# Find indices of nonzero elements in each row
#      nonzero_i <- which(mat[i, ] != 0)
#      nonzero_j <- which(mat[j, ] != 0)

# Calculate intersection and union for presence/absence
#      inersection_ind <- intersect(nonzero_i, nonzero_j)
#      intersection <- length(inersection_ind)
#      union <- length(union(nonzero_i, nonzero_j))

# Calculate Jaccard similarity
#      similarity <- intersection / union

# Append the similarity to the vector
#      jaccard_similarities <- c(jaccard_similarities, similarity)
#      intersections[[length(intersections) + 1]] <- inersection_ind
#    }
#  }
#  return(list("jaccard" = jaccard_similarities, "inter" = intersections))
#}

#pairwise_jaccard_similarity_filter <- function(mat, indices) {
# Number of rows in the matrix
#  n <- nrow(mat)

# Initialize a vector to store the Jaccard similarities
#  jaccard_similarities_filtered <- numeric()

#  for (i in 1:(n-1)) {
#    for (j in (i+1):n) {
# Find indices of nonzero elements in each row
#      nonzero_i <- which(mat[i, ] != 0)
#      nonzero_j <- which(mat[j, ] != 0)

# Calculate intersection and union for presence/absence
#      inersection_ind <- intersect(nonzero_i, nonzero_j)
#      intersection <- length(inersection_ind)
#      union <- length(union(nonzero_i, nonzero_j))

# Calculate Jaccard similarity
#      similarity <- intersection / union

# Append the similarity to the vector#
#      jaccard_similarities <- c(jaccard_similarities, similarity)
#      intersections[[length(intersections) + 1]] <- inersection_ind
#    }
#  }
#  return(list("jaccard" = jaccard_similarities, "inter" = intersections))
#}

#pairwise_bray_curtis <- function(mat) {
# Number of rows in the matrix
#  n <- nrow(mat)

# Initialize a vector to store the Bray-Curtis distances
#  distances <- numeric()

#  for (i in 1:(n-1)) {
#    for (j in (i+1):n) {

# Bray-Curtis calculation
#      sum_min <- sum(pmin(mat[i, ], mat[j, ]))
#      sum_total <- sum(mat[i, ]) + sum(mat[j, ])
#      distance <- 1 - (2 * sum_min / sum_total)

#      distances <- c(distances, distance)
#    }
#  }

#  return(distances)
#}

#gut_results <- list()

#x_gut <- pairwise_bray_curtis(gut_post_ABX) 
#y_gut <- pairwise_jaccard_similarity(gut_ABX)

#indices <- x_gut != 1

#x_gut <- x_gut[indices]
#y_gut <- y_gut[indices]

#indices <- y_gut != 0

#x_gut <- x_gut[indices]
#y_gut <- y_gut[indices]

#gut_results$x_gut <- x_gut
#gut_results$y_gut <- y_gut

#p_gut <- create_correlation_plot(gut_results$y_gut,
#                                 gut_results$x_gut)

#gut_results$cor_gut <- cor(gut_results$x_gut,
#                           gut_results$y_gut,
#                           method = "spearman")
#test_gut <- cor.test(gut_results$x_gut,
#                     gut_results$y_gut,
#                     method = "spearman")
#gut_results$p_value <- test_gut$p.value

#pairwise_bray_curtis <- function(mat) {
#  #Number of rows in the matrix
#  n <- nrow(mat)
#  
#  #Initialize a vector to store the Bray-Curtis distances
#  distances <- numeric()
#  
#  for (i in 1:(n-1)) {
#    for (j in (i+1):n) {
#      
#      #Bray-Curtis calculation
#      sum_min <- sum(pmin(mat[i, ], mat[j, ]))
#      sum_total <- sum(mat[i, ]) + sum(mat[j, ])
#      distance <- 1 - (2 * sum_min / sum_total)
#      
#      distances <- c(distances, distance)
#    }
#  }
#  
#  return(distances)
#}
