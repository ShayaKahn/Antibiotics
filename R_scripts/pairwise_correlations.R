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
      
      if (length(inersection_ind_ABX) == 0){
          # Find indices of nonzero elements in each row in post_ABX_cohort
          nonzero_i_post_ABX <- which(post_ABX_cohort[i, ] != 0)
          nonzero_j_post_ABX <- which(post_ABX_cohort[j, ] != 0)
        } else {
          # Find indices of nonzero elements in each row in post_ABX_cohort
          nonzero_i_post_ABX <- which(post_ABX_cohort[i, ][-inersection_ind_ABX] != 0)
          nonzero_j_post_ABX <- which(post_ABX_cohort[j, ][-inersection_ind_ABX] != 0)
        }
      
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


create_correlation_plot <- function(x, y, width, height, x_title, y_title,
                                    main_title, alpha) {
  # Function to generate the plots
  p <- ggplot(data = data.frame(x = x, y = y),
              aes(x = x, y = y)) +
    #geom_point(alpha = 0.3) +  
    geom_jitter(width = width, height = height, alpha = alpha) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1.2) +  
    theme_minimal() +
    labs(x = x_title,
         y = y_title,
         title = main_title) +  
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),  
          axis.text.x = element_text(size = 14),  
          axis.text.y = element_text(size = 14),  
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16),
          axis.line = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5, size = 20))  
  return(p)
}

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
gut_post_ABX_shuffled <- read.csv("C:/Users/USER/Desktop/Antibiotics/Gut/Data/shuffled_gut.csv",
                                    header = FALSE)
gut_post_ABX_shuffled <- as.matrix(gut_post_ABX_shuffled)
gut_post_ABX_shuffled <- filter_low_abundance(gut_post_ABX_shuffled, 0.5*1e-2)

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

## Results gut

# Real

gut_results <- list()

results_jaccard <- pairwise_jaccard_similarity(gut_ABX, gut_post_ABX)

x_gut_jaccard <- results_jaccard$"jaccard_ABX"
y_gut_jaccard <- results_jaccard$"jaccard_post_ABX"

#indices_jaccard <- y_gut_jaccard != 0 & y_gut_jaccard != 1 &
#  !is.na(y_gut_jaccard) & !is.na(x_gut_jaccard) & x_gut_jaccard != 0 &
#  x_gut_jaccard != 1

#x_gut_jaccard <- x_gut_jaccard[indices_jaccard]
#y_gut_jaccard <- y_gut_jaccard[indices_jaccard]

gut_results$jaccard_ABX <- x_gut_jaccard
gut_results$jaccard_post_ABX <- y_gut_jaccard

p_gut_jaccard <- create_correlation_plot(gut_results$jaccard_ABX,
                                         gut_results$jaccard_post_ABX,
                                         0,0,"Jaccard ABX",
                                         "Jaccard post-ABX",
                                         "Real data",1)

gut_results$cor_gut_jaccard <- cor(gut_results$jaccard_ABX,
                                   gut_results$jaccard_post_ABX,
                                   method = "pearson")
test_gut_jaccard <- cor.test(gut_results$jaccard_ABX,
                             gut_results$jaccard_post_ABX,
                             method = "pearson")
gut_results$p_value_jaccard <- test_gut_jaccard$p.value

# Shuffled

gut_Shuffled_results <- list()

results_Shuffled_jaccard <- pairwise_jaccard_similarity(gut_ABX, gut_post_ABX_shuffled)

x_gut_Shuffled_jaccard <- results_Shuffled_jaccard$"jaccard_ABX"
y_gut_Shuffled_jaccard <- results_Shuffled_jaccard$"jaccard_post_ABX"

#indices_jaccard <- y_gut_jaccard != 0 & y_gut_jaccard != 1 &
#  !is.na(y_gut_jaccard) & !is.na(x_gut_jaccard) & x_gut_jaccard != 0 &
#  x_gut_jaccard != 1

#x_gut_jaccard <- x_gut_jaccard[indices_jaccard]
#y_gut_jaccard <- y_gut_jaccard[indices_jaccard]

gut_Shuffled_results$jaccard_ABX <- x_gut_Shuffled_jaccard
gut_Shuffled_results$jaccard_post_ABX <- y_gut_Shuffled_jaccard

p_gut_Shuffled_jaccard <- create_correlation_plot(gut_Shuffled_results$jaccard_ABX,
                                                  gut_Shuffled_results$jaccard_post_ABX,
                                                  0,0,"Jaccard ABX",
                                                  "Jaccard post-ABX",
                                                  "Shuffled real data",1)

gut_Shuffled_results$cor_gut_Shuffled_jaccard <- cor(gut_Shuffled_results$jaccard_ABX,
                                                     gut_Shuffled_results$jaccard_post_ABX,
                                                     method = "pearson")
test_gut_Shuffled_jaccard <- cor.test(gut_Shuffled_results$jaccard_ABX,
                                      gut_Shuffled_results$jaccard_post_ABX,
                                      method = "pearson")
gut_Shuffled_results$p_value_jaccard <- test_gut_Shuffled_jaccard$p.value

grid.arrange(p_gut_jaccard, p_gut_Shuffled_jaccard, ncol = 2, nrow = 1)

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
                                         "Jaccard post-ABX",
                                         "Simulated data",0.3)

sim_results$cor_sim_jaccard <- cor(sim_results$jaccard_ABX,
                                   sim_results$jaccard_post_ABX,
                                   method = "spearman")
test_sim_jaccard <- cor.test(sim_results$jaccard_ABX,
                             sim_results$jaccard_post_ABX,
                             method = "spearman")
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
                                                  "Jaccard post-ABX",
                                                  "Shuffled simulated data",0.3)

sim_shuffled_results$cor_sim_shuffled_jaccard <- cor(sim_shuffled_results$jaccard_ABX,
                                                     sim_shuffled_results$jaccard_post_ABX,
                                                     method = "spearman")
test_sim_shuffled_jaccard <- cor.test(sim_shuffled_results$jaccard_ABX,
                                      sim_shuffled_results$jaccard_post_ABX,
                                      method = "spearman")
sim_shuffled_results$p_value_jaccard <- test_sim_shuffled_jaccard$p.value

grid.arrange(p_sim_jaccard, p_sim_shuffled_jaccard, ncol = 2, nrow = 1)

# spearman correlation vs interaction strength

spearman_corrs <- read.csv("C:/Users/USER/Desktop/Antibiotics/sim_results/spearman_corrs.csv", header = FALSE)
spearman_corrs <- as.vector(spearman_corrs)
interaction_strength_vector <- read.csv("C:/Users/USER/Desktop/Antibiotics/sim_results/interaction_strength_vector.csv", header = FALSE)
interaction_strength_vector<- as.vector(interaction_strength_vector)
