library(ggplot2)
library(readxl)
library(vegan)
library(gridExtra)

########## Functions ##########

pairwise_num_species <- function(ABX_cohort, post_ABX_cohort) {
  # Number of rows in the matrix
  n <- nrow(ABX_cohort)
  
  # Initialize a vector to store the Jaccard similarities
  n_sim_ABX <- numeric()
  n_sim_post_ABX <- numeric()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Find indices of nonzero elements in each row in ABX_cohort
      nonzero_i_ABX <- which(ABX_cohort[i, ] != 0)
      nonzero_j_ABX <- which(ABX_cohort[j, ] != 0)
      
      # Calculate intersection and union for presence/absence in ABX_cohort
      inersection_ind_ABX <- intersect(nonzero_i_ABX, nonzero_j_ABX)
      
      num_species_ABX <- length(inersection_ind_ABX) 
      
      # Append the similarity to the vector
      n_sim_ABX <- c(n_sim_ABX, num_species_ABX)
      
      nonzero_i_post_ABX <- which(post_ABX_cohort[i, ][-inersection_ind_ABX] != 0)
      nonzero_j_post_ABX <- which(post_ABX_cohort[j, ][-inersection_ind_ABX] != 0)
      
      # Calculate intersection and union for presence/absence in ABX_cohort
      inersection_ind_post_ABX <- intersect(nonzero_i_post_ABX, nonzero_j_post_ABX)
      
      num_species_post_ABX <- length(inersection_ind_post_ABX) / length(post_ABX_cohort[i, ][-inersection_ind_ABX])
      
      # Append the similarity to the vector
      n_sim_post_ABX <- c(n_sim_post_ABX, num_species_post_ABX)
      
    }
  }
  return(list("species_ABX" = n_sim_ABX,
              "species_post_ABX" = n_sim_post_ABX))
  
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
}

########## Load data ##########

# Eran Elinav

spo_baseline <- read_excel("C:/Users/USER/Desktop/Antibiotics/Eran Elinav/Data/Spo_baseline_cohort.xlsx",
                           col_names = FALSE)
spo_baseline <- as.matrix(spo_baseline)
spo_baseline <- spo_baseline[, -1]
spo_ABX <- read_excel("C:/Users/USER/Desktop/Antibiotics/Eran Elinav/Data/Spo_ABX_cohort.xlsx",
                      col_names = FALSE)
spo_ABX <- as.matrix(spo_ABX)
spo_ABX <- spo_ABX[, -1]
spo_post_ABX <- read_excel("C:/Users/USER/Desktop/Antibiotics/Eran Elinav/Data/Spo_post_ABX_cohort.xlsx",
                           col_names = FALSE)
spo_post_ABX <- as.matrix(spo_post_ABX)
spo_post_ABX <- spo_post_ABX[, -1]


# Recovery

recovery_baseline <- read_excel("C:/Users/USER/Desktop/Antibiotics/Recovery/Data/baseline_rel_abund_rarefied_appear_4.xlsx",
                                col_names = FALSE)
recovery_baseline <- as.matrix(recovery_baseline)
recovery_ABX <- read_excel("C:/Users/USER/Desktop/Antibiotics/Recovery/Data/rel_abund_rarefied_4.xlsx",
                           col_names = FALSE)
recovery_ABX <- as.matrix(recovery_ABX)
recovery_post_ABX <- read_excel("C:/Users/USER/Desktop/Antibiotics/Recovery/Data/rel_abund_rarefied_180_appear_4.xlsx",
                                col_names = FALSE)
recovery_post_ABX <- as.matrix(recovery_post_ABX)

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
column_sums_gut_ABX <- colSums(gut_ABX)
normalized_gut_ABX <- sweep(gut_ABX, 2, column_sums_gut_ABX, "/")
gut_ABX = t(normalized_gut_ABX)
column_sums_gut_post_ABX <- colSums(gut_post_ABX)
normalized_gut_post_ABX <- sweep(gut_post_ABX, 2, column_sums_gut_post_ABX, "/")
gut_post_ABX = t(normalized_gut_post_ABX)

# simulations

sim_ABX <- read.csv("C:/Users/USER/Desktop/Antibiotics/sim_results/perturbed_state.csv", header = FALSE)
sim_ABX <- as.matrix(sim_ABX)
sim_post_ABX <- read.csv("C:/Users/USER/Desktop/Antibiotics/sim_results/filtered_post_perturbed_state.csv", header = FALSE)
sim_post_ABX <- as.matrix(sim_post_ABX)
sim_post_ABX_shuffled <- read.csv("C:/Users/USER/Desktop/Antibiotics/sim_results/shuffled_state.csv", header = FALSE)
sim_post_ABX_shuffled <- as.matrix(sim_post_ABX_shuffled)

# Results recovery

recovery_results <- list()

results_species <- pairwise_num_species(recovery_ABX, recovery_post_ABX)

x_recovery_species <- results_species$"species_ABX"
y_recovery_species <- results_species$"species_post_ABX"

recovery_results$species_ABX <- x_recovery_species
recovery_results$species_post_ABX <- y_recovery_species

p_recovery_species <- create_correlation_plot(recovery_results$species_ABX,
                                              recovery_results$species_post_ABX,
                                              0,0,"# Common species ABX",
                                              "# Common species post-ABX")

recovery_results$cor_recovery_species <- cor(recovery_results$species_ABX,
                                             recovery_results$species_post_ABX,
                                             method = "spearman")
test_recovery_species <- cor.test(recovery_results$species_ABX,
                                  recovery_results$species_post_ABX,
                                  method = "spearman")
recovery_results$p_value_species <- test_recovery_species$p.value

## Results gut

gut_results <- list()

results_species <- pairwise_num_species(gut_ABX, gut_post_ABX)

x_gut_species <- results_species$"species_ABX"
y_gut_species <- results_species$"species_post_ABX"

indices_species <- y_gut_species != 0

x_gut_species <- x_gut_species[indices_species]
y_gut_species <- y_gut_species[indices_species]

gut_results$species_ABX <- x_gut_species
gut_results$species_post_ABX <- y_gut_species

p_gut_species <- create_correlation_plot(gut_results$species_ABX,
                                         gut_results$species_post_ABX,
                                         0,0,"# Common species ABX",
                                         "# Common species post-ABX")

gut_results$cor_gut_species <- cor(gut_results$species_ABX,
                                   gut_results$species_post_ABX,
                                   method = "spearman")
test_gut_species <- cor.test(gut_results$species_ABX,
                             gut_results$species_post_ABX,
                             method = "spearman")
gut_results$p_value_species <- test_gut_species$p.value

# Results spo

spo_results <- list()

results_species <- pairwise_num_species(spo_ABX, spo_post_ABX)

x_spo_species <- results_species$"species_ABX"
y_spo_species <- results_species$"species_post_ABX"

spo_results$species_ABX <- x_spo_species
spo_results$species_post_ABX <- y_spo_species

p_spo_species <- create_correlation_plot(spo_results$species_ABX,
                                         spo_results$species_post_ABX,
                                         0,0,"# Common species ABX",
                                         "# Common species post-ABX")

spo_results$cor_spo_species <- cor(spo_results$species_ABX,
                                   spo_results$species_post_ABX,
                                   method = "spearman")
test_spo_species <- cor.test(spo_results$species_ABX,
                             spo_results$species_post_ABX,
                             method = "spearman")
spo_results$p_value_species <- test_spo_species$p.value

grid.arrange(p_spo_species, p_recovery_species, p_gut_species,
             ncol = 1, nrow = 3)

###### simulations #########

## real results

sim_results <- list()

results_species <- pairwise_num_species(sim_ABX, sim_post_ABX)

x_sim_species <- results_species$"species_ABX"
y_sim_species <- results_species$"species_post_ABX"

sim_results$species_ABX <- x_sim_species
sim_results$species_post_ABX <- y_sim_species

p_sim_species <- create_correlation_plot(sim_results$species_ABX,
                                         sim_results$species_post_ABX,
                                         0,0,"# Common species ABX",
                                         "# Common species post-ABX")

sim_results$cor_sim_species <- cor(sim_results$species_ABX,
                                   sim_results$species_post_ABX,
                                   method = "spearman")
test_sim_species <- cor.test(sim_results$species_ABX,
                             sim_results$species_post_ABX,
                             method = "spearman")
sim_results$p_value_species <- test_sim_species$p.value

## shuffled results

results_species <- pairwise_num_species(sim_ABX, sim_post_ABX_shuffled)

x_sim_shuffled_species <- results_species$"species_ABX"
y_sim_shuffled_species <- results_species$"species_post_ABX"

sim_shuffled_results$species_ABX <- x_sim_shuffled_species
sim_shuffled_results$species_post_ABX <- y_sim_shuffled_species

p_sim_shuffled_species <- create_correlation_plot(sim_shuffled_results$species_ABX,
                                                  sim_shuffled_results$species_post_ABX,
                                                  0,0,"# Common species ABX",
                                                  "# Common species post-ABX")

sim_shuffled_results$cor_sim_shuffled_species <- cor(sim_shuffled_results$species_ABX,
                                                     sim_shuffled_results$species_post_ABX,
                                                     method = "spearman")
test_sim_shuffled_species <- cor.test(sim_shuffled_results$species_ABX,
                                      sim_shuffled_results$species_post_ABX,
                                      method = "spearman")
sim_shuffled_results$p_value_species <- test_sim_shuffled_species$p.value

grid.arrange(p_sim_species, p_sim_shuffled_species, ncol = 1, nrow = 2)