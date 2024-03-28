library(ggplot2)
library(gridExtra)
library(grid)

# Recovery data
SR_recovery <- c(0.09523809523809523, 0.6268656716417911, 0.6883116883116883,
                  0.36923076923076925, 0.1564625850340136, 0.2732919254658385,
                  0.6015037593984962, 0.39766081871345027, 0.366412213740458)
D_recovery <- c(0.6943919900071144, 0.6576755399645706, 0.6111901462598519,
                0.7539708731053667, 0.6902174936990457, 0.8272625649938046,
                0.5944160123289969, 0.5276614500358653, 0.6321068399629699)

pearson_cor <- cor(SR_recovery, D_recovery, method = "pearson")
print(paste("Pearson correlation coefficient:", pearson_cor))
pearson_test <- cor.test(SR_recovery, D_recovery, method = "pearson")
print(paste("Pearson correlation p-value:", pearson_test$p.value))
spearman_cor <- cor(SR_recovery, D_recovery, method = "spearman")
print(paste("Spearman correlation coefficient:", spearman_cor))
spearman_test <- cor.test(SR_recovery, D_recovery, method = "spearman")
print(paste("Spearman correlation p-value:", spearman_test$p.value))

p_recovery <- ggplot(data = data.frame(SR_recovery, D_recovery),
                     aes(x = SR_recovery, y = D_recovery)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") + 
  theme_minimal() +
  labs(x = "Species Ratio (ABX/Baseline)",
       y = "Distance from baseline (BC)") +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16),
        axis.line = element_line(color = "black"))  

# Gut data
SR_gut <- c(0.6541353383458647, 0.7156862745098039, 0.18947368421052632,
            0.6902173913043478, 0.8095238095238095, 0.55, 0.39090909090909093,
            0.6782006920415224, 0.6627218934911243, 0.8881118881118881,
            0.6741573033707865, 0.7184466019417476, 0.7754010695187166,
            0.5992217898832685, 0.7642276422764228, 0.35772357723577236,
            0.8702290076335878, 0.09036144578313253, 0.7413793103448276)
D_gut <- c(0.30173144448480443, 0.35040561811357307, 0.7376195665334787,
           0.4487226056423296, 0.45344472696452354, 0.41603099648867903,
           0.486983896355491, 0.420753117810873, 0.31250756750211894,
           0.4245065988618477, 0.5849376437825402, 0.3782540259111272,
           0.30282116478992616, 0.44472696452355004, 0.34786293740162244,
           0.7159462404649474, 0.29918876377285386, 0.5545465552730355,
           0.38200750696210195)

pearson_cor <- cor(SR_gut, D_gut, method = "pearson")
print(paste("Pearson correlation coefficient:", pearson_cor))
pearson_test <- cor.test(SR_gut, D_gut, method = "pearson")
print(paste("Pearson correlation p-value:", pearson_test$p.value))
spearman_cor <- cor(SR_gut, D_gut, method = "spearman")
print(paste("Spearman correlation coefficient:", spearman_cor))
spearman_test <- cor.test(SR_gut, D_gut, method = "spearman")
print(paste("Spearman correlation p-value:", spearman_test$p.value))

p_gut <- ggplot(data = data.frame(SR_gut, D_gut),
                     aes(x = SR_gut, y = D_gut)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") + 
  theme_minimal() +
  labs(x = "Species Ratio (ABX/Baseline)",
       y = "Distance from baseline (BC)") +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16),
        axis.line = element_line(color = "black"))


# Eran Elinav
SR_EE <- c(0.5256064690026954, 0.29473684210526313, 0.5060240963855422,
           0.6747720364741641, 0.17539863325740318, 0.11737089201877934,
           0.43956043956043955)
D_EE <- c(0.16010536698942424, 0.16939780801519388, 0.12058060932661546,
          0.16026549554189046, 0.36687968536859006, 0.36878867001855137,
          0.34374546807050366)

pearson_cor <- cor(SR_EE, D_EE, method = "pearson")
print(paste("Pearson correlation coefficient:", pearson_cor))
pearson_test <- cor.test(SR_EE, D_EE, method = "pearson")
print(paste("Pearson correlation p-value:", pearson_test$p.value))
spearman_cor <- cor(SR_EE, D_EE, method = "spearman")
print(paste("Spearman correlation coefficient:", spearman_cor))
spearman_test <- cor.test(SR_EE, D_EE, method = "spearman")
print(paste("Spearman correlation p-value:", spearman_test$p.value))

#SR_EE <- c(0.5256064690026954, 0.29473684210526313, 0.5060240963855422,
#                 0.6747720364741641, 0.17539863325740318, 0.11737089201877934,
#                 0.43956043956043955)
#D_EE <- c(0.125346312955821, 0.12199793633360434, 0.058553012699553375,
#                0.15991830841886945, 0.3603506527044332, 0.34444141995442606,
#                0.35085440728577094)

#SR_EE <- c(0.5256064690026954, 0.29473684210526313, 0.5060240963855422,
#                 0.6747720364741641, 0.17539863325740318, 0.11737089201877934,
#                 0.43956043956043955)
#D_EE <- c(0.11668415882078054, 0.12853593130626262, 0.10316974213865587,
#                0.13770937602831868, 0.349177554174149, 0.35202704502571475,
#                0.3328157308955728)

p_EE <- ggplot(data = data.frame(SR_EE, D_EE),
                aes(x = SR_EE, y = D_EE)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") + 
  theme_minimal() +
  labs(x = "Species Ratio (ABX/Baseline)",
       y = "Distance from baseline (BC)") +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16),
        axis.line = element_line(color = "black"))

grid_layout <- arrangeGrob(p_recovery, p_gut, p_EE, ncol = 3, nrow = 1)

grid.draw(grid_layout)

grid.text("a", x = unit(0.028, "npc") - unit(1, "lines"),
          y = unit(0.975, "npc"), just = "left", gp = gpar(fontsize = 25))
grid.text("b", x = unit(0.358, "npc") - unit(1, "lines"),
          y = unit(0.975, "npc"), just = "left", gp = gpar(fontsize = 25))
grid.text("c", x = unit(0.688, "npc") - unit(1, "lines"),
          y = unit(0.975, "npc"), just = "left", gp = gpar(fontsize = 25))

