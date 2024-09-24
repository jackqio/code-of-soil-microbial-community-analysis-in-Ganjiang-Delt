
rm(list=ls())
#
library(car)

#
OTU_abundance <- as.data.frame(t(read.table("otu table.txt")))
Env_value <- read.table("en fac plot.txt", header = TRUE, row.names = 1)

# CY1:31,HH32:68
env_value <- Env_value[,8:length(Env_value)]

# Extract the row name of the environmental factor data frame
env_row_names <- rownames(env_value)

# Filtering the OTU table based on row names of environmental factors
filtered_OTU_abundance <- OTU_abundance[env_row_names, ]
otu_pca <- as.data.frame(t(filtered_OTU_abundance))

df <- env_value

# multiple linear regression analysis (MLRA)
model <- lm(filtered_OTU_abundance[,1] ~ . , data = df)

# Calculating the VIF value
vif_values <- vif(model)

# Print VIF values for each environmental factor
print(vif_values)

# Check VIF values and remove factors with high colinearity
high_vif_threshold <- 5
while (any(vif_values > high_vif_threshold)) {
  # Find the factor corresponding to the largest VIF value
  max_vif_factor <- names(vif_values)[which.max(vif_values)]
  
  # Remove the factor with the largest VIF value from the data frame
  df <- df[, !(names(df) %in% max_vif_factor)]
  
  # Rerunning multiple linear regression analysis and VIF calculations
  model <- lm(filtered_OTU_abundance[,1] ~ ., data = df)
  vif_values <- vif(model)
  
  # 
  print("")
  print(vif_values)
}


#### According to the results of VIF
####PCA----
library(ggplot2)
library(ggbiplot)
library(ggnewscale)
library(dplyr)

# Delete all 0 columns
df <- df[, colSums(df != 0) > 0]

#The original data were z-score normalized
dt<-as.matrix(scale(df))
# normalization
# normalizing function
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
# Applying the Normalization Function
normalized_df <- as.data.frame(lapply(df, normalize))




# Computation of principal components
pca <- prcomp(dt)
# View the results of the principal component analysis
summary(pca)

# Extract PC Score；
dt1<-pca$x
head(dt1)
# Add grouping information;
dt1<-data.frame(dt1,Env_value$group)
head(dt1)

# plot
pcaxy <- as.data.frame(pca$x[,1:2])
pcaxy$avde <- Env_value$avde
pcaxy$Group1 <- Env_value$Group1
pcaaro <- data.frame(pc1=pca$rotation[,1])
pcaaro$pc2 <- pca$rotation[,2]
pcaaro$lab <- row.names(pcaaro)
# Calculate the mean and standard deviation for each subgroup
group_summary <- pcaxy %>%
  group_by(Group1) %>%
  summarize(mean_PC1 = mean(PC1), mean_PC2 = mean(PC2),
            sd_PC1 = sd(PC1), sd_PC2 = sd(PC2))

s_pca <- summary(pca)
cap1_exp <- paste('PCA1:', round(s_pca$importance[2,1]*100, 2), '%')
cap2_exp <- paste('PCA2:', round(s_pca$importance[2,2]*100, 2), '%')


p <- ggplot() +
  geom_point(data=subset(pcaxy,Group1=='CY'), 
             mapping=aes(PC1, PC2, color=subset(pcaxy,Group1=='CY')[['avde']], size=subset(pcaxy,Group1=='CY')[['avde']])) + 
  scale_color_gradient('CY', low = "#ffbfb8", high = "#8D0B25") +
  new_scale("color") +
  geom_point(data=subset(pcaxy,Group1=='HH'),
             mapping=aes(PC1, PC2, color=subset(pcaxy,Group1=='HH')[['avde']], size=subset(pcaxy,Group1=='HH')[['avde']])) + 
  scale_color_gradient('HH', low = "#CCFAFF", high = "#0E606B") + 
  theme_bw() +
  geom_segment(data = pcaaro, 
               aes(x = 0, xend = pc1*5, y = 0, yend = pc2*5), 
               arrow = arrow(length = unit(0.02, "npc")), 
               alpha = 1) +
  # geom_text(data = pcaaro,
  #           aes(x = pc1*10, y = pc2*10, label = lab),
  #           hjust = 0.5, vjust = 0.5,
  #           nudge_x = 0.01, nudge_y = 0.01,
  #           size = 3.5) + 
  labs(x=cap1_exp, y=cap2_exp) +
  stat_ellipse(data=pcaxy,
               geom = "polygon", level=0.9,
               linetype = 1, size=0.8,
               aes(x = PC1, y = PC2, fill=Group1),
               alpha=0.1, show.legend = T) +
  scale_size_continuous(name="avde", range = c(4, 10))
p

dt2 <- data.frame(Env_value[,1:6], dt)

#PERMANOVA
library(vegan)
#
dt_pa <- data.frame(dt2[,7:length(dt2)], group = dt2$group2)

#计算PE
res_permanova <- adonis2(dt_pa[,-length(dt_pa)]~dt_pa$group, 
                         data = dt_pa[,-length(dt_pa)], 
                         method = "euclidean")
res_permanova

## statistical significance of difference between PC1 and PC2
p_value <- adonis2(pca$x[,1:2] ~ Env_value$Group1, method = "euclidean")
p_value


####PCA statisitic test##
# 1. Load the necessary R packages
library(FactoMineR)  # 
library(factoextra)  # 

# 2. Standardized data
dt_scaled <- dt

# 3. Perform principal component analysis
pca_result <- PCA(dt_scaled, graph = FALSE)

# 4. View the results of the principal component analysis
# View the loadings of each variable (component loadings)
pca_result$var$coord  # 查看所有变量在主成分上的投影

# # 5. Visualization of PCA results
# # Loading plot for the variable
# fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# 6. Obtained variance contribution of each principal component (explained variance)
pca_result$eig  # Shows the eigenvalues of the principal components and the percentage of variance explained by them

# 7. Evaluating component load significance using Bootstrap 
# Load the boot package for boot law verification
library(boot)

# Define bootstrap sampling function to obtain the first principal component loadings for each PCA
boot_pca <- function(data, indices) {
  data_sample <- data[indices, ]  # Generate new data based on bootstrap samples
  pca_res <- PCA(scale(data_sample), graph = FALSE)  # Implementation of PCA
  return(pca_res$var$coord[, 1])  # Returns the loading of the first principal component
}

# Execute Bootstrap
set.seed(123)
boot_results <- boot(data = dt_scaled, statistic = boot_pca, R = 500)

# Calculate the 95% confidence interval for the load
boot.ci(boot_results, type = "perc")







####2.Get confidence intervals for each factor####
# Assuming your dataframe is named dt and has been normalized to dt_scaled

# 
library(FactoMineR)
library(factoextra)
library(boot)

dt_scaled <- dt
# 
pca_result <- PCA(dt_scaled, graph = FALSE)


# 
boot_pca <- function(data, indices) {
  data_sample <- data[indices, ]
  pca_res <- PCA(data_sample, graph = FALSE)
  return(pca_res$var$coord[, 1])  # Return the loading of the first principal component, modify it yourself if needed
}

# 
set.seed(123)
boot_results <- boot(data = dt_scaled, statistic = boot_pca, R = 500)

# 
p <- ncol(dt_scaled)
variable_names <- colnames(dt_scaled)

# Create a data frame to store the confidence interval results
ci_df <- data.frame(
  Variable = variable_names,
  Lower = numeric(p),
  Upper = numeric(p),
  Significant = logical(p)
)

# Loop to obtain confidence intervals for each variable and determine significance
for (i in 1:p) {
  ci <- boot.ci(boot_results, type = "bca", index = i)#The ci$*** below will have to be changed in response to the parameter change.
  ci_df$Lower[i] <- ci$bca[4]  # lower limit
  ci_df$Upper[i] <- ci$bca[5]  # limit
  # Determine whether the confidence interval contains zero
  ci_df$Significant[i] <- !(ci_df$Lower[i] <= 0 & ci_df$Upper[i] >= 0)
}

# 查看结果
print(ci_df)
