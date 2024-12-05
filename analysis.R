#Microglia-Morphology-Analysis
#Author: Alexandra Lawson
#Dataset taken from the Ciernia Lab
#Description: pca and clustering analysis to investigate microglial morphology and how they respond to immune activation, followed by graphical analysis and linear mixed model stats testing with post hoc analysis

#setting working directory 
setwd("/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis") #set working directory

#loading packages 
library(dplyr)
library(MicrogliaMorphologyR)
library(factoextra)
library(ppclust)
set.seed(1)


#loading in the two sample datasets
data_1xLPS <- MicrogliaMorphologyR::data_1xLPS_mouse
data_2xLPS <- MicrogliaMorphologyR::data_2xLPS_mouse
#removing the two columns that are extra in the second dataset (for ease of merging later)
data_2xLPS_mutated <- data_2xLPS[, -c(5, 6)]
#combining the two datasets 
combined_data <- rbind(data_1xLPS, data_2xLPS_mutated)
#performing a log transformation on the relevant columns from the dataset 
normalized_combined_data <-transform_log(combined_data, 1, start=7, end=33) 


#saving dataset as separate file 
write.csv(combined_data, "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/data.csv", row.names = FALSE)

#plotting the dataset to explore PCA values 
pcadata_elbow(normalized_combined_data, featurestart=7, featureend=33)

#calculating pca data
pca_data <- pcadata(normalized_combined_data, featurestart=7, featureend=33,
                    pc.start=1, pc.end=10)

#normalizing first 3 PCA's
pca_data_scale <- transform_scale(pca_data, start=1, end=3) # scale pca data as input for k-means clustering

#kmeans sampling + kmeans clustering
kmeans_input <- pca_data_scale[1:3]
sampling <- kmeans_input[sample(nrow(kmeans_input), 5000),] #sample 5000 random rows for cluster optimization
fviz_nbclust(sampling, kmeans, method = 'silhouette', nstart=25, iter.max=50) # 4 clusters
data_kmeans <- kmeans(kmeans_input, centers=4)

#Creating a new data frame that contains the first 2 PCs and original dataset, then renaming the data_kmeans$cluster column to simply say "Cluster"
pca_kmeans <- cbind(pca_data[1:2], combined_data, as.data.frame(data_kmeans$cluster)) %>%
  rename(Cluster=`data_kmeans$cluster`)

#display cluster features in a heat map + plot 
clusterfeatures(pca_kmeans_w_cluster, featurestart=9, featureend=35)
plot <- clusterplots(pca_kmeans, "PC1", "PC2")

#view PCA feature correlations and how they explain each morphometric component
pcfeaturecorrelations(pca_data, pc.start=1, pc.end=3, 
                      feature.start=17, feature.end=43, 
                      rthresh=0.75, pthresh=0.05, 
                      title="Correlation between PCs and features")

#calculating percentage of microglia in each cluster and saving in a new dataframe
cp <- clusterpercentage(pca_kmeans, "Cluster", Sex, Treatment, MouseID)
cp$Treatment <- factor(cp$Treatment, levels=c("PBS","2xLPS", "LPS"))

#mutating cluster column to assign classically understood shapes to each group
cp <- cp %>% mutate(Cluster = 
                      case_when(Cluster=="4" ~ "Ramified",
                                Cluster=="3" ~ "Ameboid",
                                Cluster=="2" ~ "Hypertrophic",
                                Cluster=="1" ~ "Rod-Like"))

#mutating pca_kmeans dataframe with cluster assignments for ease of analysis 
pca_kmeans_w_cluster <- pca_kmeans %>% mutate(Cluster = 
                                                case_when(Cluster=="4" ~ "Ramified",
                                                          Cluster=="3" ~ "Ameboid",
                                                          Cluster=="2" ~ "Hypertrophic",
                                                          Cluster=="1" ~ "Rod-Like"))
                                                

#Ensuring relevant columns are factors for statistical analysis 
stats_input <- cp
stats_input$MouseID <- factor(stats_input$MouseID)
stats_input$Cluster <- factor(stats_input$Cluster)
stats_input$Treatment <- factor(stats_input$Treatment)

#linear mixed model with 2 post-hoc outputs 
stats_output <- stats_cluster.animal(data = stats_input, 
                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                       posthoc1 = "~Treatment|Cluster|Sex", 
                                       posthoc2 = "~Sex|Cluster|Treatment", adjust = "holm")


#saving stats outputs 
treatment_cluster_stats <- stats_output[[2]]
write.csv(treatment_cluster_stats, "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/statistical-results/hard-cluster-post-hoc1.csv", row.names = FALSE)
sex_cluster_stats <- stats_output[[3]]
write.csv(sex_cluster_stats, "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/statistical-results/hard-cluster-post-hoc2.csv", row.names = FALSE)


significance_data <- stats_output[[2]]   # Statistical significance information for treatment vs. cluster vs. sex

#New column to flag significant bars
significance_data <- significance_data %>%
  mutate(
    outline_color = ifelse(Significant == "significant", "red", "black")  # Red outline for significant bars, black for others 
  )

#Creating the plot
plot_hard_cluster <- ggplot(significance_data %>% filter(Cluster=="Ramified"), aes(x = interaction(Cluster, contrast), y = estimate, fill = Sex)) +
  geom_bar(
    stat = "identity", 
    position = position_dodge(width = 0.8), 
    aes(color = outline_color)  # Use the outline_color column for the bar borders
  ) +
  geom_errorbar(aes(
    ymin = estimate - SE,
    ymax = estimate + SE
  ),
  position = position_dodge(width = 0.8), width = 0.2) +
  labs(
    title = "Comparison of Treatment Effects by Cluster, Treatment, and Sex",
    x = "Cluster and Treatment Comparison",
    y = "Estimated Difference",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)  #Rotate x-axis labels for readability
  ) +
  scale_fill_manual(
    values = c("F" = "white", "M" = "lightgrey")  #Custom colors for Female and Male
  ) +
  scale_color_identity()  #Use color directly from the outline_color column

# Display the plot
print(plot_hard_cluster)
#Saving the plot 
ggsave(filename = "hard-cluster-plot.png", plot = plot_hard_cluster, path = "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/plots-and-images", width = 8, height = 6, dpi = 300)

#Fuzzy clustering code lines commented out because the process is very time intensive
#The output was saved in an RDS file for later use 
#fuzzy cluster analysis
#data_kmeans <- fcm(kmeans_input, centers=4, nstart=25)
#View(data_kmeans)
#saveRDS(data_kmeans, "fuzzy_clustering.rds")

#Load the fuzzy cluster object from the saved location 
fuzzy_clustering <- readRDS("/Users/alexlawson/Downloads/fuzzy_clustering.rds")
#Extracting the columns with membership percentages
memberships <- fuzzy_clustering$u
membership_df <- as.data.frame(memberships)


#Creating a new data frame that contains the first 2 PCs and original dataset + renaming the data_kmeans$cluster column to simply say "Cluster"
fuzzy_cluster_data <- cbind(pca_data[1:2], combined_data, membership_df, pca_kmeans_w_cluster$Cluster) %>%
  rename(Cluster=`pca_kmeans_w_cluster$Cluster`)

View(fuzzy_cluster_data)

fuzzy_cluster_data_filtered <- fuzzy_cluster_data %>% 
  filter(`Cluster 1` > 0.70|
           `Cluster 2` > 0.70|
           `Cluster 3` > 0.70|
           `Cluster 4` > 0.70)
#Viewing the cluster features for analysis 
clusterfeatures(fuzzy_cluster_data_filtered, featurestart=9, featureend=35)

#Calculating the cluster percentages
cp_fuzzy <- clusterpercentage(fuzzy_cluster_data_filtered, "Cluster", MouseID, Treatment, Sex)
View(cp_fuzzy)


#Linear mixed model for fuzzy data
stats_input_fuzzy <- cp_fuzzy
stats_input_fuzzy$MouseID <- factor(stats_input_fuzzy$MouseID)
stats_input_fuzzy$Cluster <- factor(stats_input_fuzzy$Cluster)
stats_input_fuzzy$Treatment <- factor(stats_input_fuzzy$Treatment)
stats_output_fuzzy <- stats_cluster.animal(data = stats_input_fuzzy, 
                                     model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                     posthoc1 = "~Treatment|Cluster|Sex", 
                                     posthoc2 = "~Sex|Cluster|Treatment", adjust = "holm")

#Saving the output from the fuzzy cluster analysis 
stats_treatment_cluster_fuzzy <- stats_output_fuzzy[[2]]
write.csv(stats_treatment_cluster_fuzzy, "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/statistical-results/fuzzy-cluster-post-hoc1.csv", row.names = FALSE)
stats_sex_fuzzy <- stats_output_fuzzy[[3]]
write.csv(stats_sex_fuzzy, "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/statistical-results/fuzzy-cluster-post-hoc2.csv", row.names = FALSE)



#Ncolumn to flag significant bars
significance_data_fuzzy <- stats_output_fuzzy[[2]] %>%
  mutate(
    outline_color = ifelse(Significant == "significant", "red", "black")  # Red outline for significant bars
  )

#Creating the plot
plot_fuzzy <- ggplot(significance_data_fuzzy, aes(x = interaction(Cluster, contrast), y = estimate, fill = Sex)) +
  geom_bar(
    stat = "identity", 
    position = position_dodge(width = 0.8), 
    aes(color = outline_color)  # Use the outline_color column for the bar borders
  ) +
  geom_errorbar(aes(
    ymin = estimate - SE,
    ymax = estimate + SE
  ),
  position = position_dodge(width = 0.8), width = 0.2) +
  labs(
    title = "Comparison of Treatment Effects by Cluster, Treatment, and Sex",
    x = "Cluster and Treatment Comparison",
    y = "Estimated Difference",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  ) +
  scale_fill_manual(
    values = c("F" = "white", "M" = "lightgrey")  # Custom colors for Female and Male
  ) +
  scale_color_identity()  # Use color directly from the outline_color column

#Display the plot
print(plot_fuzzy)

#Saving the plot 
ggsave(filename = "fuzzy-cluster-plot.png", plot = plot_fuzzy, path = "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/plots-and-images", width = 8, height = 6, dpi = 300)
cluster_column <- pca_kmeans_w_cluster$Cluster

#Creating a plot to visualize the ramified cluster data
significance_data_fuzzy <- stats_output_fuzzy[[2]] %>%
  mutate(
    outline_color = ifelse(Significant == "significant", "red", "black")  # Red outline for significant bars
  )


# Create the plot
plot_fuzzy_ramified <- ggplot(significance_data_fuzzy %>% filter(Cluster == "Ramified"), aes(x = interaction(Cluster, contrast), y = estimate, fill = Sex)) +
  geom_bar(
    stat = "identity", 
    position = position_dodge(width = 0.8), 
    aes(color = outline_color)  # Use the outline_color column for the bar borders
  ) +
  geom_errorbar(aes(
    ymin = estimate - SE,
    ymax = estimate + SE
  ),
  position = position_dodge(width = 0.8), width = 0.2) +
  labs(
    title = "Comparison of Treatment Effects by Cluster, Treatment, and Sex",
    x = "Cluster and Treatment Comparison",
    y = "Estimated Difference",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  ) +
  scale_fill_manual(
    values = c("F" = "white", "M" = "lightgrey")  # Custom colors for Female and Male
  ) +
  scale_color_identity()  # Use color directly from the outline_color column

#Displaying the plot
print(plot_fuzzy_ramified)

#Extracting only the ramified column from the original hard cluster set 
ramified_data <- cbind(normalized_combined_data, cluster_column) %>% filter(cluster_column=="Ramified")
#calculating pca data for the ramified data 
pca_data_ramified <- pcadata(ramified_data, featurestart=7, featureend=33,
                    pc.start=1, pc.end=10)
#normalizing first 4 PCA's
pca_data_scale_ramified <- transform_scale(pca_data_ramified, start=1, end=4) # scale pca data as input for k-means clustering
#kmeans sampling + kmeans clustering for the ramified cluster with silhouette + wss to determine optimal number of clusters 
kmeans_input_ramified <- pca_data_scale_ramified[1:4]
sampling_ramified <- kmeans_input_ramified[sample(nrow(kmeans_input_ramified), 5000),] #sample 5000 random rows for cluster optimization
fviz_nbclust(sampling_ramified, kmeans, method = 'silhouette', nstart=25, iter.max=50) # 5 clusters
fviz_nbclust(sampling_ramified, kmeans, method = 'wss', nstart=25, iter.max=50)
data_kmeans_ramified <- kmeans(kmeans_input_ramified, centers=5)


#Creating a new data frame that contains the first 2 PCs and original dataset, then renaming the data_kmeans$cluster column to simply say "Cluster"
pca_kmeans_ramified <- cbind(pca_data_ramified[1:2], ramified_data, as.data.frame(data_kmeans_ramified$cluster)) %>%
  rename(Cluster=`data_kmeans_ramified$cluster`)

#display cluster features in a heat map + plot 
clusterfeatures(pca_kmeans_ramified, featurestart=9, featureend=35)
plot <- clusterplots(pca_kmeans, "PC1", "PC2")


#Linear mixed model stats for the ramified cluster 
stats_input_ramified <- clusterpercentage(pca_kmeans_ramified, "Cluster", Sex, Treatment, MouseID)
stats_input_ramified$MouseID <- factor(stats_input_ramified$MouseID)
stats_input_ramified$Cluster <- factor(stats_input_ramified$Cluster)
stats_input_ramified$Treatment <- factor(stats_input_ramified$Treatment)
stats_output_ramified <- stats_cluster.animal(data = stats_input_ramified, 
                                           model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                           posthoc1 = "~Treatment|Cluster|Sex", 
                                           posthoc2 = "~Treatment|Cluster")


#Saving the ramified cluster stats 
anova_stats_ramified <- stats_output_ramified[[1]]
write.csv(anova_stats_ramified, "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/statistical-results/anova-stats-ramified-cluster.csv", row.names = FALSE)
treatment_cluster_stats_ramified <- stats_output_ramified[[2]]
write.csv(treatment_cluster_stats_ramified, "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/statistical-results/hard-cluster-ramified-post-hoc1.csv", row.names = FALSE)
sex_cluster_stats_ramified <- stats_output_ramified[[3]]
write.csv(sex_cluster_stats_ramified, "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/statistical-results/hard-cluster-ramified-post-hoc2.csv", row.names = FALSE)

#New column to flag significant bars
significance_data_ramified <- stats_output_ramified[[2]] %>%
  mutate(
    outline_color = ifelse(Significant == "significant", "red", "black")  # Red outline for significant bars
  )

#Creating the plot 
plot_ramified_all <- ggplot(significance_data_ramified, aes(x = interaction(Cluster, contrast), y = estimate, fill = Sex)) +
  geom_bar(
    stat = "identity", 
    position = position_dodge(width = 0.8), 
    aes(color = outline_color)  # Use the outline_color column for the bar borders
  ) +
  geom_errorbar(aes(
    ymin = estimate - SE,
    ymax = estimate + SE
  ),
  position = position_dodge(width = 0.8), width = 0.2) +
  labs(
    title = "Comparison of Treatment Effects by Cluster, Treatment, and Sex",
    x = "Cluster and Treatment Comparison",
    y = "Estimated Difference",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  ) +
  scale_fill_manual(
    values = c("F" = "white", "M" = "lightgrey")  # Custom colors for Female and Male
  ) +
  scale_color_identity()  # Use color directly from the outline_color column

#Displaying the plot
print(plot_ramified_all)
#Saving the plot 
ggsave(filename = "ramified_all.png", plot = plot_ramified_all, path = "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/plots-and-images", width = 8, height = 6, dpi = 300)


#New plot to just visualize the clusters of interest 
#New column to flag significant bars and reorder levels
significance_data_ramified <- stats_output_ramified[[3]] %>%
  mutate(
    outline_color = ifelse(Significant == "significant", "red", "black"),  # Red outline for significant bars
    fill_color = ifelse(Significant == "significant", "lightblue", "gray"), # Different fill colors for significance
    Cluster_Contrast = interaction(Cluster, contrast)  # Combine Cluster and contrast
  ) %>%
  arrange(Cluster, contrast) %>%  # Order data by Cluster and then contrast
  mutate(
    Cluster_Contrast = factor(Cluster_Contrast, levels = unique(Cluster_Contrast))  # Set factor levels
  )

#Creating the plot 
plot_ramified <- ggplot(significance_data_ramified, aes(x = Cluster_Contrast, y = estimate)) +
  geom_bar(
    stat = "identity", 
    position = position_dodge(width = 0.8), 
    aes(color = outline_color, fill = fill_color)  # Add fill aesthetic
  ) +
  geom_errorbar(aes(
    ymin = estimate - SE,
    ymax = estimate + SE
  ),
  position = position_dodge(width = 0.8), width = 0.2) +
  labs(
    title = "Comparison of Treatment Effects by Cluster, Treatment",
    x = "Cluster and Treatment Comparison",
    y = "Estimated Difference"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  ) +
  scale_color_identity() +  # Use color directly from the outline_color column
  scale_fill_identity()     # Use fill directly from the fill_color column

#Printing the plot
print(plot_ramified)
#Saving the plot 
ggsave(filename = "ramified_1_and_2_plot.png", plot = plot_ramified, path = "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/plots-and-images", width = 8, height = 6, dpi = 300)




#Supervised cluster analysis purely for investigating purposes 
########################
#Supervised-analysis
normalized_supervised_data <-transform_log(combined_data, 1, start=7, end=33) 
normalized_supervised_data$Sex <- as.factor(normalized_supervised_data$Sex)
normalized_supervised_data$Group <- as.factor(normalized_supervised_data$Treatment)

#calculating pca data
pca_data_supervised <- pcadata(normalized_supervised_data, featurestart=7, featureend=33,
                    pc.start=1, pc.end=10)

#normalizing first 3 PCA's
pca_data_normalized_supervised<- transform_scale(pca_data, start=1, end=3) 

#Modifying column names so that they work in the rfmodel
kmeans_result_supervised <- kmeans(pca_data_normalized_supervised[1:3], centers = 2)
table(kmeans_result_supervised$cluster, normalized_supervised_data$Sex)
colnames(normalized_supervised_data) <- gsub(" ", "_", colnames(normalized_supervised_data))
colnames(normalized_supervised_data) <- gsub("[\\(/)]", "_", colnames(normalized_supervised_data))
colnames(normalized_supervised_data) <- gsub("'", "", colnames(normalized_supervised_data))
colnames(normalized_supervised_data) <- gsub("#", "num", colnames(normalized_supervised_data))

#Random-forests supervised analysis 
rf_model <- randomForest(Sex ~ . - Antibody - MouseID - UniqueID - ID, 
                         data = normalized_supervised_data, 
                         importance = TRUE, 
                         ntree = 500)

#Exploring random forest variables
print(rf_model)
varImpPlot(rf_model)

# Foreground pixel density
ggplot(normalized_supervised_data, aes(x = Treatment, y = Density_of_foreground_pixels_in_hull_area, fill=Sex)) +
  geom_boxplot() +
  labs(title = "Distribution of Foreground Pixel Density by Sex",
       x = "Treatment", y = "Density_of_foreground_pixels_in_hull_area", fill="Sex") +
  theme_minimal()

# average branch length 
ggplot(normalized_supervised_data, aes(x = Treatment, y = Average_branch_length, fill=Sex)) +
  geom_boxplot() +
  labs(title = "Distribution of Foreground Pixel Density by Sex",
       x = "Treatment", y = "Average_branch_length", fill="Sex") +
  theme_minimal()

#plotting M vs. F cluster percentages in experimental and treatment groups 
ggplot(cp, aes(x = Cluster, y = percentage, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Percentages of Microglia in Each Category Separated by Sex and Treatment", x = "Cluster", y = "Percentage", fill = "Sex") +
  theme_minimal() +
  facet_wrap(~ Treatment)

