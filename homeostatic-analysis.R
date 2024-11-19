#setting working directory 
setwd("/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis") #set working directory

#loading packages 
library(dplyr)
library(MicrogliaMorphologyR)
library(factoextra)
library(ppclust)
library(randomForest)
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

# Here, we are creating a new data frame that contains the first 2 PCs and original dataset, then renaming the data_kmeans$cluster column to simply say "Cluster". You can bind together as many of the PCs as you want. Binding the original, untransformed data is useful if you want to plot the raw values of any individual morphology measures downstream. 
pca_kmeans <- cbind(pca_data[1:2], combined_data, as.data.frame(data_kmeans$cluster)) %>%
  rename(Cluster=`data_kmeans$cluster`)

#display cluster features in a heat map + plot 
clusterfeatures(pca_kmeans, featurestart=9, featureend=35)
plot <- clusterplots(pca_kmeans, "PC1", "PC2")

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
                                                          Cluster=="1" ~ "Rod-Like")
                                                )


#Ensure 'Sex' and 'Treatment' are factors
stats_input <- cp
stats_input$MouseID <- factor(stats_input$MouseID)
stats_input$Cluster <- factor(stats_input$Cluster)
stats_input$Treatment <- factor(stats_input$Treatment)
stats_output <- stats_cluster.animal(data = stats_input, 
                                       model = "percentage ~ Sex*Cluster*Treatment", 
                                       posthoc1 = "~Treatment|Cluster|Sex", 
                                       posthoc2 = "~Sex|Cluster|Treatment")
stats_output[[2]]
stats_output[[3]]


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

#Performing individual stats on variables of interest
morphology_data <- combined_data %>% 
  group_by(Sex, Treatment, MouseID) %>% 
  summarise(across("Foreground pixels":"Maximum branch length", ~mean(.x))) %>% 
  gather(Measure, Value, "Foreground pixels":"Maximum branch length")

morphology_stats_input <- morphology_data 
morphology_stats_input$Treatment <- factor(morphology_stats_input$Treatment)
morphology_stats_input$Sex <- factor(morphology_stats_input$Sex)
# run stats analysis for changes in individual morphology measures
# you can specify up to two posthoc comparisons (posthoc1 and posthoc2 arguments) - if you only have one set of posthocs to run, specify the same comparison twice for both arguments. you will just get the same results in output[[2]] and output[[3]].
morphology_stats_input_filtered <- morphology_stats_input %>% filter(Measure %in% c("Foreground pixels",
                       "Average branch length"))
View(morphology_stats_input_filtered)
morphology_stats_testing_filtered <- stats_morphologymeasures.animal(data = morphology_stats_input_filtered, 
                                                 model = "Value ~ Sex*Treatment", type="lm",
                                                 posthoc1 = "~Sex|Treatment", 
                                                 posthoc2 = "~Treatment|Sex")

combined_data_gathered <- combined_data %>% gather(measure, value, 7:ncol(combined_data))
female_data <- combined_data_gathered %>% filter(Sex == "F")
male_data <- combined_data_gathered %>% filter(Sex == "M")
outliers_distributions(female_data)
outliers_distributions(male_data)









# #
# 
# #filtering out only control and ramified cells
# control_cluster_four <- pca_kmeans[, -c(1, 2)] %>% filter(Treatment == "PBS", Cluster == "4")
# 
# control_cluster_four <- control_cluster_four[ , -ncol(control_cluster_four)]
# View(control_cluster_four)
# pcadata_elbow(control_cluster_four, featurestart=7, featureend=33)
# ramified_normalized <-transform_log(control_cluster_four, 1, start=7, end=33) 
# #calculating pca data
# pca_data_ramified <- pcadata(ramified_normalized, featurestart=7, featureend=33,
#                     pc.start=1, pc.end=10)
# 
# #normalizing first 3 PCA's
# pca_data_scale_ramified <- transform_scale(pca_data_ramified, start=1, end=3) 
# #kmeans sampling + kmeans clustering
# kmeans_input_ramified <- pca_data_scale_ramified[1:3]
# sampling_ramified <- kmeans_input_ramified[sample(nrow(kmeans_input_ramified), 5000),] #sample 5000 random rows for cluster optimization
# fviz_nbclust(sampling_ramified, kmeans, method = 'silhouette', nstart=25, iter.max=50) # 4 clusters
# data_kmeans_ramified <- kmeans(kmeans_input_ramified, centers=6)
# pca_kmeans_ramified <- cbind(pca_data_ramified[1:2], control_cluster_four, as.data.frame(data_kmeans_ramified$cluster)) %>%
#   rename(Cluster=`data_kmeans_ramified$cluster`)
# 
# View(pca_kmeans_ramified)
# #display cluster features in a heat map + plot 
# clusterfeatures(pca_kmeans_ramified, featurestart=9, featureend=35)
# plot <- clusterplots(pca_kmeans_ramified, "PC1", "PC2")
# 
# cp_ramified <- clusterpercentage(pca_kmeans_ramified, "Cluster", Sex, Treatment)
# ggplot(cp_ramified, aes(x = Cluster, y = percentage, fill = Sex)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(title = "Percentages of Microglia in Each Category Separated by Sex and Treatment", x = "Cluster", y = "Percentage", fill = "Sex") +
#   theme_minimal() +
#   facet_wrap(~ Treatment)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # 
# # # Function to classify microglia morphology based on relative characteristics
# # classify_microglia_shapes <- function(data) {
# #   # Calculate mean values for each cluster
# #   cluster_means <- data %>%
# #     group_by(Cluster) %>%
# #     summarise(
# #       branches = mean(`# of branches`, na.rm = TRUE),
# #       circularity = mean(Circularity, na.rm = TRUE),
# #       density = mean(`Density of foreground pixels in hull area`, na.rm = TRUE),
# #       span_ratio = mean(`Span ratio of hull (major/minor axis)`, na.rm = TRUE),
# #       avg_branch_length = mean(`Average branch length`, na.rm = TRUE),
# #       max_span = mean(`Maximum span across hull`, na.rm = TRUE)
# #     )
# #   # Determine ranks for each feature to assign relative morphologies
# #   cluster_means <- cluster_means %>%
# #     mutate(
# #       branches_rank = dense_rank(branches),
# #       circularity_rank = dense_rank(desc(circularity)),  # Higher circularity ranks higher
# #       density_rank = dense_rank(desc(density)),
# #       span_ratio_rank = dense_rank(desc(span_ratio)),
# #       avg_branch_length_rank = dense_rank(desc(avg_branch_length)),
# #       max_span_rank = dense_rank(desc(max_span))
# #     ) %>%
# #     # Assign Morphology based on relative characteristics
# #     mutate(Morphology = case_when(
# #       # Ameboid: Low branch count, high circularity and density
# #       branches_rank == min(branches_rank) & 
# #         circularity_rank == max(circularity_rank) & 
# #         density_rank == max(density_rank) ~ "Ameboid",
# #       
# #       # Rod-like: Moderate branch count, elongated shape with high span ratio and moderate circularity
# #       branches_rank < max(branches_rank) & 
# #         span_ratio_rank == max(span_ratio_rank) &
# #         circularity_rank < max(circularity_rank) ~ "Rod-like",
# #       
# #       # Ramified: High branch count, low density, and low circularity
# #       branches_rank == max(branches_rank) & 
# #         density_rank == min(density_rank) & 
# #         circularity_rank == min(circularity_rank) ~ "Ramified",
# #       
# #       # Hypertrophic: Higher branch length and span, moderate to high branch count
# #       avg_branch_length_rank == max(avg_branch_length_rank) & 
# #         max_span_rank == max(max_span_rank) & 
# #         branches_rank >= median(branches_rank) ~ "Hypertrophic",
# #       
# #       # Default for clusters that don’t match specific criteria
# #       TRUE ~ "Unclassified"
# #     ))
# #   
# #   # Join the morphology assignments back to the original data
# #   data <- data %>%
# #     left_join(cluster_means %>% select(Cluster, Morphology), by = "Cluster")
# #   
# #   return(data)
# # }
# 
# # #Introduction
# # #Goal of the project/research question
# # #Achieve it - normalize, group, and compare using blah blah blah 
# # #I used these code segments to do this 
# # #Dont go through code line by line, summarize it - don't get lost in the lines in code
# # #always think about what that means 
# # #plots are good!!
# # #next steps 
# # 
# # #Ameboid only
# # # Extracting the "Ameboid" cluster
# # 
# # ramified_cluster <- filter(pca_kmeans, Cluster == "2")
# # ramified_cluster <- ramified_cluster[, -c(1, 2)]
# # ramified_cluster <- ramified_cluster[, !(names(ramified_cluster) %in% c("Cluster"))]
# # normalized_ramified_cluster <- transform_log(ramified_cluster, 1, start=7, end=33) 
# # #plotting the dataset to explore PCA values 
# # pcadata_elbow(normalized_ramified_cluster, featurestart=7, featureend=33)
# # pca_data_ramified <- pcadata(normalized_ramified_cluster, featurestart=7, featureend=33,
# #                     pc.start=1, pc.end=10)
# # pca_data_scale_ramified <- transform_scale(pca_data_ramified, start=1, end=4)
# # kmeans_input_ramified <- pca_data_scale_ramified[1:4]
# # 
# # sampling_ramified <- kmeans_input_ramified[sample(nrow(kmeans_input_ramified), 5000),] #sample 5000 random rows for cluster optimization
# # 
# # fviz_nbclust(sampling_ramified, kmeans, method = 'wss', nstart=25, iter.max=50) # 5 clusters
# # fviz_nbclust(sampling_ramified, kmeans, method = 'silhouette', nstart=25, iter.max=50) # 5 clusters
# # 
# # # cluster and combine with original data
# # data_kmeans_ramified <- kmeans(kmeans_input_ramified, centers=5)
# # 
# # # Here, we are creating a new data frame that contains the first 2 PCs and original dataset, then renaming the data_kmeans$cluster column to simply say "Cluster". You can bind together as many of the PCs as you want. Binding the original, untransformed data is useful if you want to plot the raw values of any individual morphology measures downstream. 
# # pca_kmeans_ramified <- cbind(pca_data_ramified[1:2], ramified_cluster, as.data.frame(data_kmeans_ramified$cluster)) %>%
# #   rename(Cluster=`data_kmeans_ramified$cluster`)
# # 
# # clusterfeatures(pca_kmeans_ramified, featurestart=9, featureend=33)
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # # Selecting relevant columns: total area, branch area, and number of branches (update column names if needed)
# # ameboid_data <- ameboid_cluster %>%
# #   select(MouseID, Sex, TotalArea = "Area", BranchArea = "Average branch length", NumBranches = "# of branches")
# # 
# # # Summary statistics by Sex for visual inspection
# # summary(ameboid_data)
# # 
# # # Performing ANOVA to assess differences between sexes for each feature
# # # Total Area
# # total_area_anova <- aov(TotalArea ~ Sex, data = ameboid_data)
# # summary(total_area_anova)
# # 
# # # Branch Area
# # branch_area_anova <- aov(BranchArea ~ Sex, data = ameboid_data)
# # summary(branch_area_anova)
# # 
# # # Number of Branches
# # num_branches_anova <- aov(NumBranches ~ Sex, data = ameboid_data)
# # summary(num_branches_anova)
# # 
# # # Post-hoc tests (if ANOVA shows significance)
# # if (summary(total_area_anova)[[1]][["Pr(>F)"]][1] < 0.05) {
# #   posthoc_total_area <- TukeyHSD(total_area_anova)
# #   print(posthoc_total_area)
# # }
# # 
# # if (summary(branch_area_anova)[[1]][["Pr(>F)"]][1] < 0.05) {
# #   posthoc_branch_area <- TukeyHSD(branch_area_anova)
# #   print(posthoc_branch_area)
# # }
# # 
# # if (summary(num_branches_anova)[[1]][["Pr(>F)"]][1] < 0.05) {
# #   posthoc_num_branches <- TukeyHSD(num_branches_anova)
# #   print(posthoc_num_branches)
# # }
# # 
# # #Ramified only
# # # Extracting the "Ramified" cluster
# # ramified_cluster <- filter(pca_kmeans, Cluster == "2")
# # 
# # # Selecting relevant columns: total area, branch area, and number of branches (update column names if needed)
# # ramified_data <- ramified_cluster %>%
# #   select(Sex, TotalArea = "Area", BranchArea = "Average branch length", NumBranches = "# of branches")
# # 
# # # Summary statistics by Sex for visual inspection
# # summary(ramified_data)
# # 
# # # Performing ANOVA to assess differences between sexes for each feature
# # # Total Area
# # total_area_anova_ramified <- aov(TotalArea ~ Sex, data = ramified_data)
# # summary(total_area_anova_ramified)
# # View(total_area_anova_ramified)
# # 
# # # Branch Area
# # branch_area_anova_ramified <- aov(BranchArea ~ Sex, data = ramified_data)
# # summary(branch_area_anova_ramified)
# # 
# # # Number of Branches
# # num_branches_anova_ramified <- aov(NumBranches ~ Sex, data = ramified_data)
# # summary(num_branches_anova_ramified)
# # 
# # 
# # # Post-hoc tests (if ANOVA shows significance)
# # if (summary(total_area_anova_ramified)[[1]][["Pr(>F)"]][1] < 0.05) {
# #   posthoc_total_area_ramified <- TukeyHSD(total_area_anova_ramified)
# #   print(posthoc_total_area_ramified)
# # }
# # 
# # if (summary(branch_area_anova_ramified)[[1]][["Pr(>F)"]][1] < 0.05) {
# #   posthoc_branch_area_ramified <- TukeyHSD(branch_area_anova_ramified)
# #   print(posthoc_branch_area_ramified)
# # }
# # 
# # if (summary(num_branches_anova_ramified)[[1]][["Pr(>F)"]][1] < 0.05) {
# #   posthoc_num_branches_ramified <- TukeyHSD(num_branches_anova_ramified)
# #   print(posthoc_num_branches_ramified)
# # }
# # #prep the code that we are going to present 
# # #one ttest or anova (review that!), which is the independent and which is the dependent variable 
# # #for loop!! 
# # 
