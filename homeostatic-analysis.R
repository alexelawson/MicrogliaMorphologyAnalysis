#setting working directory 
setwd("/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis") #set working directory

#loading packages 
library(dplyr)
library(MicrogliaMorphologyR)
library(factoextra)
library(ppclust)
set.seed(1)

#5-10lines of your own code 
#pvalue adjustment method, anova, ttest 
#read in the data and know the row data that will be enough FOR NEXT TIME 
#what do I need to check, do I need to normalize/remove outliers etc. 
#saving the directory for the data 

#loading in the two sample datasets
data_1xLPS <- MicrogliaMorphologyR::data_1xLPS_mouse
data_2xLPS <- MicrogliaMorphologyR::data_2xLPS_mouse

#filtering out only the control rows from each dataset
data_1xLPS_control <- filter(data_1xLPS, Treatment == 'PBS')
data_2xLPS_control <- filter(data_2xLPS, Treatment == 'PBS')

#removing the two columns that are extra in the second dataset (for ease of merging later)
data_2xLPS_control_mutated <- data_2xLPS_control[, -c(5, 6)]

#combining the two datasets 
combined_control_data <- rbind(data_1xLPS_control, data_2xLPS_control_mutated)

write.csv(combined_control_data, "/Users/alexlawson/Documents/GitHub/MicrogliaMorphologyAnalysis/control_data.csv", row.names = FALSE)

#performing a log transformation on the relevant columns from the dataset 
normalized_combined_control_data <-transform_log(combined_control_data, 1, start=7, end=33) 

#plotting the dataset to explore PCA values 
pcadata_elbow(normalized_combined_control_data, featurestart=7, featureend=33)

#
pca_data <- pcadata(normalized_combined_control_data, featurestart=7, featureend=33,
                    pc.start=1, pc.end=10)
pca_data_scale <- transform_scale(pca_data, start=1, end=3) # scale pca data as input for k-means clustering

kmeans_input <- pca_data_scale[1:3]
sampling <- kmeans_input[sample(nrow(kmeans_input), 5000),] #sample 5000 random rows for cluster optimization
fviz_nbclust(sampling, kmeans, method = 'silhouette', nstart=25, iter.max=50) # 4 clusters
data_kmeans <- kmeans(kmeans_input, centers=4)

# Here, we are creating a new data frame that contains the first 2 PCs and original dataset, then renaming the data_kmeans$cluster column to simply say "Cluster". You can bind together as many of the PCs as you want. Binding the original, untransformed data is useful if you want to plot the raw values of any individual morphology measures downstream. 
pca_kmeans <- cbind(pca_data[1:2], combined_control_data, as.data.frame(data_kmeans$cluster)) %>%
  rename(Cluster=`data_kmeans$cluster`)

clusterfeatures(pca_kmeans, featurestart=9, featureend=33)

plot <- clusterplots(pca_kmeans, "PC1", "PC2")


clusterfeatures(pca_kmeans, featurestart=9, featureend=35)
#cluster 1 = hypertrophic (average territory span, high branch thickness as explained by pixel density in hull)
#cluster 2 = ramified (largest territory span and branching complexity)
#cluster 3 = ameboid (lowest territory span, high circularity, smallest branch lengths)
#cluster 4 = rod-like (greatest oblongness, lowest circularity)

cp <- clusterpercentage(pca_kmeans, "Cluster", Sex)
cp <- cp %>% mutate(Cluster = 
                      case_when(Cluster=="4" ~ "Rod-like",
                                Cluster=="3" ~ "Ameboid",
                                Cluster=="2" ~ "Ramified",
                                Cluster=="1" ~ "Hypertrophic"))

ggplot(cp, aes(x = Cluster, y = percentage, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparison of Percentages by Gender",
       x = "Cluster",
       y = "Percentage (%)") +
  scale_fill_manual(values = c("F" = "black", "M" = "white"))

stats.input <- 

stats.testing <- stats_cluster.animal(data = stats.input, 
                                      model = "percentage ~ Cluster*Sex", 
                                      posthoc1 = "~Sex|Cluster", 
                                      posthoc2 = "~Sex|Cluster", adjust = "bonferroni")

cp %>% # in this example, we filter for our brain region of interest
  clusterpercentage_boxplots(Sex)





#Introduction
#Goal of the project/research question
#Achieve it - normalize, group, and compare using blah blah blah 
#I used these code segments to do this 
#Dont go through code line by line, summarize it - don't get lost in the lines in code
#always think about what that means 
#plots are good!!
#next steps 

#Ameboid only
# Extracting the "Ameboid" cluster

ramified_cluster <- filter(pca_kmeans, Cluster == "2")
ramified_cluster <- ramified_cluster[, -c(1, 2)]
ramified_cluster <- ramified_cluster[, !(names(ramified_cluster) %in% c("Cluster"))]
normalized_ramified_cluster <- transform_log(ramified_cluster, 1, start=7, end=33) 
#plotting the dataset to explore PCA values 
pcadata_elbow(normalized_ramified_cluster, featurestart=7, featureend=33)
pca_data_ramified <- pcadata(normalized_ramified_cluster, featurestart=7, featureend=33,
                    pc.start=1, pc.end=10)
pca_data_scale_ramified <- transform_scale(pca_data_ramified, start=1, end=4)
kmeans_input_ramified <- pca_data_scale_ramified[1:4]

sampling_ramified <- kmeans_input_ramified[sample(nrow(kmeans_input_ramified), 5000),] #sample 5000 random rows for cluster optimization

fviz_nbclust(sampling_ramified, kmeans, method = 'wss', nstart=25, iter.max=50) # 5 clusters
fviz_nbclust(sampling_ramified, kmeans, method = 'silhouette', nstart=25, iter.max=50) # 5 clusters

# cluster and combine with original data
data_kmeans_ramified <- kmeans(kmeans_input_ramified, centers=5)

# Here, we are creating a new data frame that contains the first 2 PCs and original dataset, then renaming the data_kmeans$cluster column to simply say "Cluster". You can bind together as many of the PCs as you want. Binding the original, untransformed data is useful if you want to plot the raw values of any individual morphology measures downstream. 
pca_kmeans_ramified <- cbind(pca_data_ramified[1:2], ramified_cluster, as.data.frame(data_kmeans_ramified$cluster)) %>%
  rename(Cluster=`data_kmeans_ramified$cluster`)

clusterfeatures(pca_kmeans_ramified, featurestart=9, featureend=33)














# Selecting relevant columns: total area, branch area, and number of branches (update column names if needed)
ameboid_data <- ameboid_cluster %>%
  select(MouseID, Sex, TotalArea = "Area", BranchArea = "Average branch length", NumBranches = "# of branches")

# Summary statistics by Sex for visual inspection
summary(ameboid_data)

# Performing ANOVA to assess differences between sexes for each feature
# Total Area
total_area_anova <- aov(TotalArea ~ Sex, data = ameboid_data)
summary(total_area_anova)

# Branch Area
branch_area_anova <- aov(BranchArea ~ Sex, data = ameboid_data)
summary(branch_area_anova)

# Number of Branches
num_branches_anova <- aov(NumBranches ~ Sex, data = ameboid_data)
summary(num_branches_anova)

# Post-hoc tests (if ANOVA shows significance)
if (summary(total_area_anova)[[1]][["Pr(>F)"]][1] < 0.05) {
  posthoc_total_area <- TukeyHSD(total_area_anova)
  print(posthoc_total_area)
}

if (summary(branch_area_anova)[[1]][["Pr(>F)"]][1] < 0.05) {
  posthoc_branch_area <- TukeyHSD(branch_area_anova)
  print(posthoc_branch_area)
}

if (summary(num_branches_anova)[[1]][["Pr(>F)"]][1] < 0.05) {
  posthoc_num_branches <- TukeyHSD(num_branches_anova)
  print(posthoc_num_branches)
}

#Ramified only
# Extracting the "Ramified" cluster
ramified_cluster <- filter(pca_kmeans, Cluster == "2")

# Selecting relevant columns: total area, branch area, and number of branches (update column names if needed)
ramified_data <- ramified_cluster %>%
  select(Sex, TotalArea = "Area", BranchArea = "Average branch length", NumBranches = "# of branches")

# Summary statistics by Sex for visual inspection
summary(ramified_data)

# Performing ANOVA to assess differences between sexes for each feature
# Total Area
total_area_anova_ramified <- aov(TotalArea ~ Sex, data = ramified_data)
summary(total_area_anova_ramified)
View(total_area_anova_ramified)

# Branch Area
branch_area_anova_ramified <- aov(BranchArea ~ Sex, data = ramified_data)
summary(branch_area_anova_ramified)

# Number of Branches
num_branches_anova_ramified <- aov(NumBranches ~ Sex, data = ramified_data)
summary(num_branches_anova_ramified)


# Post-hoc tests (if ANOVA shows significance)
if (summary(total_area_anova_ramified)[[1]][["Pr(>F)"]][1] < 0.05) {
  posthoc_total_area_ramified <- TukeyHSD(total_area_anova_ramified)
  print(posthoc_total_area_ramified)
}

if (summary(branch_area_anova_ramified)[[1]][["Pr(>F)"]][1] < 0.05) {
  posthoc_branch_area_ramified <- TukeyHSD(branch_area_anova_ramified)
  print(posthoc_branch_area_ramified)
}

if (summary(num_branches_anova_ramified)[[1]][["Pr(>F)"]][1] < 0.05) {
  posthoc_num_branches_ramified <- TukeyHSD(num_branches_anova_ramified)
  print(posthoc_num_branches_ramified)
}
#prep the code that we are going to present 
#one ttest or anova (review that!), which is the independent and which is the dependent variable 
#for loop!! 

