# Covid data clustering:
# required packages:
library(cluster)
library(factoextra)
library(gridExtra)
library(clustMixType)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(GGally)
library(readr)
library(caret)
library(lattice)
library(FeatureImpCluster)
library(gmp)
library(ClustImpute)
library(tidyverse)
library(magrittr)
library(attempt)
library(rmcorr)

#Loading the theme 

my_theme <- function(base_size = 10, base_family = "sans"){
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      axis.title = element_text(size = 10),
      panel.grid.major = element_line(color = "gray"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#f7fdff"),
      strip.background = element_rect(fill = "#001d60", color = "#00113a", size =0.5),
      strip.text = element_text(face = "bold", size = 10, color = "white"),
      legend.position = "bottom",
      legend.justification = "center",
      legend.background = element_blank(),
      panel.border = element_rect(color = "grey5", fill = NA, size = 0.5)
    )
}

theme_set(my_theme())
clust_colmap = c("#f7286d","#1faae0","#ffbf1f")

# glopal function to remove NAs
bar_missing <- function(x){
  require(reshape2)
  x %>%
    is.na %>%
    melt %>%
    ggplot(data = .,
           aes(x = Var2)) +
    geom_bar(aes(y=(..count..),fill=value),alpha=0.7)+
    scale_fill_manual(values=c("skyblue","red"),
                      name = "",
                      labels = c("Available","Missing"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle=45, vjust=0.5)) +
    labs(x = "Variables in Dataset",
         y = "Observations")+coord_flip()
}
# 

#Games.howell Post-hoc tests function to compare the output clusters
games.howell <- function(grp, obs) {
  
  #Create combinations
  combs <- combn(unique(grp), 2)
  
  # Statistics that will be used throughout the calculations:
  # n = sample size of each group
  # groups = number of groups in data
  # Mean = means of each group sample
  # std = variance of each group sample
  n <- tapply(obs, grp, length)
  groups <- length(tapply(obs, grp, length))
  Mean <- tapply(obs, grp, mean)
  std <- tapply(obs, grp, var)
  
  statistics <- lapply(1:ncol(combs), function(x) {
    
    mean.diff <- Mean[combs[2,x]] - Mean[combs[1,x]]
    
    #t-values
    t <- abs(Mean[combs[1,x]] - Mean[combs[2,x]]) / sqrt((std[combs[1,x]] / n[combs[1,x]]) + (std[combs[2,x]] / n[combs[2,x]]))
    
    # Degrees of Freedom
    df <- (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]])^2 / # Numerator Degrees of Freedom
      ((std[combs[1,x]] / n[combs[1,x]])^2 / (n[combs[1,x]] - 1) + # Part 1 of Denominator Degrees of Freedom 
         (std[combs[2,x]] / n[combs[2,x]])^2 / (n[combs[2,x]] - 1)) # Part 2 of Denominator Degrees of Freedom
    
    #p-values
    p <- ptukey(t * sqrt(2), groups, df, lower.tail = FALSE)
    
    # Sigma standard error
    se <- sqrt(0.5 * (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]]))
    
    # Upper Confidence Limit
    upper.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff + qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Lower Confidence Limit
    lower.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff - qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Group Combinations
    grp.comb <- paste(combs[1,x], ':', combs[2,x])
    
    # Collect all statistics into list
    stats <- list(grp.comb, mean.diff, se, t, df, p, upper.conf, lower.conf)
  })
  
  # Unlist statistics collected earlier
  stats.unlisted <- lapply(statistics, function(x) {
    unlist(x)
  })
  
  # Create dataframe from flattened list
  results <- data.frame(matrix(unlist(stats.unlisted), nrow = length(stats.unlisted), byrow=TRUE))
  
  # Select columns set as factors that should be numeric and change with as.numeric
  results[c(2, 3:ncol(results))] <- round(as.numeric(as.matrix(results[c(2, 3:ncol(results))])), digits = 3)
  
  # Rename data frame columns
  colnames(results) <- c('groups', 'Mean Difference', 'Standard Error', 't', 'df', 'p', 'upper limit', 'lower limit')
  
  return(results)
}

#
# read and load the data
setwd("~/Downloads/covid-data")
# read the data
covid <- readxl::read_excel("Metadata-Workbook.xlsx")
covid <- as.data.frame(covid)
attach(covid)

# testing homogeneity of variability:
# response variable one by one:
response_factor <- covid$`MIP-3α`

# predictor variable
time.point <- as.factor(covid$`Date Post Positive`)

# individual ID
subject.ID <- as.factor(covid$`Patient ID`)

# dataframe
my.dataframe <- data.frame(subject.ID ,time.point,response_factor)
neo <- matrix(my.dataframe$response_factor, nrow = 159, ncol = 3)
mauchly.test(lm(neo ~ 1), X = ~ 1)


# using rmcorr to identify the correlations and acounting for the repeated measures:

corrr<- rmcorr(
  covid$`Patient ID`,
  `IL-6`,
  `MIP-3α`,
  dataset = covid,
  CI.level = 0.95,
  CIs = c("analytic", "bootstrap"),
  nreps = 100,
  bstrap.out = F
)
## S3 method for class 'rmc'
plot(
  corrr,
  overall = F,
  palette = NULL,
  xlab = NULL,
  ylab = NULL,
  overall.col = "gray60",
  overall.lwd = 3,
  overall.lty = 2)

#
corrr<- rmcorr(
  covid$`Patient ID`,
  `IL-23`,
  `IL-12p70`,
  dataset = covid,
  CI.level = 0.95,
  CIs = c("analytic", "bootstrap"),
  nreps = 100,
  bstrap.out = F
)

plot(
  corrr,
  overall = F,
  palette = NULL,
  xlab = NULL,
  ylab = NULL,
  overall.col = "gray60",
  overall.lwd = 3,
  overall.lty = 2)

#explore the missing values 
bar_missing(covid)
#scale  only the numeric data:
scaled_df = covid %>%as.data.frame()
scaled_df[,c(2, 17:59)] <- scale(covid[,c(2, 17:59)])
head(scaled_df)

# set the cols name for plot later:
cols <- c("Sex", "CHF", "Diabetes", "Renal", "OSA", "Hypertension",  "Asthma")
scaled_df %<>% mutate_at(cols, factor)
str(scaled_df)
# prepare the matrix to use it in the k-prototypes clustering:
X_mat = scaled_df %>% select(-c(1)) # remove the label ID
set.seed(771)
# choose the optimal k using the total within sum of square method:

Es <- numeric(10)

for(i in 1:10){
  kpres <- kproto(X_mat, 
                  k = i, nstart = 25, 
                  lambda = lambdaest(X_mat),
                  verbose = FALSE)
  Es[i] <- kpres$tot.withinss}

# plot the optimal k based on effect size
tibble(Cluster = c(1:10), Es = Es) %>% 
  ggplot(aes(x = Cluster, y = Es)) + 
  geom_point(size = 3, 
             col ="red3") +
  geom_path() + 
  geom_vline(xintercept = 3, 
             linetype = 2)+
  scale_x_continuous(breaks = c(1:10))

# calculate the the optimal number of clusters based on the choosen index for k-Prototype clustering using ptbiserial and silhouette methods
k_opt <- validation_kproto(method = "ptbiserial", data = X_mat, k = 2:10, nstart = 5)
#plot the silhouette result
tibble(Cluster = c(2:10), 
       Metric = as.vector(k_opt$indices)) %>% 
  ggplot(aes(x = Cluster, 
             y = Metric)) + 
  geom_point(size = 3, 
             col ="red3") +
  geom_path() + 
  geom_vline(xintercept = 3, 
             linetype = 2)+
  scale_x_continuous(breaks = c(2:10))
#
total_withinss <- c()

for (i in 1:8) {
  kproto <- clustMixType::kproto(X_mat,
                                 k = i,
                                 nstart = 25)
  total_withinss[i] <- kproto$tot.withinss
}
total_withinss
# plot the elbow method
tibble(k = 1:length(total_withinss),
       total_error = total_withinss
) %>%
  ggplot(aes(x = k,
             y = total_error)
  ) +
  geom_point(size = 2) +
  geom_line() +
  theme_bw() +
  labs(x = "Number of Clusters",
       y = "tot.withinss") +
  geom_text(x = 3,
            y = total_withinss[3],
            label = "ELBOW",
            alpha = 0.5,
            color = "blue",
            size = 5)

set.seed(771)
#clustering the data using k=3
kpres_3 = kproto(x = X_mat,
                 k = 3,
                 lambda = lambdaest(X_mat))

covid <- na.omit(covid)
dim(covid) #157  59, 1 observation(s) was removed 
# model validation:

valid_df = covid %>% mutate(Cluster = as.factor( kpres_3$cluster))
# explore the distribution of the categorical factors among the clusters:
valid_df %>%gather(Age,CHF,Hypertension,Diabetes,
                   key = "Para",value="Value")%>%
  ggplot(aes(x=Value, fill = Cluster))+
  geom_density(alpha=0.5,col="black")+
  facet_wrap(~Para,ncol=2,scales = "free")+
  scale_fill_manual(values=clust_colmap )
# explore the distribution of age/IL-6/the categorical factors among the clusters:
valid_df %>%gather(Age,CHF,Hypertension,Diabetes, `IL-6`, Renal,
                   key = "Para",value="Value")%>%
  ggplot(aes(x = Cluster, y=Value, fill = Cluster))+
  geom_boxplot(alpha=0.5,col="black")+
  facet_wrap(~Para,ncol=2,scales = "free")+
  coord_flip()+
  scale_fill_manual(values=clust_colmap )

## explore the distribution of the Analytes across age:
par(mfcol=c(2,2))
valid_df %>% ggplot(aes(x=Age,y=SAP))+
  stat_density2d(geom="polygon",
                 aes(fill=Cluster,
                     col = Cluster,
                     alpha = ..level..))+
  geom_jitter(color="black",size=2, alpha = 0.3)+
  scale_fill_manual(values=clust_colmap)+
  scale_color_manual(values=clust_colmap)

#################################################
