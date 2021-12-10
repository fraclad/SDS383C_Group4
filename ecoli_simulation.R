library(dbscan)
library(speccalt)
library(dirichletprocess)
library(mvnfast)
library(mclust)
library(NbClust)
library(tsne)
library(scatterplot3d)
library(ggplot2)

tiff("test.tiff", units="in", width=5, height=5, res=300)

data("iris")
x <- as.matrix(iris[, 1:4])
db <- dbscan(x, eps = .4, minPts = 4)

kNNdistplot(x, k = 5 - 1)
abline(h=.7, col = "red", lty=2)

ecoli <- read.table("./data/ecoli.data")

# e coli dataset without labels
data <- ecoli[, 2:8]
# labels
labels <- ecoli[,9]

# kmeans
km = kmeans(data, centers = 8)
km_assignments <- km$cluster
km_rand <- adjustedRandIndex(labels, km_assignments)

# dbscan
kNNdistplot(data, k = 8)
db <- dbscan(data, eps = 0.3, minPts = 3)
db_assignments <- db$cluster
db_rand <- adjustedRandIndex(labels, db_assignments)

# speccalt
kern <- local.rbfdot(data)
s <- speccalt(kern) # I think these are the labels of the clusters it assigns
s_rand <- adjustedRandIndex(labels, s)

# dpmg
model <- DirichletProcessMvnormal(data)

# fgm
source("./CODE_original.R")
fgm <- FG_mixture(data, K = 10, M=25, Iter=500, alpha=1, mm=10)
index_fg = rep(1,nrow(data)) # cluster assignments 
# assignments are based on what sphere it is assigned 
L = ncol(fgm$inclusion_matrix)
for(j in 1:L) {
  index_fg[which(fgm$inclusion_matrix[,j]==1)]=j
}
fg_rand <- adjustedRandIndex(labels, index_fg)

results <- data.frame(c("k-means", "dbscan", "speccalt", "FG"),
                      c(km_rand, db_rand, s_rand, fg_rand))
colnames(results) <- c("Method", "RandIndex")
order <- c("dbscan", "k-means", "speccalt", "FG")

ggplot() + 
  geom_col(data=results, aes(x=Method, y=RandIndex)) +
  scale_x_discrete(limits = order) +
  ylab("adjusted Rand Index") +
  ggtitle("e. coli protein clustering")

ggsave("ecoli_bar.jpeg")

leaf <- read.csv("./leaf/leaf.csv")

features <- leaf[, 3:16]
labels <- leaf[,1]

# k means
km = kmeans(features, centers = 36)
km_assignments <- km$cluster
km_rand <- adjustedRandIndex(labels, km_assignments)

# dbscan
db <- dbscan(features, eps = 0.3, minPts = 3)
db_assignments <- db$cluster
db_rand <- adjustedRandIndex(labels, db_assignments)

# speccalt
kern <- local.rbfdot(features)
s <- speccalt(kern) # I think these are the labels of the clusters it assigns
s_rand <- adjustedRandIndex(labels, s)

# fg mixture
source("./SDS383c_Group4-main/CODE_original.R")
fgm <- FG_mixture(data, K = 10, M=25, Iter=500, alpha=1, mm=10)
index_fg = rep(1,nrow(data)) # cluster assignments 
# assignments are based on what sphere it is assigned 
L = ncol(fgm$inclusion_matrix)
for(j in 1:L) {
  index_fg[which(fgm$inclusion_matrix[,j]==1)]=j
}
fgm_rand <- adjustedRandIndex(labels, index_fg)

results <- data.frame(c("k-means", "dbscan", "speccalt", "FG"),
                      c(km_rand, db_rand, s_rand, fgm_rand))
colnames(results) <- c("Method", "RandIndex")
order <- c("dbscan", "k-means", "speccalt", "FG")

ggplot() + 
  geom_col(data=results, aes(x=Method, y=RandIndex)) +
  scale_x_discrete(limits = order) +
  ylab("adjusted Rand Index") +
  ggtitle("leaf dataset clustering")
ggsave("leaf.jpeg")
