library(tree)
library(MASS)
library(tidyverse)
library(stringr)
set.seed(1)

df_train <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_train.csv")
df_train$gramstain <- as.numeric(df_train$gramstain)
df_train$yield <- 1- df_train$yield

df_test <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_test.csv")
df_test$gramstain <- as.numeric(df_test$gramstain)
df_test$yield <- 1- df_test$yield

# log transform data
logdf_train <- log(df_train) %>%
   filter(!is.nan(BGE))

logdf_test <- log(df_test) %>%
  filter(!is.nan(BGE))

# simple tree
train = sample(1:nrow(logdf_train), nrow(logdf_train))
tree.bge = tree(BGE~., logdf_train, subset=train)
summary(tree.bge)
plot(tree.bge)
text(tree.bge)

cv.bge = cv.tree(tree.bge)
plot(cv.bge$size, cv.bge$dev, type='b')

# pruning
prune.bge = prune.tree(tree.bge, best=5)
plot(prune.bge)
text(prune.bge, pretty=0)

yhat = predict(tree.bge, newdata=logdf_test)
bge.test = logdf_test["BGE"]
plot(yhat,bge.test$BGE)
abline(0,1)
mean((yhat-bge.test$BGE)^2)
sqrt(mean((yhat-bge.test$BGE)^2))

library(randomForest)

# bagging m=p
bag.bge = randomForest(BGE~., data=logdf_train, subset=train, mtry=13, importance=TRUE)
yhat.bag = predict(bag.bge, newdata=logdf_test)
plot(yhat.bag, bge.test$BGE)
abline(0,1)
mean((yhat.bag-bge.test$BGE)^2)
sqrt(mean((yhat.bag-bge.test$BGE)^2))

# random forest m=sqrt(p)
rf.bge = randomForest(BGE~., data=logdf_train, subset=train, importance=TRUE)
yhat.rf = predict(rf.bge, newdata=logdf_test)
mean((yhat.rf-bge.test$BGE)^2)
sqrt(mean((yhat.rf-bge.test$BGE)^2))

# predictor importance
importance(rf.bge)
varImpPlot(rf.bge)

# boosting
library(gbm)

boost.bge = gbm(BGE~., data=logdf_train, distribution="gaussian", n.trees=500)
#svg("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/plots/batch_cue_regression_tree.svg")
summary(boost.bge)

# partial dependence
par(mfrow=c(1,2))
plot(boost.bge, i="genomesize")
plot(boost.bge, i="gc")

yhat.boost = predict(boost.bge, newdata=logdf_test, n.tress=500)
svg("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/plots/batch_cue_regression_error.svg")
plot(yhat.boost, bge.test$BGE)
abline(0,1)
dev.off()
mean((yhat.boost-bge.test$BGE)^2)

print(pretty.gbm.tree(boost.bge, i.tree = boost.bge$n.trees))

# shrinkage parameter
boost.boston = gbm(medv~., data=Boston[train,], distribution="gaussian", n.trees=500, interaction.depth=4, shrinkage=0.2, verbose=F)
yhat.boost = predict(boost.boston, newdata=Boston[-train,], n.tress=500)
mean((yhat.boost-boston.test)^2)

df_all <- bind_rows(df_train, df_test)
df_low <- df_all %>% filter(BGE < 0.25)
df_high <- df_all %>% filter(BGE > 0.25)
lm.low <- lm(log(rho)~log(yield), df_low)
summary(lm.low)
lm.low.poly <- lm(log(rho)~poly(log(yield),3), df_low)
summary(lm.low.poly)

lm.high <- lm(log(rho)~log(yield), df_high)
summary(lm.high)
lm.high.poly <- lm(log(rho)~poly(log(yield),3), df_high)
summary(lm.high.poly)

write.csv(df_low, "/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_bge_low.csv")
write.csv(df_high, "/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_bge_high.csv")

