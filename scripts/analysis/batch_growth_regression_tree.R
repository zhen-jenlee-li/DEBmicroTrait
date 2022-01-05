library(tree)
library(MASS)
library(tidyverse)
library(stringr)
set.seed(1)

df_train <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_train_growth.csv")
df_train$gramstain <- as.numeric(df_train$gramstain)
df_train$yield <- 1- df_train$yield
df_train$affinity <- df_train$Vmax/df_train$KD

df_test <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_test_growth.csv")
df_test$gramstain <- as.numeric(df_test$gramstain)
df_test$yield <- 1- df_test$yield
df_test$affinity <- df_test$Vmax/df_test$KD

# log transform data
logdf_train <- (df_train) %>%
  filter(!is.nan(rgrowth) & affinity > 0.0)

logdf_test <- (df_test) %>%
  filter(!is.nan(rgrowth) & affinity > 0.0)

# simple tree
train = sample(1:nrow(logdf_train), nrow(logdf_train))
tree.bge = tree(rgrowth~., logdf_train, subset=train)
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
bge.test = logdf_test["rgrowth"]
plot(yhat,bge.test$rgrowth)
abline(0,1)
mean((yhat-bge.test$rgrowth)^2)
sqrt(mean((yhat-bge.test$rgrowth)^2))

library(randomForest)

# bagging m=p
bag.bge = randomForest(rgrowth~., data=logdf_train, subset=train, mtry=13, importance=TRUE)
yhat.bag = predict(bag.bge, newdata=logdf_test)
plot(yhat.bag, bge.test$rgrowth)
abline(0,1)
mean((yhat.bag-bge.test$rgrowth)^2)
sqrt(mean((yhat.bag-bge.test$rgrowth)^2))

# random forest m=sqrt(p)
rf.bge = randomForest(rgrowth~., data=logdf_train, subset=train, importance=TRUE)
yhat.rf = predict(rf.bge, newdata=logdf_test)
mean((yhat.rf-bge.test$rgrowth)^2)
sqrt(mean((yhat.rf-bge.test$rgrowth)^2))

# predictor importance
importance(rf.bge)
varImpPlot(rf.bge)

# boosting
library(gbm)

boost.bge = gbm(rgrowth~., data=logdf_train, distribution="gaussian", n.trees=500)
#svg("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/plots/batch_cue_regression_tree.svg")
summary(boost.bge)

df_all <- bind_rows(logdf_train, logdf_test)
df_low <- df_all %>% filter(rgrowth < 0.017)
df_high <- df_all %>% filter(rgrowth > 0.017)
lm.low <- lm(log(rgrowth)~log(Vmax), df_low)
summary(lm.low)
lm.low <- lm(log(rgrowth)~log(KD), df_low)
summary(lm.low)

lm.high <- lm(log(rgrowth)~log(Vmax/KD), df_high)
summary(lm.high)
lm.high <- lm(log(KD)~log(rgrowth), df_high)
summary(lm.high)
lm.low <- lm(log(KD)~log(rgrowth), df_low)
summary(lm.low)

write.csv(df_low, "/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_bge_low.csv")
write.csv(df_high, "/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_bge_high.csv")

plot(df_low$rgrowth, df_low$Vmax, log='xy')
plot(df_high$rgrowth, df_high$Vmax, log='xy')

plot(df_all$rgrowth, df_all$KD, log='xy')
plot(df_all$rgrowth, df_all$Vmax/df_all$KD, log='xy')

plot(df_low$rgrowth, df_low$KD, log='xy')
plot(df_low$rgrowth, df_low$Vmax, log='xy')
plot(df_low$rgrowth, df_low$affinity, log='xy')
plot(df_low$rgrowth, df_low$rho, log='xy')

plot(df_high$rgrowth, df_high$KD, log='xy')
plot(df_high$rgrowth, df_high$Vmax, log='xy')
plot(df_high$rgrowth, df_high$affinity, log='xy')
plot(df_high$rgrowth, df_high$rho, log='xy')

plot(df_all$rgrowth, df_all$KD, log='xy')
plot(df_all$rgrowth, df_all$Vmax, log='xy')
plot(df_all$rgrowth, df_all$affinity, log='xy')
plot(df_all$rgrowth, df_all$rho, log='xy')

lm.low <- lm(log(KD)~log(rgrowth), df_low)
summary(lm.low)
lm.high <- lm(log(KD)~log(rgrowth), df_high)
summary(lm.high)
