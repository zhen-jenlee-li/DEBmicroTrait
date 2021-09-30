library("tidyverse")
library("lme4")
library("MuMIn")
library("diptest")

rm(list=ls()) 
setwd("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/scripts/analysis")

opar <- par(no.readonly = TRUE)  # Saves plot defaults

## Auxiliary functions
# Confidence Hulls
add.hull <- function(model = "", pred.frame = ""){
  CI.U <- predict(model, interval = "c", newdata=pred.frame)[, "upr"]
  CI.L <- predict(model, interval = "c", newdata=pred.frame)[, "lwr"]
  pred.frame2 <- unlist(pred.frame)
  X.Vec <- c(pred.frame2, tail(pred.frame2, 1), rev(pred.frame2),
             head(pred.frame2, 1))
  Y.Vec <- c(CI.U, tail(CI.L, 1), rev(CI.L), head(CI.U,1))
  polygon(X.Vec, Y.Vec, col = "gray90", border = NA)
}
# Levin's index
levins <- function(p_xi = ""){
  p = 0
  for (i in p_xi){
    p = p + i^2
  }
  nb = 1 / (length(p_xi) * p)
  return(nb)
}

## Batch model BGE predictions
BGE.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files_pub/isolates_batch_model_BGE.csv") %>% 
         filter(!is.nan(BGE))

## Taxonomic and Resource Variance Partitioning
# Model 1: Isolate Identity 
fm1 <- lmer(BGE ~ as.factor(isolate) + (1|monomer/ontology), data = BGE.f)
summary(fm1)$AIC
r.squaredGLMM(fm1)
qqnorm(resid(fm1))

# Model 2: Class Order
fm2 <- lmer(BGE ~ as.factor(class) + (1|monomer/ontology), data = BGE.f)
summary(fm2)$AIC
r.squaredGLMM(fm2)
qqnorm(resid(fm2))

# Model 3: Phylum Order
fm3 <- lmer(BGE ~ as.factor(phylum) + (1|monomer), data = BGE.f)
summary(fm3)$AIC
r.squaredGLMM(fm3)
qqnorm(resid(fm3))

# Model 4: Monomer Identity
fm4 <- lmer(BGE ~ as.factor(monomer) + (1|isolate/class/phylum), data = BGE.f)
summary(fm4)$AIC
r.squaredGLMM(fm4)
plot(predict(fm4), resid(fm4))
qqnorm(resid(fm4))

# Model 5: Monomer Class
fm5 <- lmer(BGE ~ as.factor(ontology) + (1|isolate/class/phylum), data = BGE.f)
summary(fm5)$AIC
r.squaredGLMM(fm5)
plot(predict(fm5), resid(fm5))
qqnorm(resid(fm5))

## Growth rate - CUE phase space
BGE.f.grouped <- BGE.f %>% group_by(response, ontology) %>% summarise(BGE_med = median(BGE), rgrowth_med = median(rgrowth))
plot(BGE.f.grouped$BGE_med, BGE.f.grouped$rgrowth_med)

# Maximum Growth Rate
growthLOW <- lm(rgrowth_med ~ BGE_med, data = BGE.f.grouped[BGE.f.grouped$rgrowth_med < 0.0407, ])
summary(growthLOW)
growthHIGH <- lm(rgrowth_med ~ BGE_med, data = BGE.f.grouped[BGE.f.grouped$rgrowth_med > 0.041, ])
summary(growthHIGH)
pred.frameLOW <- data.frame(BGE_med = seq(0.15, 0.6, by = 0.01))
pred.frameHIGH <- data.frame(BGE_med = seq(0.15, 0.6, by = 0.01))

## Enzymes
## Growth rate - enzyme rate phase space
BGE.f.grouped <- BGE.f %>% group_by(response, ontology) %>% summarise(BGE_med = median(BGE), rgrowth_med = median(rgrowth), renzyme_med = median(renzyme))
plot(BGE.f.grouped$rgrowth_med, BGE.f.grouped$renzyme_med)

## Reserve density
BGE.f.grouped <- BGE.f %>% group_by(response, ontology) %>% summarise(BGE_med = median(BGE), rgrowth_med = median(rgrowth), dreserve_med = median(rdreserve))
plot(BGE.f.grouped$BGE_med, BGE.f.grouped$dreserve_med)

# Niche breadth
Levin.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files_pub/isolates_batch_model_levin.csv") %>%
  replace(is.na(.), 0.0)
Levin.f.std <- Levin.f[,1:83]/(apply(Levin.f[,1:83], 1, sum))  

Levin.f.grouped <- BGE.f %>% group_by(isolate) %>% summarise(BGE_med = median(BGE))
Levin.f.grouped$levins <- levins(Levin.f.std)
Levin.f.grouped$response <- Levin.f$response
plot(Levin.f.grouped$levins, Levin.f.grouped$BGE_med)

BGELevin <- lm(BGE_med ~ levins, data = Levin.f.grouped)
summary(BGELevin)
pred.BGELevin <- data.frame(levins = seq(0.4, 1.0, by = 0.01))

dip.test(Levin.f.grouped$levins)
plot(density(Levin.f.grouped$levins))

kruskal.test(levins ~ response, Levin.f.grouped[Levin.f.grouped$response == "positive" | Levin.f.grouped$response == "undefined",])


porters.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files/isolates2traits.csv")
porters.f.std <- porters.f[,18:24]/(apply(porters.f[,18:24], 1, sum))  
Levin.f.grouped$levins <- levins(porters.f.std)
BGELevin <- lm(BGE_med ~ levins, data = Levin.f.grouped)
summary(BGELevin)
plot(Levin.f.grouped$levins, Levin.f.grouped$BGE_med)
pred.BGELevin <- data.frame(levins = seq(0.1, 1.0, by = 0.01))



?## Plot
setEPS()
postscript("../../plots/IsolatesBatchBGE.eps", width=4.5, height=4.5, bg = "white")
par(opar)

# BGE-rate
plot(rgrowth_med ~ BGE_med, BGE.f.grouped,  axes = F, type = "n",
     xlab = "", ylab = "",
     xlim = c(0.14, 0.6), ylim = c(0.005,0.08), las = 1,
     pch = 22, bg = "gray", lwd = 2, cex = 1.5)
add.hull(growthLOW, pred.frameLOW)
matlines(pred.frameLOW, predict(growthLOW, interval = "c", newdata=pred.frameLOW),
         lty=c(2,3,3), lwd=c(4,2,2), col="black")
add.hull(growthHIGH, pred.frameHIGH)
matlines(pred.frameHIGH, predict(growthHIGH, interval = "c", newdata=pred.frameHIGH),
         lty=c(2,3,3), lwd=c(4,2,2), col="black")
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "positive" & BGE.f.grouped$ontology == "Amino acids",],
       pch = 24, bg = "#1b9e77", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "positive" & BGE.f.grouped$ontology == "Organic acids",],
       pch = 24, bg = "#d95f02", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "positive" & BGE.f.grouped$ontology == "Nucleotides",],
       pch = 24, bg = "#7570b3", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "positive" & BGE.f.grouped$ontology == "Sugars",],
       pch = 24, bg = "#e7298a", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "positive" & BGE.f.grouped$ontology == "Auxins",],
       pch = 24, bg = "#66a61e", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "positive" & BGE.f.grouped$ontology == "Fatty acids",],
       pch = 24, bg = "#e6ab02", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "negative" & BGE.f.grouped$ontology == "Amino acids",],
       pch = 25, bg = "#1b9e77", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "negative" & BGE.f.grouped$ontology == "Organic acids",],
       pch = 25, bg = "#d95f02", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "negative" & BGE.f.grouped$ontology == "Nucleotides",],
       pch = 25, bg = "#7570b3", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "negative" & BGE.f.grouped$ontology == "Sugars",],
       pch = 25, bg = "#e7298a", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "negative" & BGE.f.grouped$ontology == "Auxins",],
       pch = 25, bg = "#66a61e", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "negative" & BGE.f.grouped$ontology == "Fatty acids",],
       pch = 25, bg = "#e6ab02", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "undefined" & BGE.f.grouped$ontology == "Amino acids",],
       pch = 22, bg = "#1b9e77", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "undefined" & BGE.f.grouped$ontology == "Organic acids",],
       pch = 22, bg = "#d95f02", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "undefined" & BGE.f.grouped$ontology == "Nucleotides",],
       pch = 22, bg = "#7570b3", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "undefined" & BGE.f.grouped$ontology == "Sugars",],
       pch = 22, bg = "#e7298a", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "undefined" & BGE.f.grouped$ontology == "Auxins",],
       pch = 22, bg = "#66a61e", lwd = 2, cex = 1.5)
points(rgrowth_med ~ BGE_med, BGE.f.grouped[BGE.f.grouped$response == "undefined" & BGE.f.grouped$ontology == "Fatty acids",],
       pch = 22, bg = "#e6ab02", lwd = 2, cex = 1.5)
axis(1, lwd = 2, labels = T, at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), las = 1)
axis(2, lwd = 2, labels = T, at = c(0, 0.005, 0.02, 0.04, 0.06, 0.08, 0.1), las = 1)
mtext(expression(paste("Maximum Growth Rate (h"^-1,")")), side = 2, line = 3, cex = 1)
mtext("Carbon Use Efficiency [-]", side = 1, line = 2.5, cex = 1)
legend("center", c("Positive", "Negative", "Undefined"), pch = c(24, 25, 22), pt.bg = "white", pt.cex = 1, pt.lwd = 2, bty = "n", cex = 0.8)
legend("left", c("Amino acids", "Organic acids", "Nucleotides", "Sugars", "Auxins", "Fatty acids"), pch = c(22), pt.bg = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02" ), pt.cex = 1, pt.lwd = 2, bty = "n", cex = 0.8)
box(lwd = 2)
dev.off() # this writes plot to folder
graphics.off() # shuts down open devices

## Plot
setEPS()
postscript("../../plots/IsolatesBatchLevin.eps", width=4.5, height=4.5, bg = "white")
par(opar)
# BGE-Levin
plot(levins ~ BGE_med, Levin.f.grouped, axes = F, type = "n",
     xlab = "", ylab = "",
     xlim = c(0.0, 1.0), ylim = c(0.2, 0.8), las = 1,
     pch = 22, bg = "gray", lwd = 2, cex = 1.5)
add.hull(BGELevin, pred.BGELevin)
matlines(pred.BGELevin, predict(BGELevin, interval = "c", newdata=pred.BGELevin),
         lty=c(2,3,3), lwd=c(4,2,2), col="black")
points(BGE_med ~ levins, Levin.f.grouped[Levin.f.grouped$response == "positive",],
       pch = 24, bg = "white", lwd = 2, cex = 1.5)
points(BGE_med ~ levins, Levin.f.grouped[Levin.f.grouped$response == "negative",],
       pch = 25, bg = "white", lwd = 2, cex = 1.5)
points(BGE_med ~ levins, Levin.f.grouped[Levin.f.grouped$response == "undefined",],
       pch = 22, bg = "white", lwd = 2, cex = 1.5)
axis(1, lwd = 2, labels = T, at = c(0.4, 0.6, 0.8, 1.0), las = 1)
axis(2, lwd = 2, labels = T, at = c(0.2, 0.4, 0.6, 0.8), las = 1)
mtext(expression(paste("Levins index [-]")), side = 1, line = 2.5, cex = 1)
mtext("Carbon Use Efficiency [-]", side = 2, line = 2.5, cex = 1)
legend("bottomright", c("Positive", "Negative", "Undefined"), pch = c(24, 25, 22), pt.bg = c("white") , pt.cex = 1, pt.lwd = 2, bty = "n", cex = 0.8)
box(lwd = 2)
dev.off() # this writes plot to folder
graphics.off() # shuts down open devices



## Batch model BGE predictions





BGE.f.grouped$levin = Levin.f$levins
growthLOW <- lm(BGE_med ~ levin, data = BGE.f.grouped)
summary(growthLOW)

   