library("tidyverse")
library("lme4")
library("MuMIn")
library("diptest")

rm(list=ls()) 
setwd("/Users/glmarschmann/.julia/dev/DEBmicroTrait/scripts/analysis")

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
BGE.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_BGE_new.csv") %>% 
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
BGE.f.grouped <- BGE.f %>% group_by(response, ontology) %>% summarise(BGE_med = median(BGE), rgrowth_med = median(rgrowth), BP_med = median(BP), BR_med = median(BR))
plot(BGE.f.grouped$BGE_med, BGE.f.grouped$rgrowth_med)
plot(BGE.f.grouped$BR_med, BGE.f.grouped$BP_med)


# Maximum Growth Rate
growthLOW <- lm(rgrowth_med ~ BGE_med, data = BGE.f.grouped[BGE.f.grouped$rgrowth_med < 0.0407, ])
summary(growthLOW)
growthLOW <- lm(BGE_med ~ rgrowth_med, data = BGE.f.grouped[BGE.f.grouped$rgrowth_med < 0.0407, ])
summary(growthLOW)
growthHIGH <- lm(rgrowth_med ~ BGE_med, data = BGE.f.grouped[BGE.f.grouped$rgrowth_med > 0.041, ])
summary(growthHIGH)
pred.frameLOW <- data.frame(BGE_med = seq(0.15, 0.6, by = 0.01))
pred.frameHIGH <- data.frame(BGE_med = seq(0.15, 0.6, by = 0.01))

# Respiration rate
respLOW <- lm(log(BR_med) ~ log(BP_med), data = BGE.f.grouped[BGE.f.grouped$rgrowth_med < 0.0407, ])
summary(respLOW)
respHIGH <- lm(log(BR_med) ~ log(BP_med), data = BGE.f.grouped[BGE.f.grouped$rgrowth_med > 0.041, ])
summary(respHIGH)

## Enzymes
## Growth rate - enzyme rate phase space
BGE.f.grouped <- BGE.f %>% group_by(response, ontology) %>% summarise(BGE_med = median(BGE), rgrowth_med = median(rgrowth), renzyme_med = median(renzyme))
plot(BGE.f.grouped$rgrowth_med, BGE.f.grouped$renzyme_med)

## Reserve density
BGE.f.grouped <- BGE.f %>% group_by(response, ontology) %>% summarise(BGE_med = median(BGE), rgrowth_med = median(rgrowth), dreserve_med = median(rdreserve))
plot(BGE.f.grouped$BGE_med, BGE.f.grouped$dreserve_med)

# Niche breadth
Levin.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_levin.csv") %>%
  replace(is.na(.), 0.0)
Levin.f.std <- Levin.f[,1:83]/(apply(Levin.f[,1:83], 1, sum))  

Levin.f.grouped <- BGE.f %>% group_by(isolate) %>% summarise(BGE_med = median(BGE), rgrowth_med = median(rgrowth))
Levin.f.grouped$levins <- levins(Levin.f.std)
Levin.f.grouped$response <- Levin.f$response
plot(Levin.f.grouped$levins, Levin.f.grouped$BGE_med)

BGELevin <- lm(BGE_med ~ levins, data = Levin.f.grouped)
summary(BGELevin)

BGELevin <- lm(rgrowth_med ~ levins, data = Levin.f.grouped)
summary(BGELevin)

pred.BGELevin <- data.frame(levins = seq(0.4, 1.0, by = 0.01))

dip.test(Levin.f.grouped$levins)
plot(density(Levin.f.grouped$levins))

kruskal.test(levins ~ response, Levin.f.grouped[Levin.f.grouped$response == "positive" | Levin.f.grouped$response == "undefined",])
write.csv(Levin.f.grouped, "/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_bge_growth_levin.csv")

porters.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files/isolates2traits.csv")
porters.f.std <- porters.f[,18:24]/(apply(porters.f[,18:24], 1, sum))  
Levin.f.grouped$levins <- levins(porters.f.std)
BGELevin <- lm(BGE_med ~ levins, data = Levin.f.grouped)
summary(BGELevin)
plot(Levin.f.grouped$levins, Levin.f.grouped$BGE_med)
pred.BGELevin <- data.frame(levins = seq(0.1, 1.0, by = 0.01))

# FCR - cell size
fcr.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_fcr.csv")
fcr.lm <- lm(log(FCR) ~ log(Vcell), data = fcr.f)
summary(fcr.lm)
plot(fcr.f$Vcell, fcr.f$FCR)

kruskal.test(FCR ~ response, fcr.f[fcr.f$response == "positive" | fcr.f$response == "undefined",])
kruskal.test(FCR ~ response, fcr.f[fcr.f$response == "positive" | fcr.f$response == "negative",])
kruskal.test(FCR ~ response, fcr.f[fcr.f$response == "undefined" | fcr.f$response == "negative",])





0.0?## Plot
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

# FCR - cell size
FCR.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_protein_synthesis.csv")
FCR.lm <- lm(log(FCR) ~ log(Vcell), data = FCR.f)
summary(FCR.lm)

kE.lm <- lm(log(k_E) ~ log(Vcell), data = FCR.f)
summary(kE.lm)

yEV.lm <- lm(log(y_EV) ~ log(Vcell), data = FCR.f)
summary(yEV.lm)


# Benchmark protein synthesis
kE.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/benchmark_protein_synthesis.csv")
kE.lm <- lm(diff ~ k_E_model, data = kE.f)
summary(kE.lm)
kE.lm <- lm(k_E_model ~ k_E, data = kE.f)
summary(kE.lm)

# CUE - genome size
genome.lm <- lm(log(BGE) ~ log(genomesize), data = BGE.f)
summary(genome.lm)


plot(BGE.f$BGE_med, BGE.f.grouped$dreserve_med)


ttest <- function(reg, coefnum, val){
  co <- coef(summary(reg))
  tstat <- (co[coefnum,1]-val)/co[coefnum,2]
  pstat <- 2 * pt(abs(tstat), reg$df.residual, lower.tail = FALSE)
  return(list = c(tstat, reg$df.residual, pstat))
}

BGE.f$Group <- NA
for (i in 1:dim(BGE.f)[1])
  if (BGE.f$rgrowth[i] < 0.0407){
    BGE.f$Group[i] <- "A"
  } else {
    BGE.f$Group[i] <- "B"
  }

mod1 <- lm(log(BP) ~ (log(BR)+ Group + ontology)^2, data = BGE.f)
summary(mod1)

modhigh <- lm(log(BP) ~ log(BR), data =  BGE.f[BGE.f$Group == "B",])
summary(modhigh)
ttest(modhigh,2,1)

modlow <- lm(log(BP) ~ log(BR), data =  BGE.f[BGE.f$Group == "A",])
summary(modlow)
ttest(modlow,2,1)


resp.f.grouped <- BGE.f %>% group_by(isolate) %>% summarise(BGE_med = median(BGE), rgrowth_med = median(rgrowth), BP_med = median(BP), BR_med = median(BR))
resp.f.grouped$Group <- NA
for (i in 1:dim(resp.f.grouped)[1])
  if (resp.f.grouped$rgrowth_med[i] < 0.0407){
    resp.f.grouped$Group[i] <- "A"
  } else {
    resp.f.grouped$Group[i] <- "B"
  }

modhigh <- lm(log(BP) ~ log(BR), data = resp.f.grouped[resp.f.grouped$Group == "B",])
summary(modhigh)


ttest(modhigh, 2, 1)


## Batch model BGE predictions
BGE.f.all <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_BGE_all.csv") %>% 
  filter(!is.nan(BGE))

BGE.f.all.high = BGE.f.all[BGE.f.all$rgrowth >= 0.055,] 
BGE.f.all.low = BGE.f.all[BGE.f.all$rgrowth < 0.0407, ]       

kruskal.test(BGE.f.all.high$Vmax[1:1134], BGE.f.all.low$Vmax[1:1134])
kruskal.test(BGE.f.all.high$KD[1:1134], BGE.f.all.low$KD[1:1134])
kruskal.test(BGE.f.all.high$rdreserve[1:1134], BGE.f.all.low$rdreserve[1:1134])
kruskal.test(BGE.f.all.high$rmaint[1:1134], BGE.f.all.low$rmaint[1:1134])


# thermodynamic efficiency
eta.f <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/exudates_eta.csv") 
uptake.f.grouped <- BGE.f.all %>% group_by(monomer) %>% summarise(up_med = median(ruptake), growth_med = median(rgrowth), BGE_med = median(BGE), yield_med = median(yield))
uptake.f.grouped$muED <- uptake.f.grouped$growth_med/uptake.f.grouped$yield_med
df_combined <- merge(eta.f, uptake.f.grouped, by="monomer")
plot(df_combined$muED, df_combined$eta, log='xy')
etareg <- lm(log(eta) ~ log(muED), data = df_combined)
summary(etareg)
etareg <- lm(log(eta) ~ log(growth_med), data = df_combined)
summary(etareg)
plot(df_combined$growth_med, df_combined$eta, log='xy')
plot(df_combined$muED, df_combined$lambda, log='xy')
etareg <- lm(log(lambda) ~ log(muED), data = df_combined)
summary(etareg)
plot(1/df_combined$yield_med, df_combined$lambda, log='xy')
etareg <- lm(log(lambda) ~ log(muED), data = df_combined)
summary(etareg)

ggplot(df_combined,aes(muED, eta)) +
  geom_point() +
  geom_smooth(method='lm') +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  labs(x = "Uptake rate [1/h]", y="Thermodynamic efficiency")
