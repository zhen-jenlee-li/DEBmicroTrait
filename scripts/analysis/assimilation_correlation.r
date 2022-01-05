library("PerformanceAnalytics")
library("tidyverse")

df <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_assimilation_sav1.csv")

df_norm <- df %>% 
  filter(transporter_density >1e-8) %>%
  mutate(transporter_density_log = log(transporter_density),
         sav_size_log = log(sav),
         Vmax_log = log(Vmax),
         KD_log = log(KD),
         affinity_log = log(affinity),
         #affinity_log = log(KD/Vmax)
         ) 

df_numeric <- df_norm[,12:16]
#svg("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/plots/test.svg")
chart.Correlation(df_numeric, histogram=TRUE, label="", method="spearman")
#dev.off()


