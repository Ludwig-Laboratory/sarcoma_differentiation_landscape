# TARGET OS prediatric sarcoma data
# survival analysis with MTL archetypes

library(tidyverse)


# load archetypes (same column names as expression data)
target_expr = read_delim("Data/TARGET_OS/target_os.txt")
arch = read_csv("Data/TARGET_OS/TG_Arch.csv", col_names = F) %>% as.matrix
colnames(arch) = colnames(target_expr)
rownames(arch) = paste("arch",1:12,sep="_")

# heatmap of archetype scores
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.25, 0.5), c("#009FFF", "#FFFF00", "#FF5100"))

Heatmap(arch, name = "NMF Score", col = col_fun,
        cluster_rows = F, row_names_side = "left",
        row_title = "Archetype", row_title_side = "left",
        column_title = "Osteosarcoma Sample", show_column_names = F, show_column_dend = F)


# load survival data
surv = read_tsv("Data/TARGET_OS/from_Xena/TARGET-OS.survival.tsv") %>% as.data.frame

# match the sample indices
match_ids = match(colnames(arch), surv$sample)

# combine survival dataframe with archetype scores
surv_arch = cbind(surv[match_ids,] %>% dplyr::select(OS, OS.time),
                  t(arch) ) %>% filter(!is.na(OS))


# assess archetype significance with Cox proportional hazard test
library(survival)
library(survminer)

# univariable Cox PH regressions on each archetype
coxph(Surv(OS.time, OS) ~ arch_1, data=surv_arch)
coxph(Surv(OS.time, OS) ~ arch_2, data=surv_arch)
coxph(Surv(OS.time, OS) ~ arch_3, data=surv_arch) #***
coxph(Surv(OS.time, OS) ~ arch_4, data=surv_arch)
coxph(Surv(OS.time, OS) ~ arch_5, data=surv_arch) #*
coxph(Surv(OS.time, OS) ~ arch_6, data=surv_arch) #*
coxph(Surv(OS.time, OS) ~ arch_7, data=surv_arch)
coxph(Surv(OS.time, OS) ~ arch_8, data=surv_arch)
coxph(Surv(OS.time, OS) ~ arch_9, data=surv_arch)
coxph(Surv(OS.time, OS) ~ arch_10, data=surv_arch)
coxph(Surv(OS.time, OS) ~ arch_11, data=surv_arch)
coxph(Surv(OS.time, OS) ~ arch_12, data=surv_arch)

# filter out low expressed archetypes for multivariable Cox PH
keep_archs = which(rowMeans(arch) > 0.05)
fit_filt = coxph(Surv(OS.time, OS) ~ ., data=surv_arch[,c(1,2,keep_archs+2)])
summary(fit_filt)
ggforest(fit_filt, cpositions = c(0.02, 0.12, 0.32))


# plot Kaplan-Meier curves for archetype 3 high/low
fit3 = survfit(Surv(OS.time, OS) ~ arch3_grp, data=surv_arch %>%
                 mutate(arch3_grp = if_else(arch_3>median(arch_3),"high","low")) %>%
                 mutate(OS.time = OS.time/365.25*12))
ggsurvplot(fit3, pval=T)

