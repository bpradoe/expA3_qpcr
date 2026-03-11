library(tidyverse)
library(ggpubr)
library(qqplotr)
library(rstatix)
library(patchwork)
set.seed(2234)


# Fatores -----------------------------------------------------------------
table_mRNA <- read.csv("RNA_table - mRNA.csv", dec = ",")
table_mRNA$condition <- factor(table_mRNA$condition, levels = c("mo", "li", "lps"))
table_mRNA$time <- factor(table_mRNA$time, levels = c("1h", "2h", "4h", "24h"))
table_mRNA$target <- factor(table_mRNA$target, levels = c("HIF1A", "LDHA",
                                                          "GLUT1", "PDK1"))
# transformar o valor de fc de ddct para log2fc
table_mRNA$log2fc <- log2(table_mRNA$ddct)
head(table_mRNA)

# qqplot ---------------------------------------------------------------------
str(table_mRNA)

ExpA3_qq_time_target_mRNA <- table_mRNA |> 
  ggplot(aes(sample = log2fc, fill = condition)) +
  stat_qq_band(alpha = 0.4) +
  stat_qq_point(shape = 21, size = 2) +
  stat_qq_line() +
  facet_wrap(~time*target, nrow = 2, scales = "free") +
  theme_minimal() +
  labs(title = "Q-Q Plot mRNA",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme(legend.position = "top")

# estatísitca básica -----------------------------------------------------------
table_mRNA_out <- na.omit(table_mRNA) |>   # remover NAs e extremos
  group_by(time, target, condition) |> 
  mutate(out = is_outlier(log2fc))
table_mRNA_out <- as.data.frame(table_mRNA_out) # tem que transformar em data.frame
head(table_mRNA_out)

# ver quais grupos não são normais para
# aplicar teste não paramétrico
shapiro_mRNA <- table_mRNA_out |> 
  filter(out == FALSE) |> 
  group_by(condition, time, target) |> 
  shapiro_test(log2fc) |> 
  filter(p < 0.05)

# homogeneidade das variâncias
levene_mRNA <- table_mRNA_out |> 
  filter(out == FALSE) |> 
  levene_test(log2fc ~ time)

anova_mrna <- table_mRNA_out |>
  filter(out == FALSE) |> 
  group_by(target, time) |> 
  anova_test(log2fc ~ condition)
anova_mrna <- as.data.frame(anova_mrna) |> 
  filter(p > 0.05) 

tukey_mRNA <- table_mRNA_out |> 
  filter(out == FALSE) |> 
  group_by(time, target) |> 
  tukey_hsd(log2fc ~ condition)|> 
  add_significance() |> 
  add_xy_position(x = "condition", step.increase = 0.1)

# daqui pra frente eu fiz os barplots e os teste por alvo


# LDHA --------------------------------------------------------------------

# filtrei pelo alvo, os outliers e o tempo
# que foi significativo no ANOVA
tukey_mrna_ldha <- table_mRNA_out |> 
  filter(target == "LDHA" & out == FALSE) |> 
  group_by(time) |> 
  tukey_hsd(log2fc ~ condition)|> 
  add_significance() |> 
  add_xy_position(x = "condition", step.increase = 0.1)

# checagem dos não-significativos
tukey_mrna_ldha |> 
  filter(p.adj > 0.05)

table_mRNA_out |> 
  filter(target == "LDHA" & out == FALSE) |> 
  ggplot(aes(x = condition, y = log2fc, fill = condition)) +
  stat_summary(fun.min = min, # coluna
               fun.max = max,
               fun.data = mean_se,
               geom = "col",
               width = 0.5,
               alpha = 1) +
  stat_summary(fun.min = min, # barra de erro
               fun.max = max,
               fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               alpha = 1)+
  geom_jitter(width = 0.1, shape = 21) +
  xlab("") +
  ylab("LDHA (log2FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#BFC6C4",   
      li  = "#8A2D3B",   
      lps = "#E8E2D8"        
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "li" = "Li",
    "lps" = "LPS"
  )) +
  stat_pvalue_manual(tukey_mrna_ldha, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 7) +
  theme_light() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16)) +
  facet_wrap(~time, scales = "fixed", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)


# HIF1A -------------------------------------------------------------------

# filtrei pelo alvo, os outliers e o tempo
# que foi significativo no ANOVA
tukey_mir_hif1a <- table_mRNA_out |> 
  filter(target == "HIF1A" & out == FALSE) |> 
  group_by(time) |> 
  tukey_hsd(log2fc ~ condition)|> 
  add_significance() |> 
  add_xy_position(x = "condition", step.increase = 0.1)

# checagem dos não-significativos
tukey_mir_hif1a |> 
  filter(p.adj > 0.05)

table_mRNA_out |> 
  filter(target == "HIF1A" & out == FALSE) |> 
  ggplot(aes(x = condition, y = log2fc, fill = condition)) +
  stat_summary(fun.min = min, # coluna
               fun.max = max,
               fun.data = mean_se,
               geom = "col",
               width = 0.5,
               alpha = 1) +
  stat_summary(fun.min = min, # barra de erro
               fun.max = max,
               fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               alpha = 1)+
  geom_jitter(width = 0.1, shape = 21) +
  xlab("") +
  ylab("HIF1A (log2FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#BFC6C4",   
      li  = "#8A2D3B",   
      lps = "#E8E2D8"        
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "li" = "Li",
    "lps" = "LPS"
  )) +
  stat_pvalue_manual(tukey_mir_hif1a, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 7) +
  theme_light() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16)) +
  facet_wrap(~time, scales = "fixed", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)


# GLUT1/SLC2A1 ------------------------------------------------------------

tukey_mrna_glut1 <- table_mRNA_out |> 
  filter(target == "GLUT1" & out == FALSE) |> 
  group_by(time) |> 
  tukey_hsd(log2fc ~ condition)|> 
  add_significance() |> 
  add_xy_position(x = "condition", step.increase = 0.1)

# checagem dos não-significativos
tukey_mrna_glut1 |> 
  filter(p.adj > 0.05)

table_mRNA_out |> 
  filter(target == "GLUT1" & out == FALSE) |> 
  ggplot(aes(x = condition, y = log2fc, fill = condition)) +
  stat_summary(fun.min = min, # coluna
               fun.max = max,
               fun.data = mean_se,
               geom = "col",
               width = 0.5,
               alpha = 1) +
  stat_summary(fun.min = min, # barra de erro
               fun.max = max,
               fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               alpha = 1)+
  geom_jitter(width = 0.1, shape = 21) +
  xlab("") +
  ylab("GLUT1 (log2FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#BFC6C4",   
      li  = "#8A2D3B",   
      lps = "#E8E2D8"        
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "li" = "Li",
    "lps" = "LPS"
  )) +
  stat_pvalue_manual(tukey_mrna_glut1, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 7) +
  theme_light() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16)) +
  facet_wrap(~time, scales = "fixed", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)


# PDK1 --------------------------------------------------------------------


tukey_mrna_pdk1 <- tukey_mRNA |> 
  filter(target == "PDK1")


table_mRNA_out |> 
  filter(target == "PDK1" & out == FALSE) |> 
  ggplot(aes(x = condition, y = log2fc, fill = condition)) +
  stat_summary(fun.min = min, # coluna
               fun.max = max,
               fun.data = mean_se,
               geom = "col",
               width = 0.5,
               alpha = 1) +
  stat_summary(fun.min = min, # barra de erro
               fun.max = max,
               fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               alpha = 1)+
  geom_jitter(width = 0.1, shape = 21) +
  xlab("") +
  ylab("PDK1 (log2FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#BFC6C4",   
      li  = "#8A2D3B",   
      lps = "#E8E2D8"        
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "li" = "Li",
    "lps" = "LPS"
  )) +
  stat_pvalue_manual(tukey_mrna_pdk1, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 7) +
  theme_light() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16)) +
  facet_wrap(~time, scales = "fixed", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)



# outras coisas -----------------------------------------------------------


kruskal_mRNA <- table_mRNA_out |> 
  group_by(time, target) |> 
  kruskal_test(ddct ~ condition) |> 
  filter(p < 0.05) # filtrar os que foram signifcativos

dunn_mRNA <- table_mRNA_out |> 
  filter(target == kruskal_mRNA$target & time == kruskal_mRNA$time) |> 
  group_by(time, target) |>                            
  dunn_test(ddct ~ condition,    
            p.adjust.method = "bonferroni") |> 
  add_significance() |> 
  add_xy_position(x = "condition", step.increase = 0.01)

dunn_LDHA <- dunn_mRNA |> 
  filter(target == "LDHA")



tukey_mrna_glut1 <- tukey_mRNA |> 
  filter(target == "GLUT1")
