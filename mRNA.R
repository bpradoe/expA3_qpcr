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
glimpse(table_mRNA)

# qqplot ---------------------------------------------------------------------
str(table_mRNA)

ExpA3_qq_time_target_mRNA <- table_mRNA |> 
  ggplot(aes(sample = ddct, fill = condition)) +
  stat_qq_band(alpha = 0.4) +
  stat_qq_point(shape = 21, size = 2) +
  stat_qq_line() +
  facet_wrap(~time*target, nrow = 2, scales = "free") +
  theme_minimal() +
  labs(title = "Q-Q Plot mRNA",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme(legend.position = "top")
ggsave("ExpB4_qqplot_time_mRNA.pdf", plot = ExpB4_qq_time_target_mRNA)

# boxplot -----------------------------------------------------------------

ExpB4_bp_time_target_mRNA <- table_mRNA |> 
  filter(ddct > 0.5 & ddct < 5) |> 
  ggplot(aes(x = group, y = ddct, fill = group)) +
  geom_jitter(width = 0.2, shape = 21, show.legend = FALSE) +
  geom_boxplot(alpha = 0.4, show.legend = FALSE) +
  facet_wrap(~time*target, scales = "free", nrow = 2) +
  labs(
    title = "Fold Change",
    fill = "Grupo"
  ) +
  xlab("") +
  ylab("FC") +
  theme_minimal() +
  theme()

ggsave("ExpB4_bp_time_target_mRNA.pdf", plot = ExpB4_bp_time_target_mRNA)


# estatísitca básica -----------------------------------------------------------
table_mRNA_out <- na.omit(table_mRNA) |>   # remover NAs e extremos
  group_by(time, target, condition) |> 
  mutate(out = is_outlier(ddct))
table_mRNA_out <- as.data.frame(table_mRNA_out) # tem que transformar em data.frame

shapiro_mRNA <- table_mRNA_out |> # normalidade
  group_by(condition, time) |> 
  shapiro_test(ddct) |> 
  filter(p < 0.05)

levene_mRNA <- table_mRNA_out |> # homogeneidade das variâncias
  levene_test(ddct ~ time)

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
# Barplot -------------------------------------------------

## LDHA

dunn_LDHA <- dunn_mRNA |> 
  filter(target == "LDHA")

barplot_LDHA <- table_mRNA_out |> 
  filter(target == "LDHA" & out == FALSE) |> 
  ggplot(aes(x = condition, y = ddct, fill = condition)) +
  stat_summary(fun.min = min, # coluna
               fun.max = max,
               fun.data = mean_se,
               geom = "col",
               width = 0.3,
               alpha = 1) +
  stat_summary(fun.min = min, # barra de erro
               fun.max = max,
               fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               alpha = 1)+
  geom_jitter(width = 0.1, shape = 21) +
  stat_pvalue_manual(dunn_LDHA, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 7, 
                     bracket.nudge.y = -0.01) +
  xlab("") +
  ylab("LDHA (FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#1A2CA3",   
      li  = "#F68048",   
      lps = "#2845D6"    
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "li" = "Li",
    "lps" = "LPS"
  )) +
  theme_light() +
  theme(legend.position = "none",
        strip.text = element_text(size = 7, margin = margin(1, 1, 1, 1))) +
  facet_wrap(~time, scales = "fixed", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)

## HIF1A
dunn_HIF1A <- dunn_mRNA |> 
  filter(target == "HIF1A")

barplot_HIF1A <- table_mRNA_out |> 
  filter(target == "HIF1A" & out == FALSE) |> 
  ggplot(aes(x = condition, y = ddct, fill = condition)) +
  stat_summary(fun.min = min, # coluna
               fun.max = max,
               fun.data = mean_se,
               geom = "col",
               width = 0.3,
               alpha = 1) +
  stat_summary(fun.min = min, # barra de erro
               fun.max = max,
               fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               alpha = 1)+
  geom_jitter(width = 0.1, shape = 21) +
  stat_pvalue_manual(dunn_HIF1A, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 7, 
                     bracket.nudge.y = -0.01) +
  xlab("") +
  ylab("HIF1A (FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#1A2CA3",   
      li  = "#F68048",   
      lps = "#2845D6"     
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "li" = "Li",
    "lps" = "LPS"
  )) +
  theme_light() +
  theme(legend.position = "none",
        strip.text = element_text(size = 7, margin = margin(1, 1, 1, 1))) +
  facet_wrap(~time, scales = "free", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)


## GLUT1
dunn_GLUT1 <- dunn_mRNA |> 
  filter(target == "GLUT1")

barplot_GLUT1 <- table_mRNA_out |> 
  filter(target == "GLUT1" & out == FALSE) |> 
  ggplot(aes(x = condition, y = ddct, fill = condition)) +
  stat_summary(fun.min = min, # coluna
               fun.max = max,
               fun.data = mean_se,
               geom = "col",
               width = 0.3,
               alpha = 1) +
  stat_summary(fun.min = min, # barra de erro
               fun.max = max,
               fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               alpha = 1)+
  geom_jitter(width = 0.1, shape = 21) +
  stat_pvalue_manual(dunn_GLUT1, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 7, 
                     bracket.nudge.y = -0.01) +
  xlab("") +
  ylab("GLUT1 (FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#1A2CA3",   
      li  = "#F68048",   
      lps = "#2845D6"     
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "li" = "Li",
    "lps" = "LPS"
  )) +
  theme_light() +
  theme(legend.position = "none",
        strip.text = element_text(size = 7, margin = margin(1, 1, 1, 1))) +
  facet_wrap(~time, scales = "fixed", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)

## PDK1
dunn_PDK1 <- dunn_mRNA |> 
  filter(target == "PDK1")

barplot_PDK1 <- table_mRNA_out |> 
  filter(target == "PDK1" & out == FALSE) |> 
  ggplot(aes(x = condition, y = ddct, fill = condition)) +
  stat_summary(fun.min = min, # coluna
               fun.max = max,
               fun.data = mean_se,
               geom = "col",
               width = 0.3,
               alpha = 1) +
  stat_summary(fun.min = min, # barra de erro
               fun.max = max,
               fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               alpha = 1)+
  geom_jitter(width = 0.1, shape = 21) +
  stat_pvalue_manual(dunn_PDK1, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 7, 
                     bracket.nudge.y = -0.01) +
  xlab("") +
  ylab("PDK1 (FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#1A2CA3",   
      li  = "#F68048",   
      lps = "#2845D6"     
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "li" = "Li",
    "lps" = "LPS"
  )) +
  theme_light() +
  theme(legend.position = "none",
        strip.text = element_text(size = 7, margin = margin(1, 1, 1, 1))) +
  facet_wrap(~time, scales = "fixed", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)