library(tidyverse)
library(ggpubr)
library(ggbreak)
library(qqplotr)
library(rstatix)
library(patchwork)
set.seed(2234)

# fatores -----------------------------------------------------------------
table_microRNA <- read.csv("RNA_table - microRNA.csv", dec = ",")
table_microRNA$condition <- factor(table_microRNA$condition, levels = c("mo", "Li", "lps"))
table_microRNA$time <- factor(table_microRNA$time, levels = c("1h", "2h", "4h", "24h"))
table_microRNA$target <- factor(table_microRNA$target, levels = c("mir_372", "mir_373"))

# transformar o valor de fc de ddct para log2fc
table_microRNA$log2fc <- log(table_microRNA$ddct, base = exp(2))
head(table_microRNA)

# qqplot ---------------------------------------------------------------------
str(table_microRNA)

ExpA3_qq_time_target_mRNA <- table_microRNA |>
  filter(condition != "Li") |> # essa condição deixa difícil de interpretar
  ggplot(aes(sample = ddct, fill = condition)) +
  stat_qq_band(alpha = 0.4) +
  stat_qq_point(shape = 21, size = 2) +
  stat_qq_line() +
  facet_wrap(~time*target, nrow = 2, scales = "free") +
  theme_minimal() +
  labs(title = "Q-Q Plot miR-372",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme(legend.position = "top")
ggsave("ExpA3_qqplot_time_microRNA.pdf", plot = ExpA3_qq_time_target_microRNA)

# estatísitca básica -----------------------------------------------------------
table_microRNA_out <- na.omit(table_microRNA) |>   # remover NAs e extremos
  group_by(time, target, condition) |> 
  mutate(out = is_outlier(log2fc))
table_microRNA_out <- as.data.frame(table_microRNA_out) 
head(table_microRNA_out)

# verificar quais não são normais para realizar teste não-paramétrico
shapiro_microRNA <- table_microRNA_out |>
  filter(out == FALSE) |> 
  group_by(condition, time, target) |> 
  shapiro_test(log2fc) |> 
  filter(p < 0.05)

# homogeneidade das variâncias
levene_microRNA <- table_microRNA_out |>
  filter(out == FALSE) |> 
  levene_test(log2fc ~ time)

# verificar quais não foram significativas para 
# filtrar no teste de Tukey
anova_mir <- table_microRNA_out |>
  filter(out == FALSE) |> 
  group_by(target, time) |> 
  anova_test(log2fc ~ condition)
anova_mir <- as.data.frame(anova_mir) |> 
  filter(p > 0.05) 


# ajustar o p-valor para comparações múltiplas
tukey_mir_372 <- table_microRNA_out |> 
  filter(target == "mir_372" & out == FALSE) |> 
  group_by(time) |> 
  tukey_hsd(log2fc ~ condition)|> 
  add_significance() |> 
  add_xy_position(x = "condition", step.increase = 0.1)

# ver quantas comparações não foram significativas
tukey_mir_372 |> 
  filter(p.adj > 0.05) |> 
  nrow()

tukey_mir_373 <- table_microRNA_out |> 
  filter(target == "mir_373" & out == FALSE) |> 
  group_by(time) |> 
  tukey_hsd(log2fc ~ condition)|> 
  add_significance() |> 
  add_xy_position(x = "condition", step.increase = 0.1)

# ver quantas comparações não foram significativas
tukey_mir_373 |> 
  filter(p.adj > 0.05) |> 
  nrow()
# barplot -------------------------------------------------


## miR-372
table_microRNA_out |> 
  filter(target == "mir_372" & out == FALSE) |> 
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
               alpha = 1) +
  geom_jitter(width = 0.1, shape = 21) +
  xlab("") +
  ylab("miR-372 (log2FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#BFC6C4",   
      Li  = "#8A2D3B",   
      lps = "#E8E2D8"        
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "Li" = "Li",
    "lps" = "LPS"
  )) +
  stat_pvalue_manual(tukey_mir_372, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 4, 
                     bracket.nudge.y = -0.01, 
                     step.group.by = "time") +
  theme_light() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16)) +
  facet_wrap(~time, scales = "fixed", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)



## miR-373
table_microRNA_out |> 
  filter(target == "mir_373" & out == FALSE) |> 
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
               alpha = 1) +
  geom_jitter(width = 0.1, shape = 21) +
  xlab("") +
  ylab("miR-373 (log2FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#BFC6C4",   
      Li  = "#8A2D3B",   
      lps = "#E8E2D8"        
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "Li" = "Li",
    "lps" = "LPS"
  )) +
  stat_pvalue_manual(tukey_mir_373, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 4, 
                     bracket.nudge.y = -0.01, 
                     step.group.by = "time") +
  theme_light() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 16))+
  facet_wrap(~time, scales = "fixed", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)

