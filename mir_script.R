library(tidyverse)
library(ggpubr)
library(ggbreak)
library(qqplotr)
library(rstatix)
library(patchwork)
set.seed(2234)

# fatores -----------------------------------------------------------------
table_microRNA <- read.csv("miR_table - principal.csv", dec = ",")
table_microRNA$condition <- factor(table_microRNA$condition, levels = c("mo", "Li", "lps"))
table_microRNA$time <- factor(table_microRNA$time, levels = c("1h", "2h", "4h", "24h"))
table_microRNA$target <- factor(table_microRNA$target, levels = c("mir_372", "mir_373"))


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
table_microRNA_out <- na.omit(table_microRNA) |>   # remover NAs e extremos
  group_by(time, target, condition) |> 
  mutate(out = is_outlier(ddct))
table_microRNA_out <- as.data.frame(table_microRNA_out) 


shapiro_microRNA <- table_microRNA_out |> # normalidade
  group_by(condition, time, target) |> 
  shapiro_test(ddct) |> 
  filter(p < 0.05)

levene_microRNA <- table_microRNA_out |> # homogeneidade das variâncias
  levene_test(ddct ~ time)

anova_mir <- table_microRNA_out |> 
  group_by(target, time) |> 
  anova_test(ddct ~ condition)
anova_mir <- as.data.frame(anova_mir) |> 
  filter(p > 0.05)

tukey_mir_372 <- table_microRNA_out |> 
  filter(target == "mir_372" & time != "24h") |> 
  group_by(time) |> 
  tukey_hsd(ddct ~ condition)|> 
  add_significance() |> 
  add_xy_position(x = "condition", step.increase = 0.1)

tukey_mir_373 <- table_microRNA_out |> 
  filter(target == "mir_373") |> 
  group_by(time) |> 
  tukey_hsd(ddct ~ condition)|> 
  add_significance() |> 
  add_xy_position(x = "condition", step.increase = 0.1)

# barplot -------------------------------------------------


## miR-372
barplot_mir_372 <- table_microRNA_out |> 
  filter(target == "mir_372") |> 
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
               alpha = 1) +
  geom_jitter(width = 0.1, shape = 21) +
  stat_pvalue_manual(tukey_mir_372, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 4, 
                     bracket.nudge.y = -0.01, 
                     step.group.by = "time") +
  xlab("") +
  ylab("miR-372 (FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#1A2CA3",   
      Li  = "#F68048",   
      lps = "#2845D6"        
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "Li" = "Li",
    "lps" = "LPS"
  )) +
  theme_light() +
  theme(legend.position = "none")+
  facet_wrap(~time, scales = "fixed", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)



## miR-373
barplot_mir_373 <- table_microRNA |> 
  filter(target == "mir_373") |> 
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
               alpha = 1) +
  geom_jitter(width = 0.1, shape = 21) +
  stat_pvalue_manual(tukey_mir_373, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE, label.size = 4, 
                     bracket.nudge.y = -0.01, 
                     step.group.by = "time") +
  xlab("") +
  ylab("miR-373 (FC)") +
  scale_fill_manual(
    values = c(
      mo  = "#1A2CA3",   
      Li  = "#F68048",   
      lps = "#2845D6"        
    )) +
  scale_x_discrete(labels = c(
    "mo" = "Mo",
    "Li" = "Li",
    "lps" = "LPS"
  )) +
  theme_light() +
  theme(legend.position = "none")+
  facet_wrap(~time, scales = "fixed", shrink = FALSE, 
             strip.position = "bottom", nrow = 1)
# patchwork ---------------------------------------------------------------

barplots_ExpA3 <- (barplot_mir_372/barplot_mir_373)
ggsave("ExpA3_barplots_mir.pdf", plot = get_last_plot(), scale = 3)


