library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(broom)
library(effectsize)
library(knitr)
library(kableExtra)
library(ggplot2)

res <- readRDS("gslid_sim_results_long.rds")

factors      <- c("NGEN", "NFAC", "NVAR", "GLOAD", "FLOAD", "GRHO", "FRHO", "NOBS", "LSKEW", "NCAT")
metrics      <- c("CC", "AB", "RB")
heywood_cut  <- 0.05
omega_cut    <- 0.14
metric_pat   <- paste0("^(", paste(metrics, collapse = "|"), ")_(general|specific)$")

long <- res %>%
  select(ROW, all_of(factors), Method, conv_rate, heywood_rate, matches(metric_pat)) %>%
  pivot_longer(matches(metric_pat), names_to = c("Metric", "Order"), names_sep = "_",
               values_to = "Value") %>%
  mutate(Order = str_to_title(Order),
         clean = !is.na(heywood_rate) & heywood_rate < heywood_cut)

clean <- long %>% filter(clean)

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

heywood_table <- map_df(factors, function(f) {
  res %>%
    group_by(Level = as.character(.data[[f]]), Method) %>%
    summarise(m = mean(heywood_rate, na.rm = TRUE), s = sd(heywood_rate, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(Factor = f, value = paste0(fmt(m), " (", fmt(s), ")"))
}) %>%
  select(Factor, Level, Method, value) %>%
  pivot_wider(names_from = Method, values_from = value)

kable(heywood_table, escape = FALSE, align = "l") %>%
  add_header_above(c(" " = 2, "Heywood case proportion" = ncol(heywood_table) - 2)) %>%
  add_footnote("Values in parentheses are standard deviations.", notation = "none") %>%
  kable_styling(full_width = FALSE, position = "center")

metric_table <- map_df(factors, function(f) {
  full <- long %>%
    group_by(Level = as.character(.data[[f]]), Method, Metric, Order) %>%
    summarise(full = mean(Value, na.rm = TRUE), .groups = "drop")
  nh <- clean %>%
    group_by(Level = as.character(.data[[f]]), Method, Metric, Order) %>%
    summarise(nh = mean(Value, na.rm = TRUE), .groups = "drop")
  left_join(full, nh, by = c("Level", "Method", "Metric", "Order")) %>%
    mutate(Factor = f)
}) %>%
  mutate(value = paste0(fmt(full), " (", fmt(nh), ")")) %>%
  select(Factor, Level, Method, Metric, Order, value) %>%
  pivot_wider(names_from = c(Method, Metric, Order), values_from = value, names_sep = "_")

kable(metric_table, escape = FALSE, align = "l") %>%
  add_footnote(paste0("Values in parentheses exclude conditions with a Heywood case proportion of ",
                      heywood_cut, " or more."), notation = "none") %>%
  kable_styling(full_width = FALSE, position = "center")

run_anova <- function(method, metric, order) {
  d <- clean %>%
    filter(Method == method, Metric == metric, Order == order, is.finite(Value)) %>%
    mutate(across(all_of(factors), as.factor))
  frm   <- as.formula(paste("Value ~ (", paste(factors, collapse = " + "), ")^3"))
  model <- aov(frm, data = d)

  tab <- tidy(model)
  if (!"statistic" %in% names(tab)) tab$statistic <- NA_real_
  if (!"p.value"   %in% names(tab)) tab$p.value   <- NA_real_
  tab <- rename(tab, Effect = term, SS = sumsq, MS = meansq, F = statistic, p = p.value)
  df_residual <- tab$df[tab$Effect == "Residuals"]

  omega <- omega_squared(model, partial = TRUE, ci = 0.95, alternative = "two.sided") %>%
    as.data.frame() %>%
    rename(Effect = Parameter) %>%
    mutate(CI_low  = pmax(coalesce(CI_low, 0), 0),
           CI_high = pmax(coalesce(CI_high, 0), 0))

  left_join(tab, omega, by = "Effect") %>%
    filter(Effect != "Residuals") %>%
    mutate(Method = method, Analysis = paste0(metric, "-", order),
           df_residual = df_residual,
           order = str_count(Effect, ":") + 1,
           F_fmt = fmt(F, 2),
           p_fmt = case_when(p < .001 ~ "< .001", p < .01 ~ "< .01", p < .05 ~ "< .05",
                             TRUE ~ fmt(p)),
           Omega_CI = paste0(fmt(Omega2_partial, 2), " [", fmt(CI_low, 2), ", ",
                             fmt(CI_high, 2), "]")) %>%
    arrange(order, desc(Omega2_partial)) %>%
    select(Method, Analysis, Effect, df, df_residual, SS, MS, F_fmt, p_fmt, Omega_CI,
           Omega2_partial)
}

anova_grid   <- expand.grid(method = unique(clean$Method), metric = metrics,
                            order = c("General", "Specific"), stringsAsFactors = FALSE)
anova_tables <- pmap_df(anova_grid, run_anova) %>%
  mutate(SS = round(SS, 2), MS = round(MS, 2))

kable(select(anova_tables, -Omega2_partial), escape = FALSE, align = "l",
      col.names = c("Method", "Analysis", "Effect", "df1", "df2", "SS", "MS", "F", "p",
                    "partial omega-squared [95% CI]")) %>%
  kable_styling(full_width = FALSE, position = "center")

main_text_table <- anova_tables %>%
  filter(Omega2_partial >= omega_cut) %>%
  select(-Omega2_partial)

kable(main_text_table, escape = FALSE, align = "l",
      col.names = c("Method", "Analysis", "Effect", "df1", "df2", "SS", "MS", "F", "p",
                    "partial omega-squared [95% CI]")) %>%
  collapse_rows(columns = 1:2, valign = "top") %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  footnote(general = paste0("Only effects with partial omega-squared of ", omega_cut,
                            " or more are reported. Complete ANOVA tables are provided in the ",
                            "Supplementary Material."))

writexl::write_xlsx(clean, "gslid_results_noheywood.xlsx")

method_linetype <- c(ML = "dashed", ULS = "dashed")
method_shape    <- c(ML = 16, ULS = 17)

plot_theme <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(size = 9),
        axis.title   = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.text  = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        strip.text   = element_text(size = 10, face = "bold"))

line_plot <- function(d, facet, ylab, ylim, ref = NULL) {
  p <- ggplot(d, aes(x = as.factor(NCAT), y = Value, group = Method, color = Method)) +
    geom_line(aes(linetype = Method), linewidth = 1) +
    geom_point(aes(shape = Method), size = 2.5) +
    facet_grid(as.formula(paste("~ Order +", facet)), labeller = label_both) +
    scale_linetype_manual(values = method_linetype) +
    scale_shape_manual(values = method_shape) +
    coord_cartesian(ylim = ylim) +
    labs(x = "Number of Categories", y = ylab) +
    plot_theme
  if (!is.null(ref)) p <- p + geom_hline(yintercept = ref, color = "black", linewidth = 1,
                                         linetype = "dashed")
  p
}

cc_specific <- clean %>%
  filter(Metric == "CC", Order == "Specific") %>%
  group_by(Method, Order, FLOAD, NCAT) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

ab_specific <- clean %>%
  filter(Metric == "AB", Order == "Specific") %>%
  group_by(Method, Order, NFAC, NCAT) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

print(line_plot(cc_specific, "FLOAD", "Tucker's Congruence Coefficient", c(0.50, 1), ref = 0.95))
print(line_plot(ab_specific, "NFAC", "Absolute Bias", c(0.00, 0.15)))
