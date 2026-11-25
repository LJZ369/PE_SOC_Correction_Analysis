# Author: Zhang Linjie

# ============================================================
# 0. Load packages and data
# ============================================================

library(tidyverse)
library(broom)
library(gridExtra)
library(multcomp)
library(relaimpo)
library(ggpubr)

# ============================================================
# 1. Linear regression analysis
# ============================================================

# Data preprocessing
soc.data.for.lm <- soc.data |> 
  mutate(
    Separation = factor(Separation, levels = c('no', 'mix', 'yes')),
    method = factor(method),
    delta.y = PE - CK
  )


# ------------------------------------------------------------
# 1.1 Effects of separation status
# ------------------------------------------------------------

# Linear regression
fit.lm <- lm(delta.y ~ Separation * mp.con, data = soc.data.for.lm)

summary(fit.lm)

# Multiple comparisons
glht(
  fit.lm,
  rbind(
    c(0, 0, 0, 0, -1, 1),
    c(0, 0, 0, 0, 1, 0),
    c(0, 0, 0, 0, 0, 1),
    c(0, 0, 0, 1, 1, 0)
  ),
  alternative = c("greater")
) |> summary()

# Visualization
soc.data.for.lm |> 
  ggplot(aes(x = mp.con, y = delta.y, colour = Separation)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = 'lm', se = F) +
  scale_color_manual(
    values = c("#0072B2", "#D55E00", "#EE0000FF")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = 'black', size = 12),
    axis.title = element_text(color = 'black', size = 15)
  ) +
  labs(
    y = 'Apparent SOC change (g/kg)',
    x = 'Added PE (w/w)'
  )


# ------------------------------------------------------------
# 1.2 Effects of measurement methods
# ------------------------------------------------------------

# Linear regression using fully separated samples
fit.lm2 <- lm(
  delta.y ~ method * mp.con,
  data = subset(soc.data.for.lm, Separation == 'yes')
)

summary(fit.lm2)

# Slope test
glht(
  fit.lm2,
  rbind(
    c(0, 0, 0, 1)
  ),
  alternative = c("greater")
) |> summary()

# Visualization
subset(soc.data.for.lm, Separation == 'yes') |> 
  ggplot(aes(x = mp.con, y = delta.y)) +
  geom_point(aes(colour = method), alpha = 0.7) +
  geom_smooth(aes(colour = method), method = 'lm', se = F) +
  scale_color_manual(values = c("gray70", "black")) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = 'black', size = 12),
    axis.title = element_text(color = 'black', size = 15)
  ) +
  labs(
    y = 'Apparent SOC change (g/kg)',
    x = 'Added PE (w/w)'
  )


# ============================================================
# 2. Calculation of the correction coefficient
# ============================================================

# Select studies with exposure durations of no more than 7 days
# and using the potassium dichromate oxidation method
soc.data.pe <- soc.data |> 
  filter(
    Separation == 'yes' &
      mp.time <= 7 &
      method == 'k2cr2o7'
  )

# Calculate theoretical PE-derived SOC and the difference
soc.data.pe <- soc.data.pe |> 
  mutate(
    PE.theory = (CK + 856.2772 * mp.con) / (1 + mp.con),
    delta.PE = PE.theory - PE
  )

# Nonlinear regression
fit <- nls(
  delta.PE ~ k * mp.con / (1 + mp.con),
  data = soc.data.pe,
  start = list(k = 1)
)

summary(fit)

# Confidence interval of the fitted parameter
confint(fit, "k", level = 0.95)

# Calculate model R-squared
y_obs <- soc.data.pe$delta.PE
y_pred <- predict(fit)

R2 <- 1 - sum((y_obs - y_pred)^2) /
  sum((y_obs - mean(y_obs))^2)

R2

# Calculate the correction coefficient
k = 1 - coef(fit)["k"] / 856.2772

k

# Generate fitted values
y.plot <- predict(
  fit,
  newdata = data.frame(
    mp.con = seq(0, 0.07, length.out = 200)
  )
)

plot.df <- data.frame(
  mp.con = seq(0, 0.07, length.out = 200),
  y = y.plot
)

# Visualization
plot.df |> 
  ggplot() +
  geom_line(
    aes(x = mp.con, y = y),
    linewidth = 1.5,
    color = "#0072B2"
  ) +
  geom_point(
    data = soc.data.pe,
    aes(x = mp.con, y = delta.PE),
    size = 3,
    color = "#0072B2",
    alpha = 0.6
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = 'black', size = 12),
    axis.title = element_text(color = 'black', size = 15)
  ) +
  labs(
    y = 'Apparent SOC change (g/kg)',
    x = 'Added PE (w/w)'
  )


# ============================================================
# 3. Correction of SOC changes
# ============================================================

# Exclude partially separated studies
soc.data.pe.corrct <- soc.data |> 
  filter(Separation != 'mix')

# Calculate corrected PE content and SOC changes
soc.data.pe.corrct <- soc.data.pe.corrct |> 
  mutate(
    PE.correct = case_when(
      Separation == 'yes' & method == 'toc' ~
        PE + (PE - 856.2772) * mp.con,
      
      Separation == 'yes' & method == 'k2cr2o7' ~
        PE + (PE - k * 856.2772) * mp.con,
      
      T ~ PE
    )
  ) |> 
  mutate(
    delta.y = PE - CK,
    delta.y.correct = PE.correct - CK
  )


# ============================================================
# 4. Comparison of apparent and corrected SOC changes
# ============================================================

# Paired t-test
t.test(
  soc.data.pe.corrct$delta.y,
  soc.data.pe.corrct$delta.y.correct,
  paired = TRUE
)

# Descriptive statistics
soc.data.pe.corrct$delta.y |> mean()
soc.data.pe.corrct$delta.y |> sd()

soc.data.pe.corrct$delta.y.correct |> mean()
soc.data.pe.corrct$delta.y.correct |> sd()

# Reshape data for visualization
soc.data.pe.corrct |> 
  dplyr::select(delta.y, delta.y.correct) |> 
  mutate(id = row_number()) |> 
  pivot_longer(
    cols = c(delta.y, delta.y.correct),
    names_to = "group",
    values_to = "value"
  ) |> 
  mutate(
    group = factor(
      group,
      levels = c("delta.y", "delta.y.correct")
    )
  ) |> 
  ggplot(aes(x = group, y = value, group = id)) +
  geom_line(
    color = "grey70",
    linewidth = 0.6
  ) +
  geom_point(
    aes(color = group),
    size = 3,
    alpha = 0.3,
    stroke = 1
  ) +
  scale_color_manual(
    values = c(
      "delta.y" = "#0072B2",
      "delta.y.correct" = "#D55E00"
    )
  ) +
  stat_summary(
    aes(group = group),
    fun = mean,
    geom = "point",
    color = "black",
    shape = 18,
    size = 5
  ) +
  stat_summary(
    aes(group = group),
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    color = "black",
    width = 0.15,
    linewidth = 1
  ) +
  annotate(
    'text',
    x = 1.5,
    y = 42,
    label = 'Paired t-test: p < 0.001',
    size = 5
  ) +
  labs(
    x = NULL,
    y = "SOC change (g/kg)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = 'black', size = 12),
    axis.title = element_text(color = 'black', size = 15)
  )


# ============================================================
# 5. Factors influencing PE-driven SOC changes
# ============================================================

# Convert categorical variables into factors
soc.data.pe.corrct$Separation <- as.factor(
  soc.data.pe.corrct$Separation
)

soc.data.pe.corrct$method <- as.factor(
  soc.data.pe.corrct$method
)


# ------------------------------------------------------------
# 5.1 Analysis of corrected SOC changes
# ------------------------------------------------------------

# Linear regression
fit.lm3 <- lm(
  delta.y.correct ~ mp.con + mp.time + method + Separation,
  data = soc.data.pe.corrct
)

summary(fit.lm3)

# Extract regression results
fit.lm3.table.data <- tidy(fit.lm3, conf.int = TRUE) %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  transmute(
    Variable = term,
    Estimate = round(estimate, 3),
    SE = round(std.error, 3),
    `P value` = sprintf("%.3f", p.value),
    Significance
  )

# Display regression results
fit.lm3.table <- ggtexttable(
  fit.lm3.table.data,
  rows = NULL,
  theme = ttheme("light", base_size = 12)
)

fit.lm3.table


# Relative importance analysis
imp.fit.lm3 <- calc.relimp(
  fit.lm3,
  type = "lmg",
  rela = TRUE
)

imp.fit.lm3

# Extract relative importance
imp.fit.lm3.df <- data.frame(
  variable = names(imp.fit.lm3@lmg),
  importance = as.numeric(imp.fit.lm3@lmg) * 100
)

# Visualization
imp.fit.lm3.df |> 
  ggplot(
    aes(
      x = reorder(variable, importance),
      y = importance
    )
  ) +
  geom_col(
    aes(fill = reorder(variable, importance)),
    width = 0.7
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", importance)),
    hjust = 1,
    size = 5
  ) +
  scale_fill_manual(
    values = c("gray70", "gray70", "#0072B2", "#0072B2")
  ) +
  scale_y_reverse(
    expand = expansion(mult = c(0.3, 0))
  ) +
  scale_x_discrete(position = "top") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Relative importance (%)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = 'black', size = 12),
    axis.title = element_text(color = 'black', size = 15),
    legend.position = 'none'
  )


# ------------------------------------------------------------
# 5.2 Analysis of apparent SOC changes
# ------------------------------------------------------------

# Linear regression
fit.lm4 <- lm(
  delta.y ~ mp.con + mp.time + method + Separation,
  data = soc.data.pe.corrct
)

summary(fit.lm4)

# Extract regression results
fit.lm4.table.data <- tidy(fit.lm4, conf.int = TRUE) %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  transmute(
    Variable = term,
    Estimate = round(estimate, 3),
    SE = round(std.error, 3),
    `P value` = sprintf("%.3f", p.value),
    Significance
  )

# Display regression results
fit.lm4.table <- ggtexttable(
  fit.lm4.table.data,
  rows = NULL,
  theme = ttheme("light", base_size = 12)
)

fit.lm4.table


# Relative importance analysis
imp.fit.lm4 <- calc.relimp(
  fit.lm4,
  type = "lmg",
  rela = TRUE
)

imp.fit.lm4

# Extract relative importance
imp.fit.lm4.df <- data.frame(
  variable = names(imp.fit.lm4@lmg),
  importance = as.numeric(imp.fit.lm4@lmg) * 100
)

# Visualization
imp.fit.lm4.df |> 
  ggplot(
    aes(
      x = factor(
        variable,
        levels = c('Separation', 'method', 'mp.time', 'mp.con')
      ),
      y = importance
    )
  ) +
  geom_col(
    aes(
      fill = factor(
        variable,
        levels = c('Separation', 'method', 'mp.time', 'mp.con')
      )
    ),
    width = 0.7
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", importance)),
    hjust = 0,
    size = 5
  ) +
  scale_fill_manual(
    values = c("#0072B2", "gray70", "gray70", "#0072B2")
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.3))
  ) +
  labs(
    x = NULL,
    y = "Relative importance (%)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = 'black', size = 12),
    axis.title = element_text(color = 'black', size = 15),
    legend.position = 'none'
  )
