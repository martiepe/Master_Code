library(ggplot2)
library(dplyr)
library(ggplot2)

#Load data frames from files

load("varying_M_estimates.Rda")

df_summary <- df %>%
  group_by(M) %>%
  summarise(mean_time = mean(time))


ggplot(df, aes(x = M, y = time)) +
  geom_violin(fill = "gray") +
  geom_boxplot(width = 0.1) +
  geom_text(
    data = df_summary,
    aes(x = M, y = Inf, label = round(mean_time, 2)),
    vjust = 1.2
  ) +
  labs(
    y = "time (seconds)",
    x = expression(Delta[t])
  ) +
  scale_y_log10() +
  theme_bw()





load("varying_N_estimates.Rda")

df_summary <- df %>%
  group_by(N) %>%
  summarise(mean_time = mean(time))


ggplot(df, aes(x = N, y = time)) +
  geom_violin(fill = "gray") +
  geom_boxplot(width = 0.1) +
  geom_text(
    data = df_summary,
    aes(x = N, y = Inf, label = round(mean_time, 2)),
    vjust = 1.2
  ) +
  labs(
    y = "time (seconds)",
    x = expression(Delta[t])
  ) +
  scale_y_log10() +
  theme_bw()




load("varying_thin_estimates.Rda")

df_summary <- df %>%
  group_by(dt) %>%
  summarise(mean_time = mean(time))


ggplot(df, aes(x = dt, y = time)) +
  geom_violin(fill = "gray") +
  geom_boxplot(width = 0.1) +
  geom_text(
    data = df_summary,
    aes(x = dt, y = Inf, label = round(mean_time, 2)),
    vjust = 1.2
  ) +
  labs(
    y = "time (seconds)",
    x = expression(Delta[t])
  ) +
  scale_y_log10() +
  theme_bw()





load("varying_thin_estimates_fixed_Tmax.Rda")

df_summary <- df %>%
  group_by(dt) %>%
  summarise(mean_time = mean(time))


ggplot(df, aes(x = dt, y = time)) +
  geom_violin(fill = "gray") +
  geom_boxplot(width = 0.1) +
  geom_text(
    data = df_summary,
    aes(x = dt, y = Inf, label = round(mean_time, 2)),
    vjust = 1.2
  ) +
  labs(
    y = "time (seconds)",
    x = expression(Delta[t])
  ) +
  scale_y_log10() +
  theme_bw()
