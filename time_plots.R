library(ggplot2)
library(dplyr)
library(ggplot2)

#Load data frames from files

df_summary <- df %>%
  group_by(N) %>%
  summarise(mean_time = mean(time))

ggplot(df, aes(x = N, y = time)) +
  geom_violin(fill = "gray") +
  geom_boxplot(width = 0.1) +
  geom_text(
    data = df_summary,
    aes(x = N, y = mean_time, label = round(mean_time, 2)),
    vjust = -0.5
  ) +
  labs(
    y = "time (seconds)",
    x = expression(Delta[t])
  ) +
  scale_y_log10() +
  theme_bw()




