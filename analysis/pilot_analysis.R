this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(gtable)
library(lme4)
library(tidyverse)
library(lmerTest)
library(bootstrap)
library(ggpubr)
library(stringr)
library(brms)
library(BayesFactor)

`%notin%` <- Negate(`%in%`)
raw_data_path <- "thattracefull-trials.csv"
data<-read.csv(raw_data_path)
data_no_fills <- data %>%
  filter(condition %in% c("null", "that"))%>%
  mutate(block_number = (trial_sequence_total - 1) %/% 4 + 1)

data_summary <- data_no_fills %>%
  group_by(block_number, condition) %>%
  summarise(response = mean(response))

cbPalette = c("#e69d00", "#009e74","#d55e00",  "#cc79a7", "#0071b2")

satiation_plot <- ggplot(data_no_fills, aes(x = block_number, y=response, fill=condition)) +
  geom_point(data=data_summary,alpha=.9) +
  xlab("block number") +
  ylab("average acceptability")+
  geom_smooth(method=lm) +
  scale_fill_manual(values=cbPalette) +
  theme_bw()
satiation_plot



library(dplyr)
library(tidyr)

participant_contrast <- data_no_fills %>%
  group_by(workerid, block_number, condition) %>%
  summarise(mean_response = mean(response), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_response) %>%
  mutate(contrast = null - that)
participant_slopes <- participant_contrast %>%
  group_by(workerid) %>%
  summarise(
    slope = coef(lm(contrast ~ block_number))[2],
    .groups = "drop"
  )
library(ggplot2)

slope_plot <- ggplot(participant_slopes, aes(x = "", y = slope)) +
  geom_violin(fill = "grey85", trim = FALSE) +
  geom_jitter(width = .08, alpha = .7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab(NULL) +
  ylab("Slope of contrast across blocks") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

slope_plot

library(broom)

participant_slopes <- participant_contrast %>%
  group_by(workerid) %>%
  do(tidy(lm(contrast ~ block_number, data = .))) %>%
  ungroup() %>%
  filter(term == "block_number") %>%
  select(workerid, estimate, p.value)%>%
  mutate(direction = ifelse(estimate < 0, "satiating", "non-satiating"))


satiating_participants <- participant_slopes %>%
  filter(direction == "satiating")
##

data_satiating <- participant_contrast %>%
  filter(workerid %in% satiating_participants$workerid)
ggplot(data_satiating,
       aes(x = block_number, y = contrast, group = workerid)) +
  geom_line(alpha = .5) +
  geom_smooth(aes(group = 1), method = "lm", linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Block number") +
  ylab("Contrast (null − that)") +
  theme_bw()

data_with_direction <- participant_contrast %>%
  left_join(participant_slopes %>% select(workerid, direction),
            by = "workerid")
ggplot(data_with_direction,
       aes(x = block_number, y = contrast,
           group = workerid, color = direction)) +
  geom_line(alpha = .5) +
  geom_smooth(aes(group = direction),
              method = "lm", linewidth = 1.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Block number") +
  ylab("Contrast (null − that)") +
  scale_color_manual(values = c("satiating" = "#d55e00",
                                "non-satiating" = "#0071b2")) +
  theme_bw()




### Clustering analysis
# 1. Keep only target conditions and create block number
data_no_fills <- data %>%
  filter(condition %in% c("null", "that")) %>%
  mutate(block_number = (trial_sequence_total - 1) %/% 4 + 1)

# 2. Compute participant contrast by block
participant_contrast <- data_no_fills %>%
  group_by(workerid, block_number, condition) %>%
  summarise(mean_response = mean(response, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_response) %>%
  mutate(contrast = null - that) %>%
  select(workerid, block_number, contrast)

# 3. Wide format: one row per participant
trajectory_data <- participant_contrast %>%
  pivot_wider(
    names_from = block_number,
    values_from = contrast,
    names_prefix = "block_"
  ) %>%
  drop_na()

# 4. Z-score each participant's trajectory across blocks
trajectory_matrix <- trajectory_data %>%
  select(starts_with("block_")) %>%
  as.matrix()

trajectory_matrix <- t(apply(trajectory_matrix, 1, scale))
trajectory_matrix <- trajectory_matrix %>%
  filter(block_1 )

# remove weird dimensions that sometimes appear
trajectory_matrix <- matrix(
  trajectory_matrix,
  nrow = nrow(trajectory_data),
  dimnames = list(trajectory_data$workerid, colnames(trajectory_data %>% select(starts_with("block_"))))
)

# 5. Choose k by average silhouette width
sil_width <- sapply(2:4, function(k) {
  km <- kmeans(trajectory_matrix, centers = k, nstart = 50)
  ss <- silhouette(km$cluster, dist(trajectory_matrix))
  mean(ss[, 3])
})

sil_df <- data.frame(
  k = 2:4,
  silhouette = sil_width
)

ggplot(sil_df, aes(x = k, y = silhouette)) +
  geom_line() +
  geom_point(size = 3) +
  theme_bw() +
  xlab("Number of clusters") +
  ylab("Average silhouette width")

optimal_k <- sil_df$k[which.max(sil_df$silhouette)]
optimal_k

# 6. Run k-means with chosen k
set.seed(123)
km <- kmeans(trajectory_matrix, centers = optimal_k, nstart = 50)

trajectory_data <- trajectory_data %>%
  mutate(cluster = factor(km$cluster))

# 7. Join cluster labels back to long data
data_clustered <- participant_contrast %>%
  left_join(trajectory_data %>% select(workerid, cluster), by = "workerid")

# 8. Plot participant trajectories by cluster
cluster_plot <- ggplot(
  data_clustered,
  aes(x = block_number, y = contrast, group = workerid, color = cluster)
) +
  geom_line(alpha = .35) +
  geom_smooth(aes(group = cluster), method = "lm", se = FALSE, linewidth = 1.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Block number") +
  ylab("Contrast (null - that)") +
  theme_bw()

cluster_plot