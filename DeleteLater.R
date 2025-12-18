
library(tidyverse)

df <- read.csv("data/IFEDDemoData.csv")

df_miss <- df[, c("Week", "ID", "Thyroid.volume")] %>%
  rename(time="Week", id = "ID", var="Thyroid.volume") %>%
  mutate(id = as.factor(id)) %>%
  arrange(time) %>%
  group_by(time) %>%
  summarize(Nmiss = sum(is.na(var)), 
            Pctmiss = sum(is.na(var))/length(var),
            label = paste0(Nmiss, " (", round(100*Pctmiss, 2), "%) ")) %>%
  filter(Nmiss > 20)

df_pair <- df[, c("Week", "ID", "Thyroid.volume", "Ovary.volume")] %>%
  rename(time="Week", id = "ID")

var1 <- "Thyroid.volume"
var2 <- "Ovary.volume"
N <- 136
df_cor <- df_pair %>% 
  group_by(time) %>%
  summarize(cor = cor(.data[[var1]], .data[[var2]], 
                      use = "pairwise.complete.obs"),
            Npair = sum(complete.cases(.data[[var1]], .data[[var2]]))/N)
# plot
df_cor %>% 
  filter(complete.cases(.)) %>%
  ggplot()+
  geom_point(aes(x=time, y=cor, color = Npair, size = Npair))+
  geom_line(aes(x=time, y=cor))+
  labs(title = "Empirical correlation",  y = " ",
       color = "% of complete pairs", size = "% of complete pairs")+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(order = 1), 
         size = guide_legend(order = 1))

df_pair %>%
  pivot_longer(c("Thyroid.volume", "Ovary.volume")) %>%
  filter(complete.cases(.)) %>%
  ggplot(aes(x = time, y=value, colour = name, group=interaction(id, name)))+
  geom_line(alpha = 0.7, linewidth=0.5)+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")
  # scale_x_continuous(breaks = t_brk)+
  labs(title = "Variable trend", x = input$time_var, y = " ")+
  theme(legend.position = "bottom")
