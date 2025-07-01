df <- read.csv("data/AppData.csv")
df %>% filter(time==1) %>% View()

df[, c("id", "time", "TotVol_Ovary")] %>%
  filter(complete.cases(.)) %>%
  ggplot(aes(x=time, y=as.factor(id)))+
  geom_tile()+
  scale_x_continuous(breaks = tid)+
  theme(axis.text.y = element_blank())
  
sort(unique(df$time))
