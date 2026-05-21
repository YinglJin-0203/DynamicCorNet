
# Produce figures and summary statistics that will be reported in the manuscript

library(tidyverse)

# data
list.files("data")
df <- read.csv("data/IFEDDemoData.csv")
# df <- read.csv("data/IFD126005_IFED_GIRLS.csv")

#### Measurement grid ####

data.frame(
  "Physical exam" =     c(1, 2, 4, 6, 8, 12, 16, 20, 24, 28, 32, 36),
  "Sample collection" = c(NA,2, 4, 6, 8, 12, 16, 20, 24, 28, 32, 36),
  "Ultrasound" =        c(1, NA,4, NA,NA,NA, 16, NA, 24, NA, 32, NA)
) %>%
  pivot_longer(1:3, values_to = "Week")%>%
  ggplot(aes(x=Week, y=name, color = ifelse(Week == 6, "gray60", "black")))+
  geom_point(size = 3)+
  # geom_text(data = ~filter(.x, Week == 6), aes(label = "*"), 
  #                         hjust = -0.8, size = 5, color = "gray60")+
  scale_color_identity(guide = "none")+
  scale_x_continuous(breaks = c(1, 2, 4, 6, 8, 12, 16, 20, 24, 28, 32, 36))+
  labs(x="Week", y = "")
ggsave("Manuscripts/CaseStudy/ifed_grid.jpeg", height = 5, width = 10)


#### Identify heighest concentration ####

tgrid <- sort(unique(df$Week))
N <- length(unique(df$ID))

# hormones
df %>%
  select(ID, Week, Bud.bead.volume) %>%
  group_by(ID) %>%
  rename(value = Bud.bead.volume) %>%
  arrange(desc(value)) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(Week) %>%
  summarise(nvalue = n()) #%>%
  mutate(ncum = cumsum(nvalue))
  
# FSH 55 (40.5%) peaked within 4 weeks
# estradiol: more evenly, distributed 2-28 weeks
# Testosterone: 69 (50.7%) peaked between 4 weeks
# ovary: 74 (54.4%) reached largest size between week 16-24
# bude bead: 65 (47.8%) reached alrgest size between week 4



#### Identify peaks ####

##### FSH #####
df_fill <- df %>%
  select(ID, Week, Testosterone) %>%
  rename(id = ID, time = Week, value=Testosterone) %>%
  group_by(id) %>%
  group_modify(~ {
    # Interpolate values onto the common grid
    approx_data <- approx(.x$time, .x$value, xout = tgrid, rule = 2)
    data.frame(time = approx_data$x, value = approx_data$y)
  })
    

df %>% filter(ID == "115602") %>% select(Week, FSH)
df_fill_fsh %>% filter(id == "115602") %>% select(time, value)

# function to find a peak
# A simple version: a point is a peak if it's greater than the point before and after it
is_peak <- function(x) {
  res <- x > lag(x, default = first(x)) & x > lead(x, default = last(x))
  return(res)
}


# Count Peaks per Time Point
df_peak <- df_fill %>%
  group_by(id) %>%
  mutate(peak_fsh = is_peak(value)) %>%
  ungroup() 

# sujects with no peak, 1 peak and more than 1 peak?
id_vec0 <- df_peak %>% group_by(id) %>%
  summarize(have_peak = sum(peak_fsh)) %>%
  filter(have_peak==0) %>% select(id)

id_vec1 <- df_peak %>% group_by(id) %>%
  summarize(have_peak = sum(peak_fsh)) %>%
  filter(have_peak==1) %>% select(id)

id_vec2 <- df_peak %>% group_by(id) %>%
  summarize(have_peak = sum(peak_fsh)) %>%
  filter(have_peak>1) %>% select(id)

df %>% filter(ID %in% unlist(id_vec2)) %>%
  ggplot(aes(x=Week, y=FSH, group = ID))+
  geom_line(na.rm =T)


df_fill_fsh %>% filter(id=="208302")
df %>% filter(ID=="208302")



  
df_peak %>% 
  select(ID, Week, FSH, peak_fsh) %>%
  filter(ID == "123402")


# 5. Visualize the Results
ggplot(peak_summary, aes(x = time, y = total_peaks)) +
  geom_line(color = "steelblue", size = 1) +
  labs(title = "Frequency of Peaks Across Trajectories",
       x = "Time Point",
       y = "Number of Subjects at a Peak") +
  theme_minimal()


#### Measurement grid ####
N <- length(unique(df$IFEDID))
df %>% select(IFEDID, Visit, AgeInDaysIA, AgeInDaysHorm, AgeInDaysUltra) %>%
  group_by(Visit) %>%
  summarise_at(c("AgeInDaysIA", "AgeInDaysHorm", "AgeInDaysUltra"),
               function(x){sum(!is.na(x)) > 0}) %>%
  rename(
    "Week" = Visit,
    "Physical exam" = AgeInDaysIA, 
    "Sample collection" = AgeInDaysHorm,
    "Ultrasound" = AgeInDaysUltra) %>%
  pivot_longer(2:4) %>%
  filter(value) %>%
  ggplot(aes(x=Week, y=name, color = ifelse(Week == 6, "gray60", "black")))+
  geom_point()+
  geom_text(data = ~filter(.x, Week == 6), aes(label = "*"), 
            hjust = -0.8, size = 5, color = "gray60")+
  scale_color_identity(guide = "none")+
  scale_x_continuous(breaks = unique(df$Visit))+
  theme_minimal()+
  labs(y="")
ggsave("Manuscripts/CaseStudy/ifed_grid.jpeg", width = 7, height = 3, bg = "white")


#### Complete pairs ####
# for bivariate analysis

# ovary and thyroid has two complete pairs at week 32 (11th visit), 5 complete pairs at week 24 (9th visit) 
df %>% select(ends_with("volume"), Week) %>%
  group_by(Week) %>%
  summarise(npair = sum(!is.na(Ovary.volume) & !is.na(Bud.bead.volume)))

#### Correlation of smoothed trajectories ####
library(mgcv)


df1 <- df %>% select(ID, Week, Estradiol) 
colnames(df1) <- c("id", "time", "var1")
head(df1)
df1 %>%
  group_by(id) %>% 
  summarise(n = sum(!is.na(var1))) %>% 
  arrange(n)

# no obesrvations: 341902, 656102
# 1 observation: 123402, 460602
tuniq <- sort(unique(df1$time))
span <- 3*mean(diff(tuniq))

df1 %>%
  filter(id==341902) %>%
  loess(var1 ~ time, data = ., span = span)
  

span <- 3*mean(diff(tuniq))
df_smth <- df1 %>%
  group_by(id) %>%
  mutate(n=sum(!is.na(var1))) %>%
  filter(n>3) %>%
  group_modify(~{
    fit1 <- loess(var1 ~ time, data = .x, span = span)
    pred1 <- predict(fit1, newdata = data.frame(time = tuniq))
    # fit2 <- loess(var2 ~ time, data = .x, span = span)
    # pred2 <- predict(fit2, newdata = data.frame(time = tuniq))
    data.frame(time = tuniq, pred1 = pred1)
  })

df_cor <- df_smth %>% 
  group_by(time) %>%
  summarize(cor = cor(.data$pred1, .data$pred2, 
                      use = "pairwise.complete.obs"))


df_pred %>%
  ggplot(aes(x=time, y=pred, group = id)) +
  geom_line()

# LOESS
fit1 <- loess(FSH ~ Week, data = df1)
pred1 <- predict(fit1, newdata = data.frame(Week = sort(tuniq)))
df1$FSH
summary(pred1)
plot(df1$Week, df1$FSH)
lines(sort(tuniq), pred1)
colnames(df1) <- c("id", "time", "var1", "var2")
head(df1)

