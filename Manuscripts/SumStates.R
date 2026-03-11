
# Produce figures and summary statistics that will be reported in the manuscript

library(tidyverse)

# data
list.files("data")
df <- read.csv("data/IFEDDemoData.csv")

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