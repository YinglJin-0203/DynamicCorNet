library(tidyverse)

df_org <- read.csv("Data/IFEDDemoData.csv", check.names = F)
df <- df_org %>% select( -Length, -`Weight for age`, -`Height for age`, 
                        -`Weight for height`, -`Bud bead diameter`) 
write.csv(df, file = "Data/IFEDDemoData.csv",row.names = F )

?write.csv
