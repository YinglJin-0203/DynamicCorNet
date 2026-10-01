# Here is the code used for data simulation
set.seed(123)  # for reproducibility

df <- read_rds("Manuscripts/Data/sim_smooth_p10_T10.rds")

lapply(df$obs, dim)
df <- lapply(1:length(df$obs), 
             function(x){data.frame(id = 1:100, time = x, df$obs[[x]])})
df <- bind_rows(df)

# introduce minor missing
for(i in 1:10){
  miss_idx <- sample(nrow(df), round(0.1 * nrow(df)))
  df[miss_idx, paste0("X", i)] <- NA
  
}



# introduce major missing
# at time = 6, X6 and X9 had over 90% of missing
N <- nrow(df[df$time==6, ])
miss_id <- sample(N, round(0.9*N))
df[df$time==6 & df$id %in% miss_id, "X6"] <- NA


miss_id <- sample(N, round(0.9*N))
df[df$time==6 & df$id %in% miss_id, "X9"] <- NA

View(df[df$time==6, ])

write.csv(df, "SampleData/Sim_Example.csv", row.names = F)
