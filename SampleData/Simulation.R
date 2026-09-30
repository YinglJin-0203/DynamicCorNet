# Here is the code used for data simulation
set.seed(123)  # for reproducibility

df <- read_rds("Manuscripts/Data/sim_smooth_p10_T10.rds")

lapply(df$obs, dim)
df <- lapply(1:length(df$obs), 
             function(x){data.frame(id = 1:100, time = x, df$obs[[x]])})
df <- bind_rows(df)

# introduce missing
for(i in 1:10){
  miss_idx <- sample(nrow(df), round(0.1 * nrow(df)))
  df[miss_idx, paste0("X", i)] <- NA
  
}

write.csv(df, "Sim_Example.csv", row.names = F)
