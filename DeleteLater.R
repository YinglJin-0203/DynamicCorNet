
adj_mat <- GetAdjMat(df_sim_temp2 %>% select(-id), 
                     mds_type = "Dynamic")
SplMDS_stress_t(1, )