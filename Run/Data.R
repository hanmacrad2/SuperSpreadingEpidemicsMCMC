#****************************************************************
# DATA
#****************************************************************


#****************************************************************
# CANDADIAN DATA
#****************************************************************
data_type = 'Canadian'; seed_count = 1;

canadaX = c(1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 2, 3, 0, 2, 2,
            5, 7, 9, 7, 3, 4, 1, 4, 5, 7, 7, 7, 7, 3, 3, 5, 3, 5, 7, 4, 4, 2,
            3, 1, 1, 1, 0, 0, 2, 1, 3, 2, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1,
            1, 1, 0, 1, 0, 0, 0, 1, 0, 2, 0, 3, 2, 2, 1, 2, 3, 4, 5, 5, 4, 6,
            6, 4, 8, 5, 6, 7, 5, 9, 1, 2, 3, 1, 1, 2, 0, 0, 0, 2, 0, 0, 0, 1)

#************************
#DATA I Extreme - No SS
canada_ss = rep(0, length(canadaX))
sim_data_canadaX1 = list(canadaX, canada_ss)

#************************
#DATA II Extreme SS; [1 0]
canada_bool = canadaX > 1
canada_ss = as.integer(canada_bool)
canada_ns = canadaX - canada_ss
sim_data_canadaX2 = list(canada_ns, canada_ss)

#PLOT
plot.ts(canadaX,
        main = 'SARs (2003) Canadian Outbreak',
        ylab = 'infection count')
