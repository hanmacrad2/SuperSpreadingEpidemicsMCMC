#POSTERIOR PREDICTIVE CHECKS

#FUNCTION
zigzag <- function(xs) {
  sum(abs(diff(xs)))
}

OUTER_FOLDER = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_comparison/model_evidence/NZ_DATA_WAIT_21_SUBSET_I/"
num_days = length(data_wait_08_21_sub1)
  
#1. MCMC SAMPLES FROM MODEL
i = 1
model_type = 'BASELINE'; print(model_type)
CURRENT_FOLDER = paste0(OUTER_FOLDER, model_type, '/run_', run, '/')
mcmc_ssir =  readRDS(file = paste0(CURRENT_FOLDER, 'mcmc_', tolower(model_type), '_', i,'.rds'))

#2. SIMULATE DATA FROM MODEL (store in matrix)
matrix_data <- matrix(0, nrow = num_reps, ncol = num_days)

#UPPER BOUNS + ZIG-ZAG
upper_bounds <- apply(matrix_data, 2, quantile, probs = 0.975) |> unlist()
mean_est <- apply(matrix_data, 2, mean) |> unlist()
lower_bounds <- apply(matrix_data, 2, quantile, probs = 0.025) |> unlist()


par(mfrow = c(1, 2))
hist(post_pred_samp_zz, xlim = c(0, max(max(post_pred_samp_zz), zigzag(true_data) * 2)))
abline(v = quantile(post_pred_samp_zz, probs = c(0.025, 0.975)), col = 'red')
abline(v = zigzag(true_data))
