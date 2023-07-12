#INSPECT SSIR MODEL & LIKELIHOOD

#LIKLIHOOD SSIR
infect_curve_ga1 = GET_INFECT_GAMMA_CURVE(EPI_DATA)
eta_dim = length(EPI_DATA)-1 #epidemic_data[1:(length(epidemic_data)-1)]
eta = EPI_DATA[1:eta_dim]

ssir_params = c(-1.2, 0.15)
LOG_LIKE_SSIR(EPI_DATA, infect_curve_ga1, ssir_params, eta)
