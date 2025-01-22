#2.SIMULATE FROM SSE MODEL

#PARAMS
list_params = list(k = rep(0.05, 9), r0 = rep(1.2, 9))

#1. R0, k = 0.1 (R0: 0.9, 1.5, 2.5 )
list_params_sse_r0 = list(r0 = c(0.9, 1.5, 2.5),
                          k = rep(0.1, 3))

#i. r0_param = 0.9; k = 0.1
r0_param = 0.9; k = 0.1

SSE_DATA_r01 = SIMULATE_EPI_SSE(r0 = r0_param, k = k)
plot.ts(SSE_DATA_r01)

#ii. r0_param = 1.5; k = 0.1
r0_param = 1.5; k = 0.1
SSE_DATA_r04 = SIMULATE_EPI_SSE(r0 = r0_param, k = k)
plot.ts(SSE_DATA_r044)

#iii. r0_param = 2.5; k = 0.1
r0_param = 2.5; k = 0.1
SSE_DATA_r07 = SIMULATE_EPI_SSE(r0 = r0_param, k = k)
plot.ts(SSE_DATA_r07)