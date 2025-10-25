params <- list(
  a1 = 0.2,  b1w = -0.005,
  a2 = -0.5, b21 = 1.2,
  a3 = 0.2,  b31 = -0.8, b32 = 0.9,
  a4 = -0.3, b41 = 0.5,  b42 = -0.5, b43 = 1.1,
  a5 = 10,   b51 = 2,    b52 = -1.5, b53 = 1,  b54 = -1, s5 = 3,
  a6 = 5,    b61 = 1,    b62 = 0.5,  b63 = -1, b64 = 0.8, b65 = 0.5, s6 = 4
)

saveRDS(params, file = "params/params-pop-1.rds")

missingness_params <- list(
  # MNAR models
  nr_int_1 = -0.2,
  nr_beta_1 = c(0.0, 0.2, 0.2, 0.1, 0.1, 0.05),
  
  nr_int_2 = 0.5,
  nr_beta_2 = c(0.1, 0.7, 0.1, 0.1, -0.2, -0.1),
  
  # MAR models
  nr_int_3 = -2.0,
  nr_beta_3 = c(0.6, 0.4, 0.2, 0.1, -0.05),
  
  nr_int_4 = -1.8,
  nr_beta_4 = c(0.3, 0.2, 0.1, 0.1, 0.05),
  
  nr_int_5 = -2.1,
  nr_beta_5 = c(0.7, 0.2, 0.3, -0.3, 0.02),
  
  nr_int_6 = -2.3,
  nr_beta_6 = c(0.4, 0.4, 0.1, -0.1, 0.1)
)

saveRDS(missingness_params, file = "params/missingness_params-1.rds")
