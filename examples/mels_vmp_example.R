data(riesby, package = "melsvmp")
  
riesby_vmp <- mels_vmp(y = "hamd", 
                         beta_formula = ~ week + endog + endweek, 
                         alpha_formula = ~ endog, 
                         tau_formula = ~ week + endog, 
                         id = "id",
                         data = riesby)
  
summary(riesby_vmp)