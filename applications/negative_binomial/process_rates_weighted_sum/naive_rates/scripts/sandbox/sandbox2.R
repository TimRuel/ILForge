library(glmmTMB)

fit_nb_phi <- glmmTMB(
  Y ~ 0 + process + offset(log(t)),
  family = nbinom2(),        # variance = mu + phi*mu^2
  dispformula = ~ 0 + process, # φ varies by process
  data = data
)
summary(fit_nb_phi)

theta_hat <- exp(fixef(fit_nb_phi)$cond)  # conditional mean (θ)
phi_hat   <- exp(fixef(fit_nb_phi)$disp)  # process-specific φ
