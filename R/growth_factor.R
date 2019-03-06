pre_file <- system.file('growth_factor_sheets/growth_factor', 'Values-Pre.csv', package="ode2scm")
post_file <- system.file('growth_factor_sheets/growth_factor', 'Values-Post.csv', package="ode2scm")
PRE <- as.matrix(read.csv(pre_file, header = TRUE))
POST <- as.matrix(read.csv(post_file, header = TRUE))

gf_ode <- function(states, rates){
  transition_function <- function(t, states, rates) {
    with(as.list(c(states, rates)), {
      dEGFR <- 0
      dIGFR <- 0
      dSOS_inactive <- -SOS_activation_by_EGFR * SOS_inactive * EGFR + -SOS_activation_by_IGFR * SOS_inactive * IGFR + SOS_deactivation * SOS_active
      dSOS_active <- SOS_activation_by_EGFR * SOS_inactive * EGFR + SOS_activation_by_IGFR * SOS_inactive * IGFR + -SOS_deactivation * SOS_active
      dRas_inactive <- -Ras_activation_by_SOS * Ras_inactive * SOS_active + Ras_deactivation * Ras_active
      dRas_active <- Ras_activation_by_SOS * Ras_inactive * SOS_active + -Ras_deactivation * Ras_active
      dPI3K_inactive <- -PI3K_activation_by_EGFR * PI3K_inactive * EGFR + -PI3K_activation_by_IGFR * PI3K_inactive * IGFR + -PI3K_activation_by_Ras * PI3K_inactive * Ras_active + PI3K_deactivation * PI3K_active
      dPI3K_active <- PI3K_activation_by_EGFR * PI3K_inactive * EGFR + PI3K_activation_by_IGFR * PI3K_inactive * IGFR + PI3K_activation_by_Ras * PI3K_inactive * Ras_active + -PI3K_deactivation * PI3K_active
      dAKT_inactive <- -AKT_activation_by_PI3K * AKT_inactive * PI3K_active + AKT_deactivation * AKT_active
      dAKT_active <- AKT_activation_by_PI3K * AKT_inactive * PI3K_active + -AKT_deactivation * AKT_active
      dRaf_inactive <- -Raf_activation_by_Ras * Raf_inactive * Ras_active + Raf_deactivation_by_phosphotase * Raf_active + Raf_deactivation_by_AKT * AKT_active * Raf_active
      dRaf_active <- Raf_activation_by_Ras * Raf_inactive * Ras_active + -Raf_deactivation_by_phosphotase * Raf_active + -Raf_deactivation_by_AKT * AKT_active * Raf_active
      dMek_inactive <- -Mek_activation_by_Raf * Mek_inactive * Raf_active + Mek_deactivation * Mek_active
      dMek_active <- Mek_activation_by_Raf * Mek_inactive * Raf_active - Mek_deactivation * Mek_active
      dErk_inactive <- -Erk_activation_by_Mek * Erk_inactive * Mek_active + Erk_deactivation * Erk_active
      dErk_active <- Erk_activation_by_Mek * Erk_inactive * Mek_active - Erk_deactivation * Erk_active

      list(c(dEGFR, dIGFR, dSOS_inactive, dSOS_active, dRas_inactive, dRas_active, dPI3K_inactive, dPI3K_active, dAKT_inactive, dAKT_active, dRaf_inactive, dRaf_active, dMek_inactive, dMek_active, dErk_inactive, dErk_active))
    })
  }
  attr(transition_function, 'rates') <- rates
  return(transition_function)
}

gf_sde <- function(states, rates){
  sde <- list()

  sde$Pre <- PRE
  sde$Post <- POST

  sde$h <- function(states, t, parameters=rates){
    with(as.list(c(states, parameters)), {
      out <- c(
        SOS_activation_by_EGFR * SOS_inactive * EGFR,
        SOS_activation_by_IGFR * SOS_inactive * IGFR,
        SOS_deactivation * SOS_active,
        Ras_activation_by_SOS * Ras_inactive * SOS_active,
        Ras_deactivation * Ras_active,
        PI3K_activation_by_EGFR * PI3K_inactive * EGFR,
        PI3K_activation_by_IGFR * PI3K_inactive * IGFR,
        PI3K_activation_by_Ras * PI3K_inactive * Ras_active,
        PI3K_deactivation * PI3K_active,
        AKT_activation_by_PI3K * AKT_inactive * PI3K_active,
        AKT_deactivation * AKT_active,
        Raf_activation_by_Ras * Raf_inactive * Ras_active,
        Raf_deactivation_by_phosphotase * Raf_active,
        Raf_deactivation_by_AKT * AKT_active * Raf_active,
        Mek_activation_by_Raf * Mek_inactive * Raf_active,
        Mek_deactivation * Mek_active,
        Erk_activation_by_Mek * Erk_inactive * Mek_active,
        Erk_deactivation * Erk_active
      )
      return(out)
    })
  }
  transition_function <- StepGillespie(sde)
  return(transition_function)
}
