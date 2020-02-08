mapk_ode <- function(states, rates){
  transition_function <- function(t, states, rates) {
    with(as.list(c(states, rates)), {
      dE1 <- 0
      dRaf <- -raf_activate * Raf * E1 + raf_deactivate * PRaf
      dPRaf <- raf_activate * Raf * E1 - raf_deactivate * PRaf
      dMek <- -mek_activate * Mek * PRaf + mek_deactivate * PMek
      dPMek <- mek_activate * Mek * PRaf +
        mek_deactivate * PPMek -
        mek_deactivate * PMek -
        mek_activate * PMek * PRaf
      dPPMek <- mek_activate * PMek * PRaf - mek_deactivate * PPMek
      dErk <- -erk_activate * Erk * PPMek + erk_deactivate * PErk
      dPErk <- erk_activate * Erk * PPMek +
        erk_deactivate * PPErk -
        erk_deactivate * PErk -
        erk_activate * PErk * PPMek
      dPPErk <- erk_activate * PErk * PPMek - erk_deactivate * PPErk

      list(c(dE1, dRaf, dPRaf, dMek, dPMek, dPPMek, dErk, dPErk, dPPErk))
  })
    }
    attr(transition_function, 'rates') <- rates
    return(transition_function)
}

mapk_sde <- function(states, rates, interventions = NULL){
    sde <- list()

    sde$Pre <- matrix(
        c(
            1, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 1, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1
        ), nrow=10, ncol=9, byrow=T
    )

    sde$Post <- matrix(
        c(
            1, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 0, 1, 0
        ), nrow=10, ncol=9, byrow=T
    )

    sde$h <- function(states, t, parameters=rates, interventions = NULL){
      with(as.list(c(states, parameters, interventions)), {
        out <- c(
          raf_activate * Raf * E1,
          raf_deactivate * PRaf,
          mek_activate * PRaf * Mek,
          mek_deactivate * PMek,
          mek_activate * PRaf * PMek,
          mek_deactivate * PPMek,
          erk_activate * PPMek * Erk,
          erk_deactivate * PErk,
          erk_activate * PPMek * PErk,
          erk_deactivate * PPErk
        )
        return(out)
      })
    }
    transition_function <- StepGillespie(sde)
    return(transition_function)
}
