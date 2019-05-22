mapk_ode_no_double <- function(states, rates){
  transition_function <- function(t, states, rates) {
    with(as.list(c(states, rates)), {
      dE1 <- 0
      dRaf <- -raf_activate * Raf * E1 + raf_deactivate * PRaf
      dPRaf <- raf_activate * Raf * E1 - raf_deactivate * PRaf
      dMek <- -mek_activate * Mek * PRaf + mek_deactivate * PMek
      dPMek <- mek_activate * Mek * PRaf - mek_deactivate * PMek
      dErk <- -erk_activate * Erk * PMek + erk_deactivate * PErk
      dPErk <- erk_activate * Erk * PMek - erk_deactivate * PErk

      list(c(dE1, dRaf, dPRaf, dMek, dPMek, dErk, dPErk))
    })
  }
  attr(transition_function, 'rates') <- rates
  return(transition_function)
}

mapk_sde_no_double <- function(states, rates){
  sde <- list()

  sde$Pre <- matrix(
    c(
      1,	1,	0,	0,	0,	0,	0,
      0,	0,	1,	0,	0,	0,	0,
      0,	0,	1,	1,	0,	0,	0,
      0,	0,	0,	0,	1,	0,	0,
      0,	0,	0,	0,	1,	1,	0,
      0,	0,	0,	0,	0,	0,	1
    ), nrow=6, ncol=7, byrow=T
  )

  sde$Post <- matrix(
    c(
      1,	0,	1,	0,	0,	0,	0,
      0,	1,	0,	0,	0,	0,	0,
      0,	0,	1,	0,	1,	0,	0,
      0,	0,	0,	1,	0,	0,	0,
      0,	0,	0,	0,	1,	0,	1,
      0,	0,	0,	0,	0,	1,	0
    ), nrow=6, ncol=7, byrow=T
  )

  sde$h <- function(states, t, parameters=rates){
    with(as.list(c(states, parameters)), {
      out <- c(
        raf_activate * Raf * E1,
        raf_deactivate * PRaf,
        mek_activate * PRaf * Mek,
        mek_deactivate * PMek,
        erk_activate * PMek * Erk,
        erk_deactivate * PErk
      )
      return(out)
    })
  }
  transition_function <- StepGillespie(sde)
  return(transition_function)
}
