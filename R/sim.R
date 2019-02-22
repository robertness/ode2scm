ode_sim <- function(transition_function, initial_states, times){
  initial_states <- structure(as.numeric(initial_states), names = names(initial_states))
  rates <- attr(transition_function, 'rates')
  rates <- structure(as.numeric(rates), names = names(rates))
  as_tibble(
    deSolve::ode(
      y = initial_states,
      times = times,
      func = transition_function,
      parms = rates
    )
  )
}

sde_sim <- function(transition_function, initial_states, times){
  initial_states <- structure(as.numeric(initial_states), names = names(initial_states))
  t_delta <- times[2] - times[1]
  out <- as_tibble(
    smfsb::simTs(initial_states, times[1], times[length(times)], t_delta, transition_function)
  )
  out$time <- times
  out <- out[, c('time', setdiff(names(out), 'time'))]
  return(out)
}

sim_from_random_seed <- function(seed, initial_states, sample_time, stoc_transition_function){
  initial_states <- structure(as.numeric(initial_states), names = names(initial_states))
  set.seed(seed)
  out <- smfsb::simSample(n=1, initial_states, t0=0, sample_time, stoc_transition_function)
  list(
    Raf=as.numeric(out[1, 'PRaf']),
    Mek=as.numeric(out[1, 'PPMek']),
    Erk=as.numeric(out[1, 'PPErk'])
  )
}
