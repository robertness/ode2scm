# Case Study 1: Mapk Pathway


```r
library(ode2scm)
library(ggplot2)
# Connect Python virtual environment

use_python('/anaconda3/envs/pyro/bin/python',required = TRUE)
source("R/mapk_no_double_phos.R")
```

# Modeling and simulating the pathway

Here are some rates and initial states that we are going to use for the MAPK model.  Note that the initial states determine the range of possible values at any time point including after steady state.  For example, the sum of Mek, PMek, and PPMek is constant.  

It is important that the steady states are not close to the boundary values.  For example, the steady state of PPMek should not be 100 or too close to 100 if Mek + PMek + PPMek == 100.  Otherwise the formulation of the SCM we are using will not work.  So when running sensitivity tests, you want to check the steady state values you get are not on the boundaries.

I use Raf, Mek, and Erk to refer to MAPKKK, MAPKK, and MAPK respectively.


```r
# Exp 1
rates <- list(
  raf_activate=0.1,
  raf_deactivate=0.1,
  mek_activate=0.1,
  mek_deactivate=2.0,
  erk_activate=0.1,
  erk_deactivate=1.0
)

# Exp 2
rates <- list(
  raf_activate=0.2,
  raf_deactivate=0.3,
  mek_activate=0.2,
  mek_deactivate=3.0,
  erk_activate=0.2,
  erk_deactivate=1.5
)
if(FALSE){
  
  # Exp 3
rates <- list(
  raf_activate=0.1,
  raf_deactivate=0.3,
  mek_activate=0.5,
  mek_deactivate=5.0,
  erk_activate=0.3,
  erk_deactivate=4.0
)

}
initial_states <-  list(E1=1, Raf=100, PRaf=0, Mek=100, PMek=0, Erk=100, PErk=0)
```

The following performs the ODE and SDE simulation. See the models in R/mapk.R and simuilation functions for explanation.


```r
times <- seq(0, 30, by = .1)
det_transition_func <- mapk_ode_no_double(initial_states, rates)
ode_out <- ode_sim(det_transition_func, initial_states, times)
ode_out <- cbind(ode_out, times)

# Multiplying rates by 20 to speed up Gillespie
faster_rates <- lapply(rates, `*`, 1)
stoc_transition_func <- mapk_sde_no_double(initial_states, faster_rates)
sde_out <- sde_sim(stoc_transition_func, initial_states, times)
sde_out <- cbind(sde_out, times)
```

The structure of the SCM is based on the analytical solution to steady state.  So here I calculate the steady state values using the steady state analytical solution for these reactions, just to make sure they are the same as the values from the ODE sim (and that the values of the SDE sim are centered around this value.)


```r
g1 <- function(a) a / (a + 1)

# Totals depends on initial states, total Raf is Raf + PRaf, total Mek is Mek + PMek, etc
# Writing out explicitly to keep things simple.
totals <- with(initial_states, {
  list(Raf=Raf+PRaf, Mek=Mek + PMek, Erk=Erk + PErk)
})

E1 <- initial_states$E1
Raf <- totals$Raf * g1(E1 * rates$raf_activate / rates$raf_deactivate)
Mek <- totals$Mek * g1(Raf * rates$mek_activate / rates$mek_deactivate)
Erk <- totals$Erk * g1(Mek * rates$erk_activate / rates$erk_deactivate)

steady_states <- list(Raf=Raf, Mek=Mek, Erk=Erk)
```

Plot everything and make sure they align and are not close to boundaries.


```r
cols <- c(MAP3K = '#CC79A7', MAP2K = '#0072B2', MAPK = '#D55E00')

ggplot() + 
  geom_line(aes(x=ode_out$times, y=ode_out$PRaf, color='MAP3K'), size = 1.2) + 
  geom_line(aes(x=sde_out$times, y=sde_out$PRaf, color='MAP3K')) + 
  geom_abline(aes(color='MAP3K'),intercept = Raf, slope = 0,color = cols['MAP3K'][[1]], linetype="dashed",  size = 1.2) +
  
  geom_line(aes(x=ode_out$times, y=ode_out$PMek, color='MAP2K'), size = 1.2) + 
  geom_line(aes(x=sde_out$times, y=sde_out$PMek, color='MAP2K')) + 
  geom_abline(aes(color='MAP2K'), intercept = Mek, slope = 0, color = cols['MAP2K'][[1]],linetype="dashed",  size = 1.2) +
  
  geom_line(aes(x=ode_out$times, y=ode_out$PErk, color='MAPK'), size = 1.2) + 
  geom_line(aes(x=sde_out$times, y=sde_out$PErk, color='MAPK')) + 
  geom_abline(aes(color='MAPK'), intercept = Erk, slope = 0, color = cols['MAPK'][[1]],linetype="dashed",  size = 1.2) +
  
  #labs(title='Simulation of MapK biochemical Process', x='Time', y='Prolieferated Protein abundence') +
  #ylab("% of Phosphorylated Protein") +
  theme(plot.title=element_text(size=20, face='bold', hjust=0.5),
        axis.text.x = element_text(size=35),
        axis.text.y = element_text(size=35),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom",
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  scale_color_manual("Proteins:", 
                     values=cols)
```

![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-29-1.png)

# How to run a sensitivity analysis for counterfactual inference on SCMs derived from mechanistic models.

We know we can do a counterfactual query on a SCM.  Our goal is to then simulate a counterfactual with the ODE and SDE, and make sure that the counterfactual query on the SCM returns a result that matches the simulation.

Variables for the sensitivity analysis:

* Rates.  Generate these randomly in advance, but only keep those that don't cause any of the proteins to hit the boundaries at steady state.
* ODE vs SDE. The ODE is the base case, the SDE ensures that counterfactual inference is robust enough to deal with stochasticity. 

## Pseudocode for ODE analysis

```
ode_causal_effects = [] # list of floats
scm_causal_effects = [] # list of lists, each element is samples from distribution

for rate_set in rate_sets:
  # With ODE model
  run ODE with rate_set, define raf0, mek0, erk0 as steady states of PRaf, PPMek, and PPErk
  change the Raf activation rate to a valid new number
  run the ODE with updated rate set, define raf1, erk1 as steady state PRaf/PPErk
  append erk1 - erk0 to ode_causal_effects
  
  # With SCM
  observation_model = condition(scm_model, Raf=raf0, Mek=mek0, Erk=erk0
  updated_noise <- infer(observation_model, noise, target=noise)
  counterfactual_model = do(scm_model, Raf=raf1)
  counterfactual_dist = infer(counterfactual_model, updated_noise, target=Erk)
  samples = []
  for i in 1:100:
    append (counterfactual_dist.sample() - raf0) to samples
  append samples to scm_causal_effects
```

## Instance of ODE analysis with one set of rates

The sensitivity analysis would do this across a set of rates.


```r
last_time <- 50
times <- seq(0, last_time, by = .1)
#source_python("mapk_scm_no_double_phos.py")

# Simulate a steady state observation
ode_out_1 <- ode_sim(mapk_ode_no_double(initial_states, rates), initial_states, times)
n <- which(ode_out_1$time == 50)
observation <- list(
  Raf = ode_out_1$PRaf[n],
  Mek = ode_out_1$PMek[n],
  Erk = ode_out_1$PErk[n]
)

# Simulate a counterfactual. In this case, reduce Ras activate to a 1/3 of original value
# Collect Raf outcome and Erk outcome.
intervention_rates <- rates
intervention_rates$raf_activate <- rates$raf_activate / 3
ode_out_2 <- ode_sim(mapk_ode_no_double(initial_states, intervention_rates), initial_states, times)

raf_intervention <- list(Raf=ode_out_2$PRaf[n])

ode_causal_effect <- ode_out_1$PErk[n] - ode_out_2$PErk[n]

scm_causal_effect_samples <- scm_erk_counterfactual(
  rates,
  totals,
  observation,
  raf_intervention
)
```


The SCM produces samples, the ODE provides ground truth.

## Pseudocode for SDE analysis

```
  sde_causal_effects = [] # list of lists
  scm_causal_effects = [] # list of lists, each element is samples from distribution
  for rate_set in rate_sets:
    # With ODE and SDE model
    initialize a set of random seeds as Seeds
    samples = [] # list of floats
    for seed in seeds:
      set_seed(seed)
      run SDE with rate_set, define raf0, mek0, erk0 as steady states of PRaf, PPMek, and PPErk
      change the Raf activation rate to a valid new number
      run the ODE with updated rate set, define raf1 as ODE steady state for PRaf
      run the SDE with updated rate set, define erk1 as SDE steady state for PPErk
      append erk1 - erk0 to samples
    append samples to sde_causal_effects
    
    # With SCM
    observation_model = condition(scm_model, Raf=raf0, Mek=mek0, Erk=erk0
    updated_noise <- infer(observation_model, noise, target=noise)
    counterfactual_model = do(scm_model, Raf=raf1)
    counterfactual_dist = infer(counterfactual_model, updated_noise, target=Erk)
    samples = []
    for i in 1:100:
      append (counterfactual_dist.sample() - raf0) to samples
    append samples to scm_causal_effects
```

## Instance of SDE analysis with one set of rates

The sensitivity analysis would do this across a set of rates.


```r
time_point <- 50

# Simulate a counterfactual with the ODE. In this case, reduce Raf activate to a 1/3 of original value
# Collect Raf outcome and Erk outcome.
intervention_rates <- rates
intervention_rates$raf_activate <- rates$raf_activate / 3
ode_out <- ode_sim(mapk_ode_no_double(initial_states, intervention_rates), initial_states, times)
raf_intervention <- list(Raf=ode_out$PRaf[n])


 sde_causal_effects <- NULL
 scm_causal_effect_samples <- NULL
 seeds <- seq(100000,101000,1)
 ptm <- proc.time()
 for(seed in seeds){
   # Simulate a steady state observation
   set.seed(seed)
   sde_out_1 <- sde_sim(mapk_sde_no_double(initial_states, rates), initial_states, times)
   n <- which(sde_out_1$time == 50)
   observation <- list(
     Raf = sde_out_1$PRaf[n],
     Mek = sde_out_1$PMek[n],
     Erk = sde_out_1$PErk[n]
   )
   set.seed(seed)
   sde_out_2 <- sde_sim(mapk_sde_no_double(initial_states, intervention_rates), initial_states, times)
   sde_causal_effects <- c(sde_causal_effects, sde_out_1$PErk[n] - sde_out_2$PErk[n])
 }
 t = proc.time() - ptm
 
 scm_causal_effect_samples <- c(
   scm_causal_effect_samples,
   scm_erk_counterfactual(
     rates,
     totals,
     observation,
     raf_intervention
   )
 )
```


```r
hist_data <- data.frame("scm" = scm_causal_effect_samples)
sde_data <- data.frame("sde" = sde_causal_effects)
model_color <- c(SCM = rgb(0, 0, 1, 0.5), SDE = rgb(1, 0, 0, 0.5))
ggplot() + 
  geom_histogram(data = hist_data, aes(y = (..count..)/sum(..count..), x=scm, color='SCM', bins=10) , fill=model_color['SCM'][[1]]) +
  geom_histogram(data = sde_data, aes(y = (..count..)/sum(..count..), x=sde, color='SDE', bins=10), fill=model_color['SDE'][[1]]) +
  geom_vline(xintercept = ode_causal_effect, color='#009E73', size=1.5) + 
  #labs(title='Causal Effect', x='Causal Effect(MAPK)', y='Freq (%)') +
  coord_cartesian(ylim=c(0,0.3)) +
  theme(plot.title=element_text(size=20, face='bold', hjust=0.5),
        axis.text.x = element_text(size=35),
        axis.text.y = element_text(size=35),
        axis.title.x = element_text(size=35),
        axis.title.y = element_blank(),
        #legend.text=element_text(size=25),
        #legend.title=element_text(size=25),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom",
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  scale_color_manual("Proteins:", 
                     values=model_color)
#> Warning: Ignoring unknown aesthetics: bins

#> Warning: Ignoring unknown aesthetics: bins
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk unnamed-chunk-32](figure/unnamed-chunk-32-1.png)
