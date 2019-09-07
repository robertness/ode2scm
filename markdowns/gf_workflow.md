# Case Study 2: IGF
```r
library(ode2scm)
library(ggplot2)
# Connect Python virtual environment
#use_virtualenv('venv', require=TRUE)
use_python('/anaconda3/envs/pyro/bin/python', require=TRUE)
```


```r
rate_file <- system.file('growth_factor_sheets/growth_factor', 'Values-Rates.csv', package="ode2scm")
rates_raw <- as.list(read.csv(rate_file))
#> Warning in read.table(file = file, header = header, sep = sep, quote
#> = quote, : incomplete final line found by readTableHeader on '/
#> Library/Frameworks/R.framework/Versions/3.4/Resources/library/ode2scm/
#> growth_factor_sheets/growth_factor/Values-Rates.csv'
initial_states <-  list(
  EGFR=37,
  IGFR=5,
  SOS_inactive=100,
  SOS_active=0,
  Ras_inactive=100,
  Ras_active=0,
  PI3K_inactive=100,
  PI3K_active=0,
  AKT_inactive=100,
  AKT_active=0,
  Raf_inactive=100,
  Raf_active=0,
  Mek_inactive=100,
  Mek_active=0,
  Erk_inactive=100,
  Erk_active=0
)
```


```r
rates <- list(
  SOS_activation_by_EGFR=.01,
  SOS_activation_by_IGFR=.01,
  SOS_deactivation=.5,
  Ras_activation_by_SOS=.01,
  Ras_deactivation=.5,
  PI3K_activation_by_EGFR=.01,
  PI3K_activation_by_IGFR=.01,
  PI3K_activation_by_Ras=.01,
  PI3K_deactivation=.5,
  AKT_activation_by_PI3K=.01,
  AKT_deactivation=.5,
  Raf_activation_by_Ras=.01,
  Raf_deactivation_by_AKT=.01,
  Raf_deactivation_by_phosphotase=.3,
  Mek_activation_by_Raf=.05,
  Mek_deactivation=.5,
  Erk_activation_by_Mek=.05,
  Erk_deactivation=.5
)
```


```r
g <- function(a, b) a / (a + b)

ss <- with(rates, {
      EGFR <- 37
      IGFR <- 5
      SOS_active <- 100 * g(
        SOS_activation_by_EGFR * EGFR + SOS_activation_by_IGFR * IGFR, SOS_deactivation 
      )
      Ras_active <- 100 * g(Ras_activation_by_SOS * SOS_active, Ras_deactivation)
      PI3K_active <- 100 * g(
        PI3K_activation_by_EGFR * EGFR +
        PI3K_activation_by_IGFR * IGFR + 
        PI3K_activation_by_Ras * Ras_active,
        PI3K_deactivation
      )
      AKT_active <- 100 * g(AKT_activation_by_PI3K * PI3K_active, AKT_deactivation)
      Raf_active <- 100 * g(
        Raf_activation_by_Ras * Ras_active,
        Raf_deactivation_by_phosphotase + Raf_deactivation_by_AKT * AKT_active
      )
      Mek_active <- 100 * g(Mek_activation_by_Raf * Raf_active, Mek_deactivation)
      Erk_active <- 100 * g(Erk_activation_by_Mek * Mek_active, Erk_deactivation)
      list(
        EGFR=EGFR,
        IGFR=IGFR,
        SOS_active=SOS_active,
        Ras_active=Ras_active,
        PI3K_active=PI3K_active,
        AKT_active=AKT_active,
        Raf_active=Raf_active,
        Mek_active=Mek_active,
        Erk_active=Erk_active
      )
    }
)
```



```r
slow_rates <- lapply(rates, `/`, 10)
fast_rates <- lapply(rates, `*`, 10)

times <- seq(0, 100, by = .01)

det_transition_func <- gf_ode(initial_states, slow_rates)
ode_out <- ode_sim(det_transition_func, initial_states, times)
#tail(ode_out)[endsWith(names(ode_out), '_active')]

stoc_transition_func <- gf_sde(initial_states, fast_rates)
sde_out <- sde_sim(stoc_transition_func, initial_states, times)

cols <- c(SOS = "#999999", 
          RAS = "#E69F00", 
          PI3K = "#56B4E9", 
          AKT = "#009E73", 
          Raf = "#0072B2", 
          Mek = "#D55E00", 
          Erk = "#CC79A7")
ggplot() + 
  geom_line(aes(x=times, y=ode_out$SOS_active, color='SOS'), size = 1.2) + 
  geom_line(aes(x=times, y=sde_out$SOS_active, color='SOS'), alpha=I(0.4)) + 
  geom_abline(aes(color='SOS'),intercept = ss$SOS_active, slope = 0,color = cols['SOS'][[1]], linetype="dashed",  size = 1.2) +
  
  geom_line(aes(x=times, y=ode_out$Ras_active, color='RAS'), size = 1.2) + 
  geom_line(aes(x=times, y=sde_out$Ras_active, color='RAS'), alpha=I(0.4)) + 
  geom_abline(aes(color='RAS'), intercept = ss$Ras_active, slope = 0, color = cols['RAS'][[1]],linetype="dashed",  size = 1.2) +
  
  geom_line(aes(x=times, y=ode_out$PI3K_active, color='PI3K'), size = 1.2) + 
  geom_line(aes(x=times, y=sde_out$PI3K_active, color='PI3K'), alpha=I(0.4)) + 
  geom_abline(aes(color='PI3K'), intercept = ss$PI3K_active, slope = 0, color = cols['PI3K'][[1]],linetype="dashed",  size = 1.2) +
  
  #labs(title='Simulation of MapK biochemical Process', x='Time', y='Prolieferated Protein abundence') +
  #ylab("% of Phosphorylated Protein") +
  xlab("time") + 
  theme(plot.title=element_text(size=20, face='bold', hjust=0.5),
        axis.text.x = element_text(size=35),
        axis.text.y = element_text(size=35),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=25),
        #legend.title=element_text(size=25),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom",
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  scale_color_manual("Proteins:", 
                     values=cols)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)


```r
ggplot()+
geom_line(aes(x=times, y=ode_out$AKT_active, color='AKT'), size = 1.2) + 
  geom_line(aes(x=times, y=sde_out$AKT_active, color='AKT'), alpha=I(0.4)) + 
  geom_abline(aes(color='AKT'), intercept = ss$AKT_active, slope = 0, color = cols['AKT'][[1]],linetype="dashed",  size = 1.2) +
  
  geom_line(aes(x=times, y=ode_out$Raf_active, color='Raf'), size = 1.2) + 
  geom_line(aes(x=times, y=sde_out$Raf_active, color='Raf'), alpha=I(0.4)) + 
  geom_abline(aes(color='Raf'), intercept = ss$Raf_active, slope = 0, color = cols['Raf'][[1]],linetype="dashed",  size = 1.2) +
  
  geom_line(aes(x=times, y=ode_out$Mek_active, color='Mek'), size = 1.2) + 
  geom_line(aes(x=times, y=sde_out$Mek_active, color='Mek'), alpha=I(0.4)) + 
  geom_abline(aes(color='Mek'), intercept = ss$Mek_active, slope = 0, color = cols['Mek'][[1]],linetype="dashed",  size = 1.2) +
  
  geom_line(aes(x=times, y=ode_out$Erk_active, color='Erk'), size = 1.2) + 
  geom_line(aes(x=times, y=sde_out$Erk_active, color='Erk'), alpha=I(0.4)) + 
  geom_abline(aes(color='Erk'), intercept = ss$Erk_active, slope = 0, color = cols['Erk'][[1]],linetype="dashed",  size = 1.2) +
  
  theme(plot.title=element_text(size=20, face='bold', hjust=0.5),
        axis.text.x = element_text(size=35),
        axis.text.y = element_text(size=35),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=25),
        #legend.title=element_text(size=25),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom",
        panel.background = element_rect(fill = 'white', colour = 'black')) +
  scale_color_manual("Prot:", 
                     values=cols)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

# Counterfactual

The counterfactual intervention on Ras.
 

```r
source_python(system.file("python", "gf_scm.py", package = "ode2scm"))

# Reaches equilibrium by 80
time_point <- 80

# Simulate a counterfactual with the ODE. In this case, reduce Ras rate to a 1/3 of original value
# Collect Ras outcome and Erk outcome.
intervention_rates <- rates
intervention_rates$Ras_activation_by_SOS <- rates$Ras_activation_by_SOS / 6
ode_out_2 <- ode_sim(gf_ode(initial_states, intervention_rates), initial_states, times)
n <- which(ode_out$time == time_point)
ras_intervention <- list(Ras=ode_out$Ras_active[n])


ode_causal_effect <- ode_out$Erk_active[n] - ode_out_2$Erk_active[n]

seed <- 010203

sde_causal_effects <- NULL
scm_causal_effect_samples <- NULL

ptm <- proc.time()

seeds <- seq(100000, 100020, by=1)

ptm <- proc.time()
for(seed in seeds){
  # Simulate a steady state observation
  
  set.seed(seed)
  
  sde_out_1 <- sde_sim(gf_sde(initial_states, rates), initial_states, times)
  
  observation <- list(
    SOS = sde_out_1$SOS_active[n],
    Ras = sde_out_1$Ras_active[n],
    PI3K = sde_out_1$PI3K_active[n],
    AKT = sde_out_1$AKT_active[n],
    Raf = sde_out_1$Raf_active[n],
    Mek = sde_out_1$Mek_active[n],
    Erk = sde_out_1$Erk_active[n]
  )
  
  set.seed(seed)
  
  sde_out_2 <- sde_sim(gf_sde(initial_states, intervention_rates), initial_states, times)
  
  sde_causal_effect <- c(sde_causal_effect, sde_out_1$Erk_active[n] - sde_out_2$Erk_active[n])
  
  totals <- as.list(rep(100, 7))
  names(totals) <- c('SOS', 'Ras', 'PI3K', 'AKT', 'Raf', 'Mek', 'Erk')
  
  ras_intervention <- list(Ras = 30.)
  
}
t = proc.time() - ptm
scm_causal_effect_samples <- c(
  scm_causal_effect_samples,
  scm_ras_erk_counterfactual(
    rates,
    totals,
    observation,
    ras_intervention
  )
)
```


```r
hist_data <- data.frame("scm" = scm_causal_effect_samples)
sde_data <- data.frame("sde" = sde_causal_effect)

model_color <- c(SCM = rgb(0, 0, 1, 0.5), SDE = rgb(1, 0, 0, 0.5))


ggplot() + 
  geom_histogram(binwidth=2.8, data = hist_data, aes(y = (..count..)/sum(..count..), x=scm, color='SCM') , fill=model_color['SCM'][[1]]) +
  #geom_vline(xintercept = 3.638310000000004, color=rgb(0, 0, 1, 0.5), size=1.5) + 
  geom_histogram(binwidth=2.8, data = sde_data, aes(y = (..count..)/sum(..count..), x=sde, color='SDE'), fill=model_color['SDE'][[1]]) +
  geom_vline(xintercept = ode_causal_effect, color='#009E73', size=1.5) + 
  labs(title='Causal Effect', x='Causal Effect(ERK)', y='Freq (%)') +
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
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)
