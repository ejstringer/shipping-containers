# sample from an uniform distribution
stats::runif(1, min = 0, max = 1)

# sample from an exponential distribution
stats::rexp(1)

# sample from a normal distribution
stats::rnorm(1)

# sample from a log normal distribution
stats::rlnorm(1)

calculate_return <- function(alpha) {
  risk_free_rate <- 1.03
  risky_rate <- rnorm(1) * 0.05 + 1
  (1 - alpha) * risk_free_rate + alpha * risky_rate
}


# monte carlo -------------------------------------------------------------

#install.packages('tidyverse')
library(tidyverse)

RUNS <- 1000
DECISION.STEPS <- 12
map(1:10, ~ sample(1:10, 1))
replicate(DECISION.STEPS, runif(1))
simulations <- map(1:RUNS, ~ replicate(DECISION.STEPS, runif(1) %>%
                                       calculate_return())) %>%
  set_names(paste0("sim", 1:RUNS)) %>%
  map(~ accumulate(., ~ .x * .y)) %>%
  map_dfr(~ tibble(value = .x, step = 1:DECISION.STEPS), .id = "simulation")

simulations %>%
  ggplot(aes(x = step, y = value)) +
  geom_line(aes(color = simulation)) +
  theme(legend.position = "none") +
  ggtitle("Simulations of returns from asset allocation")

summary_values <- simulations %>%
  group_by(step) %>%
  summarise(mean_return = mean(value),
            max_return = max(value),
            min_return = min(value)) %>%
  gather("series", "value", -step)

summary_values %>%
  ggplot(aes(x = step, y = value)) +
  geom_line(aes(color = series)) +
  ggtitle("Mean values from simulations")

