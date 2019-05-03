
# Install packages
library(dplyr)
library(lme4)
library(devtools)
#devtools::install_github("pcdjohnson/GLMMmisc")
library(GLMMmisc)
library(MASS)
library(data.table)
library(ggplot2)
library(gridExtra)

# Set working directory
setwd("/Users/samanthaedds/Desktop/Biostat_619")

# Read in dataset
seizuredata = fread("convulsions_data.csv", header = TRUE)

# Make variables factors
seizuredata$treat <- ifelse(seizuredata$treatment == 1, 0, 1)
seizuredata <- dplyr::select(seizuredata, -treatment)
seizuredata <- dplyr::rename(seizuredata, treatment = treat)

# Make variables factors
seizuredata$treatment <- as.factor(seizuredata$treatment)

# Create after treatment variable
seizuredata$after_t <- ifelse(seizuredata$week == 0, 0, 1)
seizuredata$after_t <- as.factor(seizuredata$after_t)

# Create a function to run simulations and bootstrapped confidence intervals
set.seed(135)
# Run glmm from lme4
model <- glmer(convulsions ~ after_t + treatment + (after_t * treatment) + (1 + after_t | patient_id), 
               data = seizuredata,
               family = poisson)
sim.pval <- function(...){
  
  # Fit original data (specifications taken from lme4)
  sim_dat <- sim.glmm(mer.fit = model)
  # Use simulated response data to refit glmm
  fit <- glmer(response ~ after_t + treatment + (after_t * treatment) +  (1 + after_t | patient_id),
               family = "poisson", data = sim_dat)
  
  # Non-parametric bootstrap estimations for confidence interval bounds
  boot.est <- sapply(1:50, function(...) {
    # Fit 50 estimates for one round of resampling
    sim.fit <- update(fit, data = sim.glmm(fit))
    (1 - exp(fixef(sim.fit)[4] - fixef(sim.fit)[2]))
  })
  # Compute 95% confidence interval based on bootstrap estimations
  ci95 <- quantile(boot.est, c(0.025, 0.975))
}

# Calculate 50 confidence intervals
sim.many.pvals <- sapply(1:50, sim.pval)

# Make a decision based on the desired effect size
decision <- ifelse(sim.many.pvals[1, ]<.15 & .15<sim.many.pvals[2, ], 0, 1)

# Count the proportion of rejections over the total, this is power
power <- sum(decision) / length(decision)
power

# Output 95% Confidence Intervals
form = sprintf("ci95.csv")
write.csv(sim.many.pvals, file = form)

# Make a subset to graph
treatment = dplyr::filter(seizuredata, treatment == 1)
control = dplyr::filter(seizuredata, treatment == 0)

# Choose how many
r = 0.002
total = nrow(seizuredata)

# Randomly sample from treatment and control
sub_set = bind_rows(treatment[sample(x = seq_len(nrow(treatment)), floor(r * total), TRUE),], control[sample(seq_len(nrow(control)), floor((r) * total), TRUE),])
sub_seizuredata <- dplyr::filter(seizuredata, patient_id %in% sub_set$patient_id)

# Plot control and treatment
ctrl <- ggplot(sub_seizuredata) + 
  geom_line(data = sub_seizuredata%>% filter(treatment == 0),
            mapping = aes(x = week, y = convulsions, color = patient_id)) +
  labs(x = "Week", y = "Number of Convulsions", title = "Control Group")

treat <- ggplot(sub_seizuredata) + 
  geom_line(data = sub_seizuredata%>% filter(treatment == 1), 
            mapping = aes(x = week, y = convulsions, color = patient_id)) +
  labs(x = "Week", y = "Number of Convulsions", title = "Treatment Group")
 
grid.arrange(ctrl, treat, ncol = 2)
