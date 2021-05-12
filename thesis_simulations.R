#### Data ####

library(tidyverse)

contact_matrix <- read_csv("contact_matrix.csv") %>%
  as.matrix()
next_generation_matrix <- read_csv("next_generation_matrix.csv") %>%
  as.matrix()
probabilities_d_l_r <- read_csv("probabilities_d_l_r.csv") %>%
  select(-1) %>%
  as.matrix()

pi <- read_csv("pi.csv") %>%
  as.matrix()

prob_death <- probabilities_d_l_r[1,] %>%
  as.numeric()
prob_long <- probabilities_d_l_r[2,] %>%
  as.numeric()
prob_recovery <- probabilities_d_l_r[3,] %>%
  as.numeric()

#### Final size functions ####

library(pracma)
### Code inspired by the code used in Britton, Ball and Trapman (2020) as ###
### referred to in main text. ###

# Computes the joint probability generating function of the 
# backward branching process for a multi-type SIR model. The
# arguments are a column vector x of the final fractions of
# infected for each group of individuals, a column vector v
# of the fractions of initially vaccinated individuals in
# each group, a column vector pi of the fractions of the 
# population present in each group, a scalar iota of the mean
# infectious period for all groups and a K x K contact matrix 
# where each column represents the group of the infected individual
# making contacts with the individuals from the groups represented
# by the rows.
backprop_pgf <- function(x, v, pi, iota, contact_matrix) {
  exp(iota * contact_matrix %*% (-pi * (1 - x) * (1 - v)))
}

# Computes the extinction probability equation of the backward
# branching process for a multi-type SIR model. The arguments 
# are a column vector x of the final fractions of infected for 
# each group of individuals, a column vector v of the fractions
# of initially vaccinated individuals in each group, a column 
# vector pi of the fractions of the population present in each
# group, a scalar iota of the mean infectious period for all
# groups and a K x K contact matrix where each column represents 
# the group of the infected individual making contacts with the 
# individuals from the groups represented by the rows.
backprop_pgf_root <- function(x, v, pi,iota, contact_matrix) {
  exp(iota * contact_matrix %*% (-pi * (1 - x) * (1 - v))) - x
}

# Computes the final size equation for a multi-type SIR model 
# where R0 > 1.The arguments are a column vector v of the 
# fractions of initially vaccinated individuals in each group, 
# a column vector pi of the fractions of the population present 
# in each group, a scalar iota of the mean infectious period for
# all groups and a K x K contact matrix where each column represents
# the group of the infected individual making contacts with the
# individuals from the groups represented by the rows.
final_fraction <- function(v, pi, iota, contact_matrix) {
  x <- matrix(rep(0, nrow(contact_matrix)), 
              nrow = nrow(contact_matrix), 
              ncol = 1)
  for (i in 1:50) {
    x <- backprop_pgf(x, v, pi, iota, contact_matrix)
  }
  prob_extinction <- fsolve(backprop_pgf_root, 
                            x,
                            v = v,
                            pi = pi,
                            iota = iota, 
                            contact_matrix = contact_matrix)
  z <- rep(0, times = nrow(contact_matrix))
  for (i in 1:nrow(contact_matrix)) {
    z[i] <- 1 - prob_extinction$x[i]
  }
  return(z)
}

# Computes the probabilities of recovery, death and long-term illness
# for each group of individuals given the row vector of final fractions of
# infected for each group.
probabilities <- function(z) {
  rbind(t(z * prob_recovery), 
        t(z * prob_death), 
        t(z * prob_long))
}

# Computes the final percentages of infected, recovered, dead and long-term
# ill of each type in the entire group given the fractions of initially 
# vaccinated as a column vector v, the fractions of individuals in each group
# as a column vector pi, the mean infectious period iota assumed to be equal 
# for all types, and the K x K contact matrix.
final_sizes <- function(v, pi, iota, contact_matrix) {
  # Fractions of each group
  z <- (1 - v) * final_fraction(v, pi, iota, contact_matrix)
  probs_z <- probabilities(z)
  
  return(100 * rbind(t(z), probs_z))
}

#### Computations ####

### Final sizes for an unvaccinated population. ###
v <- as.matrix(rep(0, 10), nrow = 10, ncol = 1)
final_sizes(v, pi, 4, contact_matrix)

### Final sizes when all groups have the same fraction c vaccinated. ###

# Computes the final percentages of infected, recovered, dead, and long-term ill
# among each group when each group has the same fraction c vaccinated. The 
# function takes as arguments the efficacy of the vaccine, the maximum fraction of
# available vaccines c_max where c_max = 0.1, 0.2, ..., 1, the vector pi of the
# fractions of individuals in each group, the mean infectious period iota assumed
# to be the same for each group, and the K x K contact matrix.
vacc_same_frac <- function(efficacy, c_max, pi, iota, contact_matrix) {
  final_inf <- matrix(0, nrow = 10 * c_max, ncol = 10)
  probs_death <- matrix(0, nrow = 10 * c_max, ncol = 10)
  probs_long <- matrix(0, nrow = 10 * c_max, ncol = 10)
  probs_recov <- matrix(0, nrow = 10 * c_max, ncol = 10)
  for (i in seq(0.1, c_max, by = 0.1)) {
    v_values <- as.matrix(c(i,i,i,i,i,i,i,i,i,i), nrow = 10, ncol = 1) 
    z <- final_fraction(v_values * efficacy, pi, iota = 4, contact_matrix)
    
    # Fractions of each group.
    final_inf[i * 10,] <- as.matrix((1 - v_values * efficacy) * z, nrow = 1, ncol = 10)
    probs_recov[i * 10,] <- probabilities(final_inf[i * 10,])[1,]
    probs_death[i * 10,] <- probabilities(final_inf[i * 10,])[2,]
    probs_long[i * 10,] <- probabilities(final_inf[i * 10,])[3,]
  }
  final_sizes <- 100 * rbind(signif(final_inf, 6), 
                             signif(probs_recov, 6), 
                             signif(probs_death, 6), 
                             signif(probs_long, 6))
  return(final_sizes)
}

vacc_same_frac(1, 0.7, pi, iota = 4, contact_matrix)
vacc_same_frac(0.9, 0.7, pi, iota = 4, contact_matrix)

### Final sizes when those with the highest probability of getting infected early ###
### in the outbreak are vaccinated first. ###

# Dominant eigenvector normalized to 1.
dominant_eigen_vec <- eigen(next_generation_matrix)$vectors[,1]
dominant_eigen_vec_norm <- dominant_eigen_vec / sum(dominant_eigen_vec) %>%
  abs() %>%
  c()

# Probabilities of getting infected in early stages of epidemic.
initial_fracs <- abs(dominant_eigen_vec_norm)
prob_infection <- initial_fracs / pi

# Computes the fractions vaccinated in each group when vaccinating in 
# descending order of contribution to the spread according to the entries of 
# the vector prob_infection. The maximum fraction of vaccinated in the population 
# c_max with possible values 0.1, 0.2, ... 1 is the only argument.
v_values_prob_infection <- function(c_max) {
  pi_reordered_prob_infection <- c(pi[10], pi[5], pi[7], pi[2], pi[6],
                                   pi[1], pi[8], pi[3], pi[4], pi[9])
  last_index_vacc <- which(cumsum(pi_reordered_prob_infection) <= c_max) %>%
    length()
  vaccinated <- sum(pi_reordered_prob_infection[1:last_index_vacc])
  frac_last_vaccinated <- (c_max - vaccinated) / pi_reordered_prob_infection[last_index_vacc + 1] %>%
    as.numeric()
  
  v_values_reordered <- c()
  for (i in 1:10) {
    if (i <= last_index_vacc) {
      v_values_reordered[i] <- 1
    }
    else {
      v_values_reordered[i] <- 0
    }
  }
  v_values_reordered[last_index_vacc + 1] <- frac_last_vaccinated
  # Changing back to original order.
  v_values <- c(v_values_reordered[6], 
                v_values_reordered[4], 
                v_values_reordered[8], 
                v_values_reordered[9], 
                v_values_reordered[2], 
                v_values_reordered[5], 
                v_values_reordered[3], 
                v_values_reordered[7], 
                v_values_reordered[10], 
                v_values_reordered[1])
  
  return(v_values)
}

# Table
v_values_prob_inf_tbl <- c()
for (i in 1:7) {
  v_values_prob_inf_tbl <- rbind(v_values_prob_inf_tbl, v_values_prob_infection(i / 10))
}
v_values_prob_inf_tbl %>%
  signif(2)


# Computes the final sizes in each group when vaccinating in order of 
# contribution to the early spread given the maximum fraction c_max of
# vaccinated in the population with possible values 0.1, 0,2, ..., 1, the
# efficacy of the vaccine, the vector pi of the fractions of individuals in 
# each group, the mean infectious period iota assumed to be the same for each 
# group and the K x K contact matrix as arguments.
vacc_order_prob_infection <- function(c_max, efficacy, pi, iota, contact_matrix) {
  v_values <- v_values_prob_infection(c_max)
  # Fractions among each group
  final_sizes <- final_sizes(v_values * efficacy, pi, iota, contact_matrix)
  
  return(final_sizes)
}

# Computes the tables of final sizes in the text, given the efficacy of the
# vaccine, the maximum fraction of vaccinated in the population c_max with 
# possible values 0.1, 0.2, ..., 1, the vector pi of the fractions of individuals
# in each group, the mean infectious period iota assumed to be the same for each 
# group and the K x K contact matrix as arguments.
create_table_prob_infection <- function(efficacy, c_max, pi, iota, contact_matrix) {
  matrix <- matrix(0, nrow = 4 * c_max * 10, ncol = 10)
  for (j in 1:4) {
    fracs <- matrix(0, nrow = c_max * 10, ncol = 10)
    for (i in 1:(c_max * 10)) {
      fracs[i,] <- vacc_order_prob_infection(i / 10, efficacy, pi, iota, contact_matrix)[j,]
    }
    matrix <- rbind(matrix, fracs)
  }
  table <- matrix %>%
    tail(-4 * 10 * c_max) %>%
    signif(4)
  return(table)
}

create_table_prob_infection(1, 0.7, pi, iota = 4, contact_matrix)
create_table_prob_infection(0.9, 0.7, pi, iota = 4, contact_matrix)

### Final sizes when vaccinating in descending order of vulnerability. ###

# Computes the fractions of vaccinated in each group when vaccinating in 
# descending order of vulnerability, given the maximum fraction vaccinated in
# the population c_max with possible values 0.1, 0.2, ..., 1.
v_values_order_vulnerable <- function(c_max) {
  pi_reordered <- c(pi[5], pi[10], pi[4], pi[9], pi[3], 
                               pi[8], pi[2], pi[7], pi[1], pi[6])
  last_index_vacc <- which(cumsum(pi_reordered) <= c_max) %>%
    length()
  vaccinated <- sum(pi_reordered[1:last_index_vacc])
  frac_last_vaccinated <- (c_max-vaccinated) / pi_reordered[last_index_vacc + 1] %>%
    as.numeric()
  
  v_values_reordered <- c()
  for (i in 1:10) {
    if (i <= last_index_vacc) {
      v_values_reordered[i] <- 1
    }
    else {
      v_values_reordered[i] <- 0
    }
  }
  v_values_reordered[last_index_vacc + 1] <- frac_last_vaccinated
  
  # Changing back to original order
  v_values <- c(v_values_reordered[9], 
                v_values_reordered[7], 
                v_values_reordered[5], 
                v_values_reordered[3], 
                v_values_reordered[1], 
                v_values_reordered[10], 
                v_values_reordered[8],
                v_values_reordered[6], 
                v_values_reordered[4], 
                v_values_reordered[2])
  
  return(v_values)
}

# Table
v_values_vuln_tbl <- c()
for (i in 1:7) {
  v_values_vuln_tbl <- rbind(v_values_vuln_tbl, v_values_order_vulnerable(i / 10))
}
v_values_vuln_tbl %>%
  signif(2)

# Computes the final sizes when vaccinating in descending order of vulnerability,
# given the maximum fraction vaccinated in the population c_max with possible 
# values 0.1, 0.2, ..., 1, the efficacy of the vaccine, the row vector pi and 
# the same vector reordered according to the strategy as pi_reordered, the mean 
# infectious period iota, and the K x K contact matrix.
vacc_order_vulnerable <- function(c_max, efficacy, pi, iota, contact_matrix) {
  v_values <- v_values_order_vulnerable(c_max)
  # Fractions among each group
  z <- (1 - v_values * efficacy) * final_fraction(v_values * efficacy, pi, iota = 4, contact_matrix)
  probs_z <- probabilities(z)
  
  return(rbind(z, probs_z))
}

# Computes the tables of final sizes in the text, given the efficacy of the 
# vaccine, the maximum fraction of vaccinated in the population c_max with 
# possible values 0.1, 0.2, ..., 1, the row vector pi and the same vector 
# reordered according to the strategy as pi_reordered, the mean infectious 
# period iota, and the K x K contact matrix.
create_table_vulnerable <- function(efficacy, c_max, pi, iota, contact_matrix) {
  df <- matrix(0, nrow = 4 * c_max * 10, ncol = 10)
  for (j in 1:4) {
    fracs <- matrix(0, nrow = c_max * 10, ncol = 10)
    for (i in 1:(c_max*10)) {
      fracs[i,] <- vacc_order_vulnerable(i / 10, efficacy, pi, iota, contact_matrix)[j,]
    }
    df <- rbind(df, fracs)
  }
  table <- 100 * df %>%
    tail(-4 * 10 * c_max) %>%
    signif(4)
  return(table)
}

create_table_vulnerable(1, 0.7, pi, iota = 4, contact_matrix)
create_table_vulnerable(0.9, 0.7, pi, iota = 4, contact_matrix)

### Computing R_v after vaccination ###

# Computes the effective reproduction number after vaccination for each strategy,
# given the efficacy of the vaccine, the maximum fraction of vaccinated in the
# population c_max with possible values 0.1, 0.2, ..., 1, the strategy as either 
# one of the strings "same frac", "spread", or "vulnerable" for the uniform 
# vaccination, vaccinating in descending contribution to early spread and
# in descending order of vulnerable, respectively.
R_v <- function(efficacy, c_max, strategy) {
  R_vs <- c()
  for (i in seq(0.1, c_max, by = 0.1)) {
    if (strategy == "same frac") {
      v_values <- c(i,i,i,i,i,i,i,i,i,i)
    }
    if (strategy == "spread") {
      v_values <- v_values_prob_infection(i)
    }
    if (strategy == "vulnerable") {
      v_values <- v_values_order_vulnerable(i)
    }
    eigens <- diag(1 - efficacy * v_values) %*% as.matrix(next_generation_matrix) %>%
      eigen()
    R_vs <- c(R_vs, abs(eigens$values[1]))
  }
  return(signif(R_vs, 3))
}

# Table
rbind(
R_v(1, 0.7, strategy = "same frac"),
R_v(1, 0.7, strategy = "spread"),
R_v(1, 0.7, strategy = "vulnerable"),

R_v(0.9, 0.7, strategy = "same frac"),
R_v(0.9, 0.7, strategy = "spread"),
R_v(0.9, 0.7, strategy = "vulnerable"))

### Comparing the final sizes of each strategy with those in the unvaccinated population. ###
unvacc <- data.frame(matrix(NA, nrow = 7, ncol = 10))
v <- as.matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 10, ncol = 1)
z <- final_fraction(v, pi, 4, contact_matrix)
probabilities(z)

for (i in 1:7) {
  unvacc[i,] <- final_fraction(v, pi, 4, contact_matrix)
}

ratio_same_perf <- vacc_same_frac(1, 0.7, pi, iota = 4, contact_matrix)[1:7,] / (100 * unvacc)
ratio_same_aon <- vacc_same_frac(0.9, 0.7, pi, iota = 4, contact_matrix)[1:7,] / (100 * unvacc)

ratio_prob_infection_perf <- create_table_prob_infection(1, 0.7, pi, iota = 4, contact_matrix)[1:7,] / (100 * unvacc)
ratio_prob_infection_aon <- create_table_prob_infection(0.9, 0.7, pi, iota = 4, contact_matrix)[1:7,] / (100 * unvacc)

ratio_vuln_perf <- create_table_vulnerable(1, 0.7, pi, iota = 4, contact_matrix)[1:7,] / (100 * unvacc)
ratio_vuln_aon <- create_table_vulnerable(0.9, 0.7, pi, iota = 4, contact_matrix)[1:7,] / (100 * unvacc)

# Tables
rbind(ratio_same_perf, ratio_same_aon) %>%
  signif(2)

rbind(ratio_prob_infection_perf, ratio_prob_infection_aon) %>%
  signif(2)

rbind(ratio_vuln_perf, ratio_vuln_aon) %>%
  signif(2)
