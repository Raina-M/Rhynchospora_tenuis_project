# Define the expected proportions (Mendelian 1:1:1:1 ratio)
expected_probs <- c(0.25, 0.25, 0.25, 0.25)

# --- Interpreting the Output ---
# X-squared: The calculated Chi-Square statistic
# df: Degrees of freedom (3)
# p-value: If p < 0.05, reject Mendelian segregation. 
#          If p > 0.05, the data follows Mendelian laws.

# ----------------- REC ---------------------
# Observed genotype combinations of chr1_chr2 in pollens
observed_counts <- c(h1_h1 = 385, h1_h2 = 511, h2_h1 = 710, h2_h2 = 208)

# Chi-Square Goodness-of-Fit Test
chisq.test(x = observed_counts, p = expected_probs)

# Chi-squared test for given probabilities
# 
# data:  observed_counts
# X-squared = 295.61, df = 3, p-value < 2.2e-16


# ----------------- JGV-16 ---------------------
observed_counts <- c(h1_h1 = 229, h1_h2 = 301, h2_h1 = 177, h2_h2 = 197)

chisq.test(x = observed_counts, p = expected_probs)

# Chi-squared test for given probabilities
# 
# data:  observed_counts
# X-squared = 39.274, df = 3, p-value = 1.518e-08


# ----------------- JGV-17 ---------------------
observed_counts <- c(h1_h1 = 194, h1_h2 = 279, h2_h1 = 146, h2_h2 = 182)

chisq.test(x = observed_counts, p = expected_probs)

# Chi-squared test for given probabilities
# 
# data:  observed_counts
# X-squared = 47.524, df = 3, p-value = 2.688e-10


# ----------------- PECP-47 ---------------------
observed_counts <- c(h1_h1 = 920, h1_h2 = 514, h2_h1 = 399, h2_h2 = 67)

chisq.test(x = observed_counts, p = expected_probs)

# Chi-squared test for given probabilities
# 
# data:  observed_counts
# X-squared = 782.71, df = 3, p-value < 2.2e-16


# ----------------- PECP-48 ---------------------
observed_counts <- c(h1_h1 = 425, h1_h2 = 217, h2_h1 = 197, h2_h2 = 36)

chisq.test(x = observed_counts, p = expected_probs)

# Chi-squared test for given probabilities
# 
# data:  observed_counts
# X-squared = 349.32, df = 3, p-value < 2.2e-16