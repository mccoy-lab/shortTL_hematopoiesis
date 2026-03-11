# Generate a uniform distribution of telomere lengths
# between the 1st and 99th population percentiles.
# Each value represents the TL at birth for a single simulated individual.

# Draw from a uniform distribution between 8500 and 13400 kb TL
set.seed(804)
lower_bound <- 8500
upper_bound <- 13400
pop_telo_lens <- runif(num_sims,
                       min = lower_bound,
                       max = upper_bound)

# Write pop_telo_lens to file
write.table(pop_telo_lens,
            file = "pop_telo_lens.txt",
            col.names = FALSE,
            row.names = FALSE)