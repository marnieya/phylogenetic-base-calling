// an explanation of the files in this folder + useful commands

// author: monica

sim_results_verr                                simulation results with a variable error rate throughout sites
sim_results_verr_compare                        simulation results with a constant error rate of 0.005 (0.5%)
get_pileup_from_sim.sh                          simulate raw reads from target fasta

1. Checking intended vs. actual depth in simulations:

        for depth in 1 2 3 4 5 10 20; do ./get_actual_cov_from_pileup.sh 'sim_results_verr' ${depth} 'verr'; done
        for depth in 1 2 3 4 5 10 20; do ./get_actual_cov_from_pileup.sh 'sim_results_verr_compare' ${depth} 'e0.005'; done
