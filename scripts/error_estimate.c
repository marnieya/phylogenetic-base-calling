#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define MAX_LINE_LENGTH 100000
#define MAX_READS 10000
#define MAX_SITES 50000
#define HASH_SIZE 50000

// Data structures
typedef struct {
    double prob[4];  // P_A, P_C, P_G, P_T
} PriorProbs;

typedef struct {
    double lik[4];  // Likelihood for A, C, G, T
} ReadLikelihood;

typedef struct {
    int position;
    int num_reads;
    ReadLikelihood *reads;
} SiteData;

typedef struct {
    int position;
    PriorProbs probs;
} PriorEntry;

// Hash table for prior probabilities
typedef struct {
    PriorEntry *entries;
    int *occupied;
    int size;
} PriorHashTable;

// Global data
static SiteData *sites = NULL;
static int num_sites = 0;
static PriorHashTable *prior_table = NULL;

// Hash table functions
PriorHashTable* create_hash_table(int size) {
    PriorHashTable *table = malloc(sizeof(PriorHashTable));
    table->entries = calloc(size, sizeof(PriorEntry));
    table->occupied = calloc(size, sizeof(int));
    table->size = size;
    return table;
}

void hash_insert(PriorHashTable *table, int position, PriorProbs probs) {
    int hash = position % table->size;
    while (table->occupied[hash]) {
        if (table->entries[hash].position == position) {
            table->entries[hash].probs = probs;
            return;
        }
        hash = (hash + 1) % table->size;
    }
    table->entries[hash].position = position;
    table->entries[hash].probs = probs;
    table->occupied[hash] = 1;
}

int hash_lookup(PriorHashTable *table, int position, PriorProbs *probs) {
    int hash = position % table->size;
    int start = hash;
    while (table->occupied[hash]) {
        if (table->entries[hash].position == position) {
            *probs = table->entries[hash].probs;
            return 1;
        }
        hash = (hash + 1) % table->size;
        if (hash == start) break;
    }
    return 0;
}

void free_hash_table(PriorHashTable *table) {
    free(table->entries);
    free(table->occupied);
    free(table);
}

// Parse PP.txt file
int parse_pp_file(const char *filename, PriorHashTable **table) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open PP file: %s\n", filename);
        return -1;
    }

    *table = create_hash_table(HASH_SIZE);
    char line[MAX_LINE_LENGTH];

    // Skip header
    if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "Error: Empty PP file\n");
        fclose(fp);
        return -1;
    }

    int count = 0;
    while (fgets(line, sizeof(line), fp)) {
        int position;
        PriorProbs probs;
        if (sscanf(line, "%d %lf %lf %lf %lf",
                   &position, &probs.prob[0], &probs.prob[1],
                   &probs.prob[2], &probs.prob[3]) == 5) {
            hash_insert(*table, position, probs);
            count++;
        }
    }

    fclose(fp);
    printf("Loaded %d prior probabilities from %s\n", count, filename);
    return 0;
}

// Clean match string using regex-like pattern matching
void clean_match_string(char *match_str) {
    char *src = match_str;
    char *dst = match_str;

    while (*src) {
        // Pattern 4: ^. (read start marker)
        if (*src == '^' && *(src+1)) {
            src += 2;  // Skip ^ and the quality character
            continue;
        }

        // Pattern 3: $ (read end marker)
        if (*src == '$') {
            src++;
            continue;
        }

        // Pattern 1: -[0-9]+[ACGTNacgtn]+ (deletion)
        if (*src == '-') {
            src++;
            // Read the number
            int del_len = 0;
            while (*src && isdigit(*src)) {
                del_len = del_len * 10 + (*src - '0');
                src++;
            }
            // Skip the deletion bases
            while (del_len > 0 && *src) {
                if (strchr("ACGTNacgtn", *src)) {
                    del_len--;
                    src++;
                } else {
                    break;
                }
            }
            continue;
        }

        // Pattern 2: +[0-9]+[ACGTNacgtn*#]+ (insertion)
        if (*src == '+') {
            src++;
            // Read the number
            int ins_len = 0;
            while (*src && isdigit(*src)) {
                ins_len = ins_len * 10 + (*src - '0');
                src++;
            }
            // Skip the insertion bases
            while (ins_len > 0 && *src) {
                if (strchr("ACGTNacgtn*#", *src)) {
                    ins_len--;
                    src++;
                } else {
                    break;
                }
            }
            continue;
        }

        // Keep other characters (., ,, A, C, G, T, etc.)
        *dst++ = *src++;
    }
    *dst = '\0';
}

// Convert phred score to error probability
double phred_to_error_prob(char phred_char) {
    int phred = (int)phred_char - 33;
    return pow(10.0, -phred / 10.0);
}

// Calculate read likelihoods for all 4 bases
void calculate_read_likelihood(char observed_base, char ref_base, char phred_char, ReadLikelihood *lik) {
    double error_prob = phred_to_error_prob(phred_char);
    double match_prob = 1.0 - error_prob;
    double mismatch_prob = error_prob / 3.0;

    // Determine what base was observed
    char obs_upper = toupper(observed_base);
    char ref_upper = toupper(ref_base);

    // If observed is . or , then it matches the reference
    if (observed_base == '.' || observed_base == ',') {
        obs_upper = ref_upper;
    }

    // Calculate likelihood for each base
    lik->lik[0] = (obs_upper == 'A') ? match_prob : mismatch_prob;  // A
    lik->lik[1] = (obs_upper == 'C') ? match_prob : mismatch_prob;  // C
    lik->lik[2] = (obs_upper == 'G') ? match_prob : mismatch_prob;  // G
    lik->lik[3] = (obs_upper == 'T') ? match_prob : mismatch_prob;  // T
}

// Parse pileup file
int parse_pileup_file(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open pileup file: %s\n", filename);
        return -1;
    }

    sites = malloc(MAX_SITES * sizeof(SiteData));
    num_sites = 0;

    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), fp)) {
        char ref_name[256];
        int position;
        char ref_base;
        int depth;
        char match_str[MAX_LINE_LENGTH];
        char quality_str[MAX_LINE_LENGTH];

        if (sscanf(line, "%s %d %c %d %s %s",
                   ref_name, &position, &ref_base, &depth,
                   match_str, quality_str) != 6) {
            continue;
        }

        // Check if this position exists in prior table
        PriorProbs probs;
        if (!hash_lookup(prior_table, position, &probs)) {
            fprintf(stderr, "Error: Pileup position %d not found in PP file\n", position);
            fclose(fp);
            return -1;
        }

        // Clean the match string
        clean_match_string(match_str);

        int match_len = strlen(match_str);
        int qual_len = strlen(quality_str);

        if (match_len != qual_len) {
            fprintf(stderr, "Warning: Length mismatch at position %d (match=%d, qual=%d)\n",
                    position, match_len, qual_len);
            // Use minimum length
            if (qual_len < match_len) {
                match_len = qual_len;
            }
        }

        if (match_len == 0) {
            continue;  // Skip sites with no reads after cleaning
        }

        // Allocate reads for this site
        sites[num_sites].position = position;
        sites[num_sites].num_reads = match_len;
        sites[num_sites].reads = malloc(match_len * sizeof(ReadLikelihood));

        // Calculate likelihoods for each read
        for (int i = 0; i < match_len; i++) {
            calculate_read_likelihood(match_str[i], ref_base, quality_str[i],
                                     &sites[num_sites].reads[i]);
        }

        num_sites++;
        if (num_sites >= MAX_SITES) {
            fprintf(stderr, "Error: Too many sites (max %d)\n", MAX_SITES);
            break;
        }
    }

    fclose(fp);
    printf("Loaded %d sites with read data from %s\n", num_sites, filename);
    return 0;
}

// Log-sum-exp trick to compute log(sum(exp(x))) stably
double log_sum_exp(double *log_values, int n) {
    // Find maximum value
    double max_val = log_values[0];
    for (int i = 1; i < n; i++) {
        if (log_values[i] > max_val) {
            max_val = log_values[i];
        }
    }

    // Handle case where all values are -infinity
    if (max_val == -INFINITY) {
        return -INFINITY;
    }

    // Compute sum of exp(x - max_val)
    double sum_exp = 0.0;
    for (int i = 0; i < n; i++) {
        sum_exp += exp(log_values[i] - max_val);
    }

    return max_val + log(sum_exp);
}

// Calculate log-likelihood for a given epsilon value (with underflow protection)
double calculate_log_likelihood(double epsilon) {
    double log_lik = 0.0;

    for (int i = 0; i < num_sites; i++) {
        int position = sites[i].position;

        // Get prior probabilities for this site
        PriorProbs probs;
        if (!hash_lookup(prior_table, position, &probs)) {
            continue;  // Should not happen, but be safe
        }

        // Calculate adjusted priors P(b, epsilon)
        double adjusted_probs[4];
        for (int b = 0; b < 4; b++) {
            adjusted_probs[b] = (probs.prob[b] + epsilon / 4.0) / (1.0 + epsilon);
        }

        // For each base b, compute: log P(b,e) + Σ_j log l_ijb
        double log_terms[4];
        for (int b = 0; b < 4; b++) {
            log_terms[b] = log(adjusted_probs[b]);

            // Sum over reads in log space
            for (int j = 0; j < sites[i].num_reads; j++) {
                double lik = sites[i].reads[j].lik[b];
                if (lik > 0) {
                    log_terms[b] += log(lik);
                } else {
                    // If likelihood is 0, this base is impossible
                    log_terms[b] = -INFINITY;
                    break;
                }
            }
        }

        // Use log-sum-exp to compute log(Σ_b exp(log_terms[b]))
        double site_log_lik = log_sum_exp(log_terms, 4);
        log_lik += site_log_lik;
    }

    return log_lik;
}

// Golden section search for optimal epsilon
double golden_section_search(double lower_bound, double upper_bound, double tolerance) {
    const double PHI = 0.618034;       // 1 / golden_ratio
    const double RESPHI = 0.381966;    // 1 - PHI

    double a = lower_bound;
    double b = upper_bound;

    // Initialize two interior points
    double x1 = a + RESPHI * (b - a);
    double x2 = a + PHI * (b - a);

    // Evaluate function at interior points
    double f1 = calculate_log_likelihood(x1);
    double f2 = calculate_log_likelihood(x2);

    printf("Initial search bounds: [%.6f, %.6f]\n", lower_bound, upper_bound);
    printf("Initial evaluations - x1=%.6f (LL=%.6f), x2=%.6f (LL=%.6f)\n", x1, f1, x2, f2);

    int step = 0;
    while ((b - a) > tolerance) {
        step++;

        if (f1 > f2) {
            // Maximum is in [a, x2]
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + RESPHI * (b - a);
            f1 = calculate_log_likelihood(x1);
        } else {
            // Maximum is in [x1, b]
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + PHI * (b - a);
            f2 = calculate_log_likelihood(x2);
        }

        if (step % 5 == 0) {
            printf("Step %d - bounds [%.8f, %.8f], interval=%.8f\n", step, a, b, b-a);
        }
    }

    double final_epsilon = (a + b) / 2.0;
    double final_ll = calculate_log_likelihood(final_epsilon);
    printf("Converged after %d steps\n", step);
    printf("Final log-likelihood: %.6f\n", final_ll);

    return final_epsilon;
}

// Cleanup
void cleanup() {
    if (sites) {
        for (int i = 0; i < num_sites; i++) {
            if (sites[i].reads) {
                free(sites[i].reads);
            }
        }
        free(sites);
    }
    if (prior_table) {
        free_hash_table(prior_table);
    }
}

// Evaluate likelihood on a grid
void evaluate_grid(double min_e, double max_e, int num_points, const char *output_file) {
    FILE *fp = fopen(output_file, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open output file: %s\n", output_file);
        return;
    }

    fprintf(fp, "epsilon,log_likelihood\n");

    double step = (max_e - min_e) / (num_points - 1);
    for (int i = 0; i < num_points; i++) {
        double epsilon = min_e + i * step;
        double ll = calculate_log_likelihood(epsilon);
        fprintf(fp, "%.10f,%.10f\n", epsilon, ll);
    }

    fclose(fp);
    printf("Grid evaluation saved to %s\n", output_file);
}

// Write corrected probabilities to file
void write_corrected_probs(const char *output_file, double epsilon) {
    FILE *fp = fopen(output_file, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open output file: %s\n", output_file);
        return;
    }

    // Write header
    fprintf(fp, "position\tPP_A\tPP_C\tPP_G\tPP_T\n");

    // Iterate through all positions in the hash table
    for (int h = 0; h < prior_table->size; h++) {
        if (prior_table->occupied[h]) {
            int position = prior_table->entries[h].position;
            PriorProbs probs = prior_table->entries[h].probs;

            // Calculate corrected probabilities
            double corrected[4];
            for (int b = 0; b < 4; b++) {
                corrected[b] = (probs.prob[b] + epsilon / 4.0) / (1.0 + epsilon);
            }

            // Write to file
            fprintf(fp, "%d\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    position, corrected[0], corrected[1], corrected[2], corrected[3]);
        }
    }

    fclose(fp);
    printf("Corrected probabilities saved to %s\n", output_file);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <PP_file> <pileup_file> [--grid min max npoints] [--output output_file]\n", argv[0]);
        return 1;
    }

    const char *pp_file = argv[1];
    const char *pileup_file = argv[2];
    const char *output_file = NULL;

    printf("=== Error Rate Estimation ===\n\n");

    // Parse PP file
    if (parse_pp_file(pp_file, &prior_table) != 0) {
        cleanup();
        return 1;
    }

    // Parse pileup file
    if (parse_pileup_file(pileup_file) != 0) {
        cleanup();
        return 1;
    }

    // Check if grid evaluation mode
    if (argc >= 7 && strcmp(argv[3], "--grid") == 0) {
        double min_e = atof(argv[4]);
        double max_e = atof(argv[5]);
        int npoints = atoi(argv[6]);

        printf("\n=== Grid Evaluation Mode ===\n");
        printf("Range: [%.6f, %.6f] with %d points\n", min_e, max_e, npoints);

        evaluate_grid(min_e, max_e, npoints, "grid_output.txt");
        cleanup();
        return 0;
    }

    // Check for output file option
    if (argc >= 5 && strcmp(argv[3], "--output") == 0) {
        output_file = argv[4];
    }

    printf("\n=== Running Golden Section Search ===\n");

    // Optimize epsilon
    double optimal_epsilon = golden_section_search(0.0, 10.0, 1e-6);

    printf("\n=== Results ===\n");
    printf("Optimal epsilon (e): %.10f\n", optimal_epsilon);

    // Write corrected probabilities if output file specified
    if (output_file != NULL) {
        printf("\n=== Writing Corrected Probabilities ===\n");
        write_corrected_probs(output_file, optimal_epsilon);
    }

    cleanup();
    return 0;
}
