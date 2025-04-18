// hyperfm.c: 2-way hypergraph partitioning using FM refinement
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>  // for INT_MAX in multi-start fallback

typedef struct {
    int *vertices;
    int size;
    int weight;
    int part_count[2];
} Net;

typedef struct {
    int *nets;
    int degree;
    int weight;      // vertex weight
    int part;
    int gain;
    int locked;
    int bucket_prev, bucket_next;
    int bucket_idx;
} Vertex;

int nv, nm, fmt;
Net *nets;
Vertex *vertices;
int max_degree;
// balance constraints bounds (sum of vertex weights)
long long total_weight;
int weight_lb, weight_ub;
// (allowed_diff unused after weight-based balance)
int allowed_diff;
int *neighbor_list;
char *touched;
int **bucket_heads;

// compute absolute
static inline int iabs(int x) { return x < 0 ? -x : x; }
// compute difference in seconds between two timespecs
static inline double diff_ts(const struct timespec *start, const struct timespec *end) {
    return (double)(end->tv_sec - start->tv_sec)
         + (double)(end->tv_nsec - start->tv_nsec) * 1e-9;
}

// comparator for greedy initial partition: sort vertices by descending degree
static int cmp_degree(const void *a, const void *b) {
    int va = *(const int*)a;
    int vb = *(const int*)b;
    return vertices[vb].degree - vertices[va].degree;
}
// comparator for sorting vertices by ascending weight
static int cmp_weight_asc(const void *a, const void *b) {
    int va = *(const int*)a;
    int vb = *(const int*)b;
    return vertices[va].weight - vertices[vb].weight;
}

// recompute gain for vertex u
int recompute_gain(int u) {
    int g = 0;
    int p = vertices[u].part;
    int q = 1 - p;
    for(int i = 0; i < vertices[u].degree; i++) {
        int e = vertices[u].nets[i];
        int cp = nets[e].part_count[p];
        int cq = nets[e].part_count[q];
        if (cp == 1 && cq > 0) g += nets[e].weight;
        if (cq == 0)        g -= nets[e].weight;
    }
    return g;
}

// remove vertex u from its bucket
void bucket_remove(int u) {
    int pid = vertices[u].part;
    int idx = vertices[u].bucket_idx;
    int prev = vertices[u].bucket_prev;
    int next = vertices[u].bucket_next;
    if (prev >= 0) vertices[prev].bucket_next = next;
    else bucket_heads[pid][idx] = next;
    if (next >= 0) vertices[next].bucket_prev = prev;
    vertices[u].bucket_prev = vertices[u].bucket_next = -1;
}

// insert vertex u into its bucket (at head)
void bucket_insert(int u) {
    int pid = vertices[u].part;
    int idx = vertices[u].bucket_idx;
    int head = bucket_heads[pid][idx];
    vertices[u].bucket_next = head;
    vertices[u].bucket_prev = -1;
    if (head >= 0) vertices[head].bucket_prev = u;
    bucket_heads[pid][idx] = u;
}

// check if moving u is legal wrt weight-balance
int is_move_legal(int u, long long *part_weight) {
    int pid = vertices[u].part;
    long long w0 = part_weight[0];
    long long w1 = part_weight[1];
    long long w = vertices[u].weight;
    long long new_w0 = w0, new_w1 = w1;
    if (pid == 0) { new_w0 -= w; new_w1 += w; }
    else          { new_w1 -= w; new_w0 += w; }
    return (new_w0 >= weight_lb && new_w0 <= weight_ub
         && new_w1 >= weight_lb && new_w1 <= weight_ub);
}

// find best move from partition pid
// find best move in partition pid given current part weights
int find_best_move(int pid, long long *part_weight) {
    int offset = max_degree;
    for(int idx = 2*offset; idx >= 0; idx--) {
        int u = bucket_heads[pid][idx];
        while(u >= 0) {
            if (!vertices[u].locked && is_move_legal(u, part_weight))
                return u;
            u = vertices[u].bucket_next;
        }
    }
    return -1;
}

// build net part counts
void build_net_counts() {
    for(int e = 0; e < nm; e++) {
        nets[e].part_count[0] = nets[e].part_count[1] = 0;
    }
    for(int v = 0; v < nv; v++) {
        int p = vertices[v].part;
        for(int i = 0; i < vertices[v].degree; i++) {
            int e = vertices[v].nets[i];
            nets[e].part_count[p]++;
        }
    }
}

// compute current cutsize
int compute_cut() {
    int cut = 0;
    for(int e = 0; e < nm; e++) {
        if (nets[e].part_count[0] > 0 && nets[e].part_count[1] > 0)
            cut += nets[e].weight;
    }
    return cut;
}

// build buckets and initial gains
void init_buckets() {
    int bucket_len = 2*max_degree + 1;
    for(int pid = 0; pid < 2; pid++) {
        for(int i = 0; i < bucket_len; i++)
            bucket_heads[pid][i] = -1;
    }
    for(int v = 0; v < nv; v++) {
        vertices[v].locked = 0;
        vertices[v].gain = recompute_gain(v);
        int idx = vertices[v].gain + max_degree;
        vertices[v].bucket_idx = idx;
        bucket_insert(v);
    }
}

// single FM refinement pass
int fm_refine() {
    // allocate static buffers to avoid repeated malloc/free
    static int *move_list = NULL;
    static int *gain_list = NULL;
    static long long *part_weight = NULL;
    if (!move_list) move_list = malloc(nv * sizeof(int));
    if (!gain_list) gain_list = malloc(nv * sizeof(int));
    if (!part_weight) part_weight = malloc(2 * sizeof(long long));
    part_weight[0] = part_weight[1] = 0;
    for (int v = 0; v < nv; v++)
        part_weight[vertices[v].part] += vertices[v].weight;

    build_net_counts();
    init_buckets();

    int best_gain = 0, gain_sum = 0, best_move_num = 0;
    int moves = 0;
    for(moves = 0; moves < nv; moves++) {
        int v0 = find_best_move(0, part_weight);
        int v1 = find_best_move(1, part_weight);
        int best_v = -1;
        if (v0 < 0 && v1 < 0) break;
        else if (v0 < 0) best_v = v1;
        else if (v1 < 0) best_v = v0;
        else best_v = (vertices[v0].gain >= vertices[v1].gain ? v0 : v1);

        int g = vertices[best_v].gain;
        move_list[moves] = best_v;
        gain_list[moves] = g;

        // apply move
        int oldp = vertices[best_v].part;
        int newp = 1 - oldp;
        vertices[best_v].locked = 1;
        bucket_remove(best_v);
        vertices[best_v].part = newp;
        // update part weights
        part_weight[oldp] -= vertices[best_v].weight;
        part_weight[newp] += vertices[best_v].weight;

        // update net counts
        for(int i = 0; i < vertices[best_v].degree; i++) {
            int e = vertices[best_v].nets[i];
            nets[e].part_count[oldp]--;
            nets[e].part_count[newp]++;
        }

        // update neighbors
        int nlist = 0;
        for(int i = 0; i < vertices[best_v].degree; i++) {
            int e = vertices[best_v].nets[i];
            for(int j = 0; j < nets[e].size; j++) {
                int u = nets[e].vertices[j];
                if (u == best_v || vertices[u].locked) continue;
                if (!touched[u]) {
                    touched[u] = 1;
                    neighbor_list[nlist++] = u;
                }
            }
        }
        for(int i = 0; i < nlist; i++) {
            int u = neighbor_list[i];
            touched[u] = 0;
            bucket_remove(u);
            vertices[u].gain = recompute_gain(u);
            vertices[u].bucket_idx = vertices[u].gain + max_degree;
            bucket_insert(u);
        }

        gain_sum += g;
        if (gain_sum > best_gain) {
            best_gain = gain_sum;
            best_move_num = moves + 1;
        }
    }

    // revert moves beyond best
    for(int i = moves - 1; i >= best_move_num; i--) {
        int v = move_list[i];
        vertices[v].part = 1 - vertices[v].part;
    }

    // buffers persist for reuse
    return best_gain;
}

int main(int argc, char **argv) {
    // parse command line
    // usage: prog [-p random|greedy] <input.hgr> <epsilon>
    // parse options: -p for initial partition method, -r for number of runs
    int initial_method = 1; // 0=random, 1=greedy default
    int num_runs = 1;
    int argi = 1;
    // -p <method>
    if (argi < argc && strcmp(argv[argi], "-p") == 0) {
        if (argi + 1 >= argc) {
            fprintf(stderr, "Missing method for -p\n");
            fprintf(stderr, "Usage: %s [-p random|greedy] [-r runs] <input.hgr> <epsilon>\n", argv[0]);
            return 1;
        }
        if (strcmp(argv[argi+1], "random") == 0) initial_method = 0;
        else if (strcmp(argv[argi+1], "greedy") == 0) initial_method = 1;
        else {
            fprintf(stderr, "Unknown partition method: %s\n", argv[argi+1]);
            fprintf(stderr, "Usage: %s [-p random|greedy] [-r runs] <input.hgr> <epsilon>\n", argv[0]);
            return 1;
        }
        argi += 2;
    }
    // -r <num_runs>
    if (argi < argc && strcmp(argv[argi], "-r") == 0) {
        if (argi + 1 >= argc) {
            fprintf(stderr, "Missing count for -r\n");
            fprintf(stderr, "Usage: %s [-p random|greedy] [-r runs] <input.hgr> <epsilon>\n", argv[0]);
            return 1;
        }
        num_runs = atoi(argv[argi+1]);
        if (num_runs < 1) num_runs = 1;
        argi += 2;
    }
    if (argc - argi < 2) {
        fprintf(stderr, "Usage: %s [-p random|greedy] [-r runs] <input.hgr> <epsilon>\n", argv[0]);
        return 1;
    }
    const char *infile = argv[argi];
    double epsilon = atof(argv[argi+1]);
    // seed random generator for multi-start
    srand((unsigned)time(NULL));
    // start timing
    struct timespec ts_start, ts_after_read, ts_after_compute, ts_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    FILE *f = fopen(infile, "r");
    if (!f) { perror("fopen"); return 1; }
    // read header
    char *line = NULL;
    size_t llen = 0;
    // Read header: hMetis-like but input is "#nets #vertices [fmt]"
    if (getline(&line, &llen, f) <= 0) {
        fprintf(stderr, "Failed to read header\n");
        return 1;
    }
    fmt = 0;
    {
        int cnt = sscanf(line, "%d %d %d", &nm, &nv, &fmt);
        if (cnt < 2) { fprintf(stderr, "Bad header\n"); return 1; }
        if (cnt < 3) fmt = 0;
    }
    // read vertex weights if any
    int *vweight = malloc(nv * sizeof(int));
    if (!vweight) { perror("malloc"); return 1; }
    if (fmt & 1) {
        for (int i = 0; i < nv; i++) {
            if (getline(&line, &llen, f) <= 0) {
                fprintf(stderr, "Unexpected EOF reading vertex weight %d\n", i);
                return 1;
            }
            vweight[i] = atoi(line);
        }
    } else {
        for (int i = 0; i < nv; i++) vweight[i] = 1;
    }
    // read nets
    nets = malloc(nm * sizeof(Net));
    char *tok;
    for(int e = 0; e < nm; e++) {
        if (getline(&line, &llen, f) <= 0) {
            fprintf(stderr, "Unexpected EOF reading net %d\n", e);
            return 1;
        }
        int cap = 8, sz = 0;
        int *verts = malloc(cap * sizeof(int));
        int w = 1;
        tok = strtok(line, " \t\n");
        if (fmt & 2) {
            if (!tok) continue;
            w = atoi(tok);
            tok = strtok(NULL, " \t\n");
        }
        while(tok) {
            int v = atoi(tok) - 1;
            if (v < 0 || v >= nv) {
                fprintf(stderr, "Vertex id out of range: %d\n", v+1);
                return 1;
            }
            if (sz >= cap) {
                cap *= 2;
                verts = realloc(verts, cap * sizeof(int));
            }
            verts[sz++] = v;
            tok = strtok(NULL, " \t\n");
        }
        nets[e].size = sz;
        nets[e].weight = w;
        nets[e].vertices = malloc(sz * sizeof(int));
        memcpy(nets[e].vertices, verts, sz * sizeof(int));
        free(verts);
    }
    fclose(f);
    // mark end of I/O (read)
    clock_gettime(CLOCK_MONOTONIC, &ts_after_read);

    // build vertex incidence lists
    int *vdeg = calloc(nv, sizeof(int));
    for(int e = 0; e < nm; e++) {
        for(int j = 0; j < nets[e].size; j++)
            vdeg[nets[e].vertices[j]]++;
    }
    vertices = malloc(nv * sizeof(Vertex));
    for(int v = 0; v < nv; v++) {
        vertices[v].degree = vdeg[v];
        vertices[v].nets = malloc(vertices[v].degree * sizeof(int));
    }
    free(vdeg);
    // fill incidence
    int *tmp = calloc(nv, sizeof(int));
    for(int e = 0; e < nm; e++) {
        for(int j = 0; j < nets[e].size; j++) {
            int v = nets[e].vertices[j];
            vertices[v].nets[tmp[v]++] = e;
        }
    }
    free(tmp);
    // assign vertex weights and compute total weight
    total_weight = 0;
    for (int v = 0; v < nv; v++) {
        vertices[v].weight = vweight[v];
        total_weight += vertices[v].weight;
    }
    free(vweight);
    // compute maximum weighted degree (sum of incident net weights) for bucket sizing
    max_degree = 0;
    for (int v = 0; v < nv; v++) {
        int wdeg = 0;
        for (int i = 0; i < vertices[v].degree; i++) {
            int e = vertices[v].nets[i];
            wdeg += nets[e].weight;
        }
        if (wdeg > max_degree)
            max_degree = wdeg;
    }

    // allocate bucket structures
    int bucket_len = 2 * max_degree + 1;
    bucket_heads = malloc(2 * sizeof(int*));
    for(int pid = 0; pid < 2; pid++)
        bucket_heads[pid] = malloc(bucket_len * sizeof(int));

    // neighbor helper
    neighbor_list = malloc(nv * sizeof(int));
    touched = calloc(nv, 1);

    // compute balance weight bounds for 2-way: [ (0.5 - eps)*W, (0.5 + eps)*W ]
    weight_lb = (int)ceil((0.5 - epsilon) * (double)total_weight);
    weight_ub = (int)floor((0.5 + epsilon) * (double)total_weight);
    int prev_cut = 0;  // best cut over runs
    // multi-start partition and FM refinement
    {
        int *best_part = malloc(nv * sizeof(int));
        if (!best_part) { perror("malloc"); return 1; }
        int best_cut = INT_MAX;
        for (int run = 0; run < num_runs; run++) {
            // initial partition for this run
            int split = -1;
            if (initial_method == 0) {
                // random initial partition with weight balance
                int *perm = malloc(nv * sizeof(int));
                if (!perm) { perror("malloc"); return 1; }
                for (int trial = 0; trial < 20; trial++) {
                    // generate random permutation
                    for (int i = 0; i < nv; i++) perm[i] = i;
                    for (int i = nv - 1; i > 0; i--) {
                        int j = rand() % (i + 1);
                        int t = perm[i]; perm[i] = perm[j]; perm[j] = t;
                    }
                    // find split point satisfying balance
                    long long cum = 0;
                    for (int i = 0; i < nv; i++) {
                        cum += vertices[perm[i]].weight;
                        if (cum >= weight_lb && cum <= weight_ub) { split = i; break; }
                        if (cum > weight_ub) break;
                    }
                    if (split >= 0) break;
                }
                if (split >= 0) {
                    for (int v = 0; v < nv; v++) vertices[v].part = 1;
                    for (int i = 0; i <= split; i++) vertices[perm[i]].part = 0;
                } else {
                    fprintf(stderr, "Warning: random partition failed to meet balance, falling back to greedy\n");
                }
                free(perm);
            }
            if (initial_method == 1 || split < 0) {
                // greedy partition
                int *order = malloc(nv * sizeof(int));
                if (!order) { perror("malloc"); return 1; }
                for (int i = 0; i < nv; i++) order[i] = i;
                qsort(order, nv, sizeof(int), cmp_degree);
                long long w0 = 0, w1 = 0;
                for (int i = 0; i < nv; i++) {
                    int v = order[i]; long long w = vertices[v].weight;
                    if ((w0 <= w1 && w0 + w <= weight_ub) || (w1 + w > weight_ub)) {
                        vertices[v].part = 0; w0 += w;
                    } else {
                        vertices[v].part = 1; w1 += w;
                    }
                }
                free(order);
            }
            // adjust balance bounds
            {
                long long w0 = 0, w1 = 0;
                for (int v = 0; v < nv; v++) {
                    if (vertices[v].part == 0) w0 += vertices[v].weight;
                    else                       w1 += vertices[v].weight;
                }
                if (w0 < weight_lb) {
                    int *plist = malloc(nv * sizeof(int));
                    int pc = 0;
                    for (int v = 0; v < nv; v++) if (vertices[v].part == 1) plist[pc++] = v;
                    qsort(plist, pc, sizeof(int), cmp_weight_asc);
                    for (int i = 0; i < pc && w0 < weight_lb; i++) {
                        int v = plist[i]; vertices[v].part = 0;
                        w0 += vertices[v].weight; w1 -= vertices[v].weight;
                    }
                    free(plist);
                } else if (w0 > weight_ub) {
                    int *plist = malloc(nv * sizeof(int));
                    int pc = 0;
                    for (int v = 0; v < nv; v++) if (vertices[v].part == 0) plist[pc++] = v;
                    qsort(plist, pc, sizeof(int), cmp_weight_asc);
                    for (int i = 0; i < pc && w0 > weight_ub; i++) {
                        int v = plist[i]; vertices[v].part = 1;
                        w0 -= vertices[v].weight; w1 += vertices[v].weight;
                    }
                    free(plist);
                }
            }
            // FM refinement for this run
            build_net_counts();
            int curr_cut = compute_cut();
            while (1) {
                fm_refine(); build_net_counts();
                int new_cut = compute_cut();
                if (new_cut < curr_cut) curr_cut = new_cut;
                else break;
            }
            // track best result
            if (curr_cut < best_cut) {
                best_cut = curr_cut;
                for (int v = 0; v < nv; v++) best_part[v] = vertices[v].part;
            }
        }
        // apply best partition
        for (int v = 0; v < nv; v++) vertices[v].part = best_part[v];
        free(best_part);
        prev_cut = best_cut;
        // mark end of compute phase
        clock_gettime(CLOCK_MONOTONIC, &ts_after_compute);
    }
    // output cut size
    printf("Final cut: %d\n", prev_cut);
    // write partition to file: <input>.part.2
    size_t fnlen = strlen(infile) + 8;
    char *outfname = malloc(fnlen);
    if (!outfname) {
        perror("malloc");
        return 1;
    }
    snprintf(outfname, fnlen, "%s.part.2", infile);
    FILE *out = fopen(outfname, "w");
    if (!out) {
        perror("fopen output");
        free(outfname);
        return 1;
    }
    for (int v = 0; v < nv; v++)
        fprintf(out, "%d\n", vertices[v].part);
    fclose(out);
    // mark end of write phase
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    // compute and print statistics
    int vert_count0 = 0, vert_count1 = 0;
    long long weight0 = 0, weight1 = 0;
    for (int v = 0; v < nv; v++) {
        if (vertices[v].part == 0) {
            vert_count0++;
            weight0 += vertices[v].weight;
        } else {
            vert_count1++;
            weight1 += vertices[v].weight;
        }
    }
    long long wdiff = llabs(weight0 - weight1);
    double wimb = (double)wdiff / (double)total_weight * 100.0;
    printf("Partition vertex counts: part0=%d, part1=%d\n", vert_count0, vert_count1);
    printf("Partition weights: part0=%lld, part1=%lld\n", weight0, weight1);
    // report actual epsilon based on achieved weight imbalance
    {
        // actual epsilon = max deviation from ideal (W/2) normalized by total weight
        double actual_eps = (double)wdiff / (2.0 * (double)total_weight);
        printf("Epsilon: %.6f\n", actual_eps);
    }
    // check balance constraints
    printf("Balance bounds (per part): [%d, %d]\n", weight_lb, weight_ub);
    printf("Partition 0 within bounds: %s\n",
           (weight0 >= weight_lb && weight0 <= weight_ub) ? "YES" : "NO");
    printf("Partition 1 within bounds: %s\n",
           (weight1 >= weight_lb && weight1 <= weight_ub) ? "YES" : "NO");
    free(outfname);
    // print timing information
    {
        double t_read = diff_ts(&ts_start, &ts_after_read);
        double t_compute = diff_ts(&ts_after_read, &ts_after_compute);
        double t_write = diff_ts(&ts_after_compute, &ts_end);
        double t_total = diff_ts(&ts_start, &ts_end);
        printf("Time (s): read=%.6f, compute=%.6f, write=%.6f, total=%.6f\n",
               t_read, t_compute, t_write, t_total);
    }
    return 0;
}