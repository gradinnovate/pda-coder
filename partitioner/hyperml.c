#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "hyperml.h"
#include <limits.h>

// Read a hypergraph from a .hgr file
Hypergraph* read_hgraph(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        perror("Failed to open input file");
        return NULL;
    }
    int num_nets = 0, num_verts = 0;
    if (fscanf(fp, "%d %d", &num_nets, &num_verts) != 2) {
        fprintf(stderr, "Invalid header in .hgr file\n");
        fclose(fp);
        return NULL;
    }
    // consume end of line
    int c;
    while ((c = fgetc(fp)) != EOF && c != '\n');
    Hypergraph* hg = malloc(sizeof(Hypergraph));
    hg->num_vertices = num_verts;
    hg->num_edges = num_nets;
    hg->vertices = calloc((size_t)num_verts, sizeof(Vertex));
    hg->edges = calloc((size_t)num_nets, sizeof(Hyperedge));
    if (!hg->vertices || !hg->edges) {
        fprintf(stderr, "Memory allocation failure\n");
        fclose(fp);
        return NULL;
    }
    // Initialize vertices
    for (int i = 0; i < num_verts; ++i) {
        hg->vertices[i].id = i;
        hg->vertices[i].weight = 1;
        hg->vertices[i].edges = NULL;
        hg->vertices[i].edge_count = 0;
    }
    // Read hyperedges
    char* line = NULL;
    size_t len = 0;
    for (int e = 0; e < num_nets; ++e) {
        ssize_t read = getline(&line, &len, fp);
        if (read <= 0) {
            fprintf(stderr, "Unexpected end of file at edge %d\n", e);
            break;
        }
        // parse integers in line
        int* pins = NULL;
        size_t pin_count = 0;
        char* tok = strtok(line, " \t\n");
        while (tok) {
            int vid = atoi(tok);
            if (vid >= 1 && vid <= num_verts) {
                pins = realloc(pins, (pin_count + 1) * sizeof(int));
                pins[pin_count++] = vid - 1;
            }
            tok = strtok(NULL, " \t\n");
        }
        hg->edges[e].id = e;
        hg->edges[e].pins = pins;
        hg->edges[e].pin_count = pin_count;
        hg->edges[e].weight = 1;
        // count incidence
        for (size_t j = 0; j < pin_count; ++j) {
            int v = pins[j];
            hg->vertices[v].edge_count++;
        }
    }
    free(line);
    // Build vertex incidence lists
    for (int i = 0; i < num_verts; ++i) {
        size_t cnt = hg->vertices[i].edge_count;
        hg->vertices[i].edges = malloc(cnt * sizeof(int));
        hg->vertices[i].edge_count = 0;
    }
    for (int e = 0; e < num_nets; ++e) {
        for (size_t j = 0; j < hg->edges[e].pin_count; ++j) {
            int v = hg->edges[e].pins[j];
            size_t idx = hg->vertices[v].edge_count;
            hg->vertices[v].edges[idx] = e;
            hg->vertices[v].edge_count++;
        }
    }
    fclose(fp);
    return hg;
}

// Free hypergraph data structures
void free_hgraph(Hypergraph* hg) {
    if (!hg) return;
    for (int i = 0; i < hg->num_edges; ++i) {
        free(hg->edges[i].pins);
    }
    for (int i = 0; i < hg->num_vertices; ++i) {
        free(hg->vertices[i].edges);
    }
    free(hg->edges);
    free(hg->vertices);
    free(hg);
}

// Greedy 2-way partition (weights may be non-unit)
int* greedy_partition(const Hypergraph* hg) {
    int n = hg->num_vertices;
    int* part = malloc(n * sizeof(int));
    if (!part) return NULL;
    // Create array of vertex IDs
    int* order = malloc(n * sizeof(int));
    for (int i = 0; i < n; ++i) order[i] = i;
    // Sort by vertex weight descending
    // Currently all weights are equal; skip sorting
    size_t w0 = 0, w1 = 0;
    for (int idx = 0; idx < n; ++idx) {
        int v = order[idx];
        if (w0 <= w1) {
            part[v] = 0;
            w0 += hg->vertices[v].weight;
        } else {
            part[v] = 1;
            w1 += hg->vertices[v].weight;
        }
    }
    free(order);
    return part;
}

// Compute partition weights
void compute_partition_weights(const Hypergraph* hg, const int* part,
                               size_t* w0, size_t* w1) {
    *w0 = *w1 = 0;
    for (int i = 0; i < hg->num_vertices; ++i) {
        if (part[i] == 0) *w0 += hg->vertices[i].weight;
        else *w1 += hg->vertices[i].weight;
    }
}

// Compute cut size
size_t compute_cut_size(const Hypergraph* hg, const int* part) {
    size_t cut = 0;
    for (int e = 0; e < hg->num_edges; ++e) {
        int p0 = 0, p1 = 0;
        for (size_t j = 0; j < hg->edges[e].pin_count; ++j) {
            int v = hg->edges[e].pins[j];
            if (part[v] == 0) p0 = 1;
            else p1 = 1;
            if (p0 && p1) {
                cut += hg->edges[e].weight;
                break;
            }
        }
    }
    return cut;
}

// Compute maximum relative balance deviation
double compute_balance_deviation(const Hypergraph* hg, const int* part) {
    size_t w0, w1;
    compute_partition_weights(hg, part, &w0, &w1);
    double total = (double)(w0 + w1);
    if (total == 0) return 0.0;
    double ideal = total / 2.0;
    double d0 = fabs((double)w0 - ideal) / total;
    double d1 = fabs((double)w1 - ideal) / total;
    return (d0 > d1 ? d0 : d1);
}

// Compute gain for moving vertex v (FM formula)
static int compute_gain(const Hypergraph* hg, int v, const int* part, const size_t* eh_count0) {
    int b = part[v];
    int gain = 0;
    for (size_t i = 0; i < hg->vertices[v].edge_count; ++i) {
        int e = hg->vertices[v].edges[i];
        size_t pc = hg->edges[e].pin_count;
        size_t c0 = eh_count0[e];
        size_t c1 = pc - c0;
        if (b == 0) {
            if (c0 == 1) gain += (int)hg->edges[e].weight;
            if (c1 == 0) gain -= (int)hg->edges[e].weight;
        } else {
            if (c1 == 1) gain += (int)hg->edges[e].weight;
            if (c0 == 0) gain -= (int)hg->edges[e].weight;
        }
    }
    return gain;
}

// Coarsen hypergraph using Heavy-Edge Matching (HEM)
Hypergraph* coarsen_hgraph(const Hypergraph* src, int* coarse_map, double epsilon) {
    int n = src->num_vertices;
    int m = src->num_edges;
    int* matched = calloc(n, sizeof(int));
    int* last_update = calloc(n, sizeof(int));
    int* score = calloc(n, sizeof(int));
    int* neighbors = malloc(n * sizeof(int));
    if (!matched || !last_update || !score || !neighbors) {
        fprintf(stderr, "Memory allocation failure in coarsening\n");
        exit(EXIT_FAILURE);
    }
    // create random order of vertices for matching
    int* order = malloc(n * sizeof(int));
    if (!order) {
        fprintf(stderr, "Memory allocation failure in coarsening (order)\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; ++i) order[i] = i;
    for (int i = n - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        int tmp = order[i]; order[i] = order[j]; order[j] = tmp;
    }
    // total vertex weight
    size_t total_w = 0;
    for (int i = 0; i < n; ++i) total_w += src->vertices[i].weight;
    size_t max_cluster = (size_t)floor((0.5 + epsilon) * (double)total_w);
    int coarse_vid = 0;
    int timestamp = 1;
    for (int oi = 0; oi < n; ++oi) {
        int v = order[oi];
        if (matched[v]) continue;
        timestamp++;
        int nbr_count = 0;
        // score neighbors
        for (size_t ei = 0; ei < src->vertices[v].edge_count; ++ei) {
            int e = src->vertices[v].edges[ei];
            for (size_t pi = 0; pi < src->edges[e].pin_count; ++pi) {
                int u = src->edges[e].pins[pi];
                if (u == v || matched[u]) continue;
                if (last_update[u] != timestamp) {
                    last_update[u] = timestamp;
                    score[u] = 1;
                    neighbors[nbr_count++] = u;
                } else {
                    score[u]++;
                }
            }
        }
        // pick best neighbor
        int best_u = -1;
        int best_score = 0;
        for (int i = 0; i < nbr_count; ++i) {
            int u = neighbors[i];
            size_t wsum = src->vertices[v].weight + src->vertices[u].weight;
            if (wsum > max_cluster) continue;
            if (score[u] > best_score) {
                best_score = score[u];
                best_u = u;
            }
        }
        if (best_u >= 0) {
            matched[v] = matched[best_u] = 1;
            coarse_map[v] = coarse_map[best_u] = coarse_vid;
            coarse_vid++;
        } else {
            matched[v] = 1;
            coarse_map[v] = coarse_vid;
            coarse_vid++;
        }
    }
    int cn = coarse_vid;
    // allocate coarse graph
    Hypergraph* dst = malloc(sizeof(Hypergraph));
    dst->num_vertices = cn;
    dst->num_edges = m;
    dst->vertices = calloc(cn, sizeof(Vertex));
    dst->edges = calloc(m, sizeof(Hyperedge));
    if (!dst->vertices || !dst->edges) {
        fprintf(stderr, "Memory allocation failure in coarse graph\n");
        exit(EXIT_FAILURE);
    }
    // init coarse vertices
    for (int i = 0; i < cn; ++i) {
        dst->vertices[i].id = i;
        dst->vertices[i].weight = 0;
        dst->vertices[i].edges = NULL;
        dst->vertices[i].edge_count = 0;
    }
    // accumulate weights
    for (int v = 0; v < n; ++v) {
        int cv = coarse_map[v];
        dst->vertices[cv].weight += src->vertices[v].weight;
    }
    // map and build coarse edges, count incidences
    int* last_visit = calloc(cn, sizeof(int));
    int cv_ts = 1;
    for (int e = 0; e < m; ++e) {
        cv_ts++;
        int pin_max = src->edges[e].pin_count;
        int* pins = malloc(pin_max * sizeof(int));
        int pc = 0;
        for (size_t j = 0; j < src->edges[e].pin_count; ++j) {
            int v = src->edges[e].pins[j];
            int cv = coarse_map[v];
            if (last_visit[cv] != cv_ts) {
                last_visit[cv] = cv_ts;
                pins[pc++] = cv;
            }
        }
        dst->edges[e].id = e;
        dst->edges[e].pins = pins;
        dst->edges[e].pin_count = pc;
        dst->edges[e].weight = src->edges[e].weight;
        for (int j = 0; j < pc; ++j) {
            dst->vertices[pins[j]].edge_count++;
        }
    }
    free(last_visit);
    // build incidence lists
    for (int i = 0; i < cn; ++i) {
        dst->vertices[i].edges = malloc(dst->vertices[i].edge_count * sizeof(int));
        dst->vertices[i].edge_count = 0;
    }
    for (int e = 0; e < dst->num_edges; ++e) {
        for (int j = 0; j < dst->edges[e].pin_count; ++j) {
            int cv = dst->edges[e].pins[j];
            int idx = dst->vertices[cv].edge_count++;
            dst->vertices[cv].edges[idx] = e;
        }
    }
    free(matched);
    free(last_update);
    free(score);
    free(neighbors);
    free(order);
    return dst;
}

// Refine partition with a single-pass FM
void refine_partition(const Hypergraph* hg, int* part, double epsilon) {
    int n = hg->num_vertices;
    int m = hg->num_edges;
    // balance bounds
    size_t w0 = 0, w1 = 0;
    for (int i = 0; i < n; ++i) {
        if (part[i] == 0) w0 += hg->vertices[i].weight;
        else w1 += hg->vertices[i].weight;
    }
    size_t W = w0 + w1;
    double low_d = (0.5 - epsilon) * (double)W;
    double high_d = (0.5 + epsilon) * (double)W;
    size_t low = (size_t)ceil(low_d);
    size_t high = (size_t)floor(high_d);
    // initial eh_count0
    size_t* eh_count0 = calloc(m, sizeof(size_t));
    for (int e = 0; e < m; ++e) {
        size_t c0 = 0;
        for (int j = 0; j < hg->edges[e].pin_count; ++j) {
            if (part[hg->edges[e].pins[j]] == 0) c0++;
        }
        eh_count0[e] = c0;
    }
    // locked and gain arrays
    char* locked = calloc(n, sizeof(char));
    int* D = malloc(n * sizeof(int));
    for (int i = 0; i < n; ++i) {
        D[i] = compute_gain(hg, i, part, eh_count0);
    }
    // FM loop
    while (1) {
        int best_v = -1;
        int best_D = INT_MIN;
        for (int v = 0; v < n; ++v) {
            if (locked[v]) continue;
            int d = D[v];
            if (d <= 0) continue;
            // check balance
            size_t nw0 = w0, nw1 = w1;
            if (part[v] == 0) { nw0 -= hg->vertices[v].weight; nw1 += hg->vertices[v].weight; }
            else { nw1 -= hg->vertices[v].weight; nw0 += hg->vertices[v].weight; }
            if (nw0 < low || nw0 > high) continue;
            if (nw1 < low || nw1 > high) continue;
            if (d > best_D) { best_D = d; best_v = v; }
        }
        if (best_v < 0) break;
        // move
        int v = best_v;
        int from = part[v];
        if (from == 0) { w0 -= hg->vertices[v].weight; w1 += hg->vertices[v].weight; }
        else { w1 -= hg->vertices[v].weight; w0 += hg->vertices[v].weight; }
        part[v] = 1 - from;
        locked[v] = 1;
        // update eh_count0 and D
        for (size_t ei = 0; ei < hg->vertices[v].edge_count; ++ei) {
            int e = hg->vertices[v].edges[ei];
            // update count
            if (from == 0) eh_count0[e]--;
            else eh_count0[e]++;
        }
        // recompute D for unlocked vertices
        for (int u = 0; u < n; ++u) {
            if (!locked[u]) D[u] = compute_gain(hg, u, part, eh_count0);
        }
    }
    // Balancing phase: ensure strict balance constraint
    // Move from side 0 to 1 if too heavy
    while (w0 > high) {
        int best_v = -1;
        int best_gain = INT_MIN;
        for (int v = 0; v < n; ++v) {
            if (part[v] != 0) continue;
            size_t wv = hg->vertices[v].weight;
            size_t new_w0 = w0 - wv;
            if (new_w0 < low || new_w0 > high) continue;
            int g = compute_gain(hg, v, part, eh_count0);
            if (g > best_gain) {
                best_gain = g;
                best_v = v;
            }
        }
        if (best_v < 0) break;
        for (size_t ei = 0; ei < hg->vertices[best_v].edge_count; ++ei) {
            int e = hg->vertices[best_v].edges[ei];
            eh_count0[e]--;
        }
        part[best_v] = 1;
        w0 -= hg->vertices[best_v].weight;
        w1 += hg->vertices[best_v].weight;
    }
    // Move from side 1 to 0 if too light
    while (w0 < low) {
        int best_v = -1;
        int best_gain = INT_MIN;
        for (int v = 0; v < n; ++v) {
            if (part[v] != 1) continue;
            size_t wv = hg->vertices[v].weight;
            size_t new_w0 = w0 + wv;
            if (new_w0 < low || new_w0 > high) continue;
            int g = compute_gain(hg, v, part, eh_count0);
            if (g > best_gain) {
                best_gain = g;
                best_v = v;
            }
        }
        if (best_v < 0) break;
        for (size_t ei = 0; ei < hg->vertices[best_v].edge_count; ++ei) {
            int e = hg->vertices[best_v].edges[ei];
            eh_count0[e]++;
        }
        part[best_v] = 0;
        w1 -= hg->vertices[best_v].weight;
        w0 += hg->vertices[best_v].weight;
    }
    free(eh_count0);
    free(locked);
    free(D);
}

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "Usage: %s <input.hgr> [epsilon]\n", argv[0]);
        return EXIT_FAILURE;
    }
    clock_t start = clock();
    // seed RNG for randomized coarsening order
    srand((unsigned)time(NULL));
    // read original hypergraph
    Hypergraph* orig = read_hgraph(argv[1]);
    if (!orig) return EXIT_FAILURE;
    // balance tolerance epsilon (optional argument)
    double epsilon = 0.03;
    if (argc == 3) {
        char* endptr = NULL;
        epsilon = strtod(argv[2], &endptr);
        if (!endptr || *endptr != '\0' || epsilon < 0.0 || epsilon > 0.5) {
            fprintf(stderr, "Invalid epsilon value '%s'\n", argv[2]);
            free_hgraph(orig);
            return EXIT_FAILURE;
        }
    }
    int threshold = 100;  // coarsest size threshold
    // multilevel hierarchy
    Hypergraph** levels = NULL;
    int** maps = NULL;
    size_t nlevels = 0;
    levels = malloc(sizeof(Hypergraph*));
    levels[0] = orig;
    nlevels = 1;
    // coarsening
    while (levels[nlevels-1]->num_vertices > threshold) {
        Hypergraph* curr = levels[nlevels-1];
        int nv = curr->num_vertices;
        int* cmap = malloc(nv * sizeof(int));
        Hypergraph* coarse = coarsen_hgraph(curr, cmap, epsilon);
        maps = realloc(maps, nlevels * sizeof(int*));
        maps[nlevels-1] = cmap;
        levels = realloc(levels, (nlevels+1) * sizeof(Hypergraph*));
        levels[nlevels] = coarse;
        nlevels++;
    }
    // initial partition on coarsest graph
    Hypergraph* coarse = levels[nlevels-1];
    int* part = greedy_partition(coarse);
    if (!part) return EXIT_FAILURE;
    refine_partition(coarse, part, epsilon);
    // uncoarsen and refine
    for (int l = (int)nlevels - 1; l > 0; --l) {
        Hypergraph* hg = levels[l-1];
        int* cmap = maps[l-1];
        int nv = hg->num_vertices;
        int* new_part = malloc(nv * sizeof(int));
        for (int i = 0; i < nv; ++i) new_part[i] = part[cmap[i]];
        free(part);
        part = new_part;
        refine_partition(hg, part, epsilon);
    }
    // compute final metrics
    size_t w0 = 0, w1 = 0;
    compute_partition_weights(orig, part, &w0, &w1);
    size_t cut = compute_cut_size(orig, part);
    double dev = compute_balance_deviation(orig, part);
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("CutSize %zu\n", cut);
    printf("Partition Sizes: %zu, %zu\n", w0, w1);
    printf("Balance Deviation: %.6f\n", dev);
    printf("Total Execution Time: %.6f\n", elapsed);
    // cleanup
    free(part);
    for (size_t i = 0; i < nlevels - 1; ++i) free(maps[i]);
    free(maps);
    for (size_t i = 0; i < nlevels; ++i) free_hgraph(levels[i]);
    free(levels);
    return EXIT_SUCCESS;
}