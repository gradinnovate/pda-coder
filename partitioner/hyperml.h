#ifndef HYPERML_H
#define HYPERML_H

#include <stddef.h>

typedef struct {
    int id;              // vertex ID (0-based)
    size_t weight;       // vertex weight
    int* edges;          // list of incident hyperedge IDs
    size_t edge_count;   // number of incident hyperedges
} Vertex;

typedef struct {
    int id;              // hyperedge ID (0-based)
    int* pins;           // list of vertex IDs (0-based)
    size_t pin_count;    // number of pins
    size_t weight;       // hyperedge weight
} Hyperedge;

typedef struct {
    int num_vertices;    // number of vertices
    int num_edges;       // number of hyperedges
    Vertex* vertices;    // array of vertices
    Hyperedge* edges;    // array of hyperedges
} Hypergraph;

// Read a hypergraph from a .hgr file. Returns NULL on failure.
Hypergraph* read_hgraph(const char* filename);
// Free hypergraph data structures.
void free_hgraph(Hypergraph* hg);
// Perform a greedy 2-way partition of the hypergraph. Returns an array of length num_vertices
// with values 0 or 1 indicating partition assignment. Caller must free().
int* greedy_partition(const Hypergraph* hg);
// Compute the cut size (total weight of hyperedges crossing the partition).
size_t compute_cut_size(const Hypergraph* hg, const int* part);
// Compute the maximum relative balance deviation of partitions (k=2).
// Deviation = max_i(|W_i - W/2|) / (W/2).
double compute_balance_deviation(const Hypergraph* hg, const int* part);
// Compute partition weights w0 and w1.
void compute_partition_weights(const Hypergraph* hg, const int* part,
                               size_t* w0, size_t* w1);
// Coarsen hypergraph using heavy-edge matching (hMETIS-style HEM).
// coarse_map should be an array of length src->num_vertices mapping each vertex to a coarse ID.
// Returns the new coarse hypergraph; caller is responsible for free_hgraph(new_hg) and free(coarse_map).
Hypergraph* coarsen_hgraph(const Hypergraph* src, int* coarse_map, double epsilon);
// Refine a 2-way partition on the hypergraph using a single pass of hypergraph FM.
// part is an array of length hg->num_vertices with values 0 or 1; modified in-place.
// epsilon is the allowed imbalance tolerance.
void refine_partition(const Hypergraph* hg, int* part, double epsilon);

#endif // HYPERML_H