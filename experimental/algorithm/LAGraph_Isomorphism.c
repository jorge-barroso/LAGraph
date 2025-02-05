//
// Created by jorge on 01/25/25.
//


#include <GraphBLAS.h>
#include <stdbool.h>

// Problem types
typedef enum {
    ISO,    // Graph isomorphism
    IND,    // Induced subgraph isomorphism
    SUB     // Subgraph isomorphism
} ProblemType;

// Node labels
typedef struct {
    int* labels;        // Node labels
    GrB_Index size;     // Number of nodes
    int num_labels;     // Number of unique labels
    int* label_freqs;   // Frequency of each label
} NodeLabels;

// BFS level information for ordering
typedef struct {
    GrB_Index* nodes;   // Nodes at this level
    GrB_Index size;     // Number of nodes
} BFSLevel;

// VF2++ state
typedef struct {
    GrB_Index* core_1;  // Mapping from G1 to G2
    GrB_Index* core_2;  // Mapping from G2 to G1
    GrB_Index size1;    // Size of G1
    GrB_Index size2;    // Size of G2
    GrB_Index depth;    // Current depth in search

    // Terminal sets
    bool* in_1;         // In-terminal set for G1
    bool* in_2;         // In-terminal set for G2
    bool* out_1;        // Out-terminal set for G1
    bool* out_2;        // Out-terminal set for G2

    // Pre-computed ordering
    GrB_Index* order;   // Node ordering for G1
} VF2State;

enum Level { FASTEST, FAST, SLOW, SLOWEST };

//#define LG_FREE_WORK                \
//{                                   \
//    GrB_free (&degrees1) ;              \
//    GrB_free (&degrees2) ;              \
//}
#define LG_FREE_WORK {}

#include "LG_internal.h"


static bool vf2pp_match_internal(const GrB_Matrix G1, const GrB_Matrix G2,
                                NodeLabels* labels1, NodeLabels* labels2,
                                VF2State* state, const ProblemType type);
static void compute_node_order(const GrB_Matrix G1, const NodeLabels* labels1, const VF2State* state);
VF2State* init_state(GrB_Matrix G1, GrB_Matrix G2);
void free_state(VF2State* state);
int compare_sizes(bool *is_isomorphic,
                 const LAGraph_Graph G1,
                 const LAGraph_Graph G2,
                 GrB_Index *size1,
                 GrB_Index *size2,
                 char *msg);
int compute_degrees(bool *is_isomorphic,
                    const LAGraph_Graph G1,
                    const LAGraph_Graph G2,
                    GrB_Vector degrees1,
                    GrB_Vector degrees2,
                    GrB_Index *nvals1,
                    GrB_Index *nvals2,
                    char *msg);
int allocate_and_extract_tuples(int64_t *values1,
                                int64_t *values2,
                                GrB_Index nvals1,
                                GrB_Index nvals2,
                                const GrB_Vector degrees1,
                                const GrB_Vector degrees2,
                                char *msg);
int compute_triangles(bool *is_isomorphic,
                    const LAGraph_Graph G,
                    GrB_Vector degrees,
                    GrB_Index *nvals,
                    char *msg);
int compare_nvals(bool *is_isomorphic,
                  GrB_Vector triangles1,
                  GrB_Vector triangles2,
                  GrB_Index *triangles_nvals1,
                  GrB_Index *triangles_nvals2,
                  GrB_Index *degrees_nvals1,
                  GrB_Index *degrees_nvals2,
                  GrB_Index size1,
                  GrB_Index size2,
                  char *msg);
int do_fastest(bool *is_isomorphic,
                  const LAGraph_Graph G1,
                  const LAGraph_Graph G2,
                  char *msg);
int do_fast(bool *is_isomorphic,
                  const LAGraph_Graph G1,
                  const LAGraph_Graph G2,
                  char *msg);
int do_slow(bool *is_isomorphic,
                  const LAGraph_Graph G1,
                  const LAGraph_Graph G2,
                  char *msg);
int do_slowest(bool *is_isomorphic,
                  const LAGraph_Graph G1,
                  const LAGraph_Graph G2,
                  char *msg);
int LAGraph_IsIsomorphic(
    bool *is_isomorphic,
    LAGraph_Graph G1,
    LAGraph_Graph G2,
    enum Level level,
    char *msg) {
    LG_CLEAR_MSG;

    (*is_isomorphic) = true;
    if (G1->kind != G2->kind) {
        return (*is_isomorphic) = false;
        return (GrB_SUCCESS);
    }

    switch (level) {
        case FASTEST:    do_fastest(is_isomorphic, G1, G2, msg); break;
        case FAST:    do_fast(is_isomorphic, G1, G2, msg); break;
        case SLOW:    do_slow(is_isomorphic, G1, G2, msg); break;
        case SLOWEST:    do_slowest(is_isomorphic, G1, G2, msg); break;
    }

//    GrB_Vector triangles1 = NULL, triangles2 = NULL;
//    GrB_Matrix L1 = NULL, U1 = NULL;
//    GRB_TRY(GrB_Matrix_new(&L1, GrB_INT64, n, n));
//    GRB_TRY(GrB_Matrix_new(&U1, GrB_INT64, n, n));
//    GRB_TRY(GrB_select(L1, NULL, NULL, GrB_TRIL, G1->A, 0, NULL));
//    GRB_TRY(GrB_select(U1, NULL, NULL, GrB_TRIU, G1->A, 0, NULL));
//
//    GrB_Index tri_nvals1, tri_nvals2;
//    GRB_TRY(GrB_Vector_nvals(&tri_nvals1, triangles1));
//    GRB_TRY(GrB_Vector_nvals(&tri_nvals2, triangles2));
//    if (tri_nvals1 != tri_nvals2) {
//        (*is_isomorphic) = false;
//        return GrB_SUCCESS;
//    }
//
//    // Initialize state
//    VF2State* state = init_state(G1, G2);
//
//    // Compute node ordering
//    compute_node_order(G1, labels1, state);
//
//    // Find mapping
//    bool success = vf2pp_match_internal(G1, G2, labels1, labels2, state, type);
//
//    if (success && mapping != NULL) {
//        *mapping = (GrB_Index*)malloc(n1 * sizeof(GrB_Index));
//        memcpy(*mapping, state->core_1, n1 * sizeof(GrB_Index));
//    }
//
//    free_state(state);
//    return success;

    return (GrB_SUCCESS);
}

// BFS-based node ordering computation
static void compute_node_order(const GrB_Matrix G1, const NodeLabels* labels1, const VF2State* state) {
    const GrB_Index n = state->size1;
    bool* visited = calloc(n, sizeof(bool));
    GrB_Index order_idx = 0;

    // Process each component
    for (GrB_Index start = 0; start < n; start++) {
        if (visited[start]) continue;

        // Find start node (highest degree, rarest label)
        GrB_Index best_start = start;
        GrB_Index max_deg = 0;
        for (GrB_Index i = start; i < n; i++) {
            if (visited[i]) continue;

            // Get degree
            GrB_Vector v;
            GrB_Vector_new(&v, GrB_BOOL, n);
            GrB_Col_extract(v, NULL, NULL, G1, GrB_ALL, n, i, NULL);
            GrB_Index deg;
            GrB_reduce(&deg, NULL, GrB_PLUS_MONOID_INT64, v, NULL);
            GrB_free(&v);

            if (deg > max_deg ||
                (deg == max_deg &&
                 labels1->label_freqs[labels1->labels[i]] <
                 labels1->label_freqs[labels1->labels[best_start]])) {
                max_deg = deg;
                best_start = i;
            }
        }

        // Do BFS from best start node
        GrB_Vector frontier;
        GrB_Vector_new(&frontier, GrB_BOOL, n);
        GrB_Vector_setElement(frontier, true, best_start);

        while (true) {
            // Get nodes at current level
            GrB_Index* level_nodes = NULL;
            GrB_Index level_size = 0;

            for (GrB_Index i = 0; i < n; i++) {
                bool in_frontier;
                GrB_Vector_extractElement(&in_frontier, frontier, i);
                if (in_frontier) {
                    level_nodes = realloc(level_nodes,
                                        (level_size + 1) * sizeof(GrB_Index));
                    level_nodes[level_size++] = i;
                }
            }

            if (level_size == 0) {
                free(level_nodes);
                break;
            }

            // Sort level nodes by (conn, deg, label_freq)
            // Using insertion sort for simplicity
            for (GrB_Index i = 1; i < level_size; i++) {
                const GrB_Index key = level_nodes[i];
                GrB_Index j = i - 1;

                while (j >= 0) {
                    const GrB_Index node1 = level_nodes[j];
                    const GrB_Index node2 = key;

                    // Get connection counts
                    GrB_Index conn1 = 0, conn2 = 0;
                    for (GrB_Index k = 0; k < order_idx; k++) {
                        bool edge1, edge2;
                        GrB_Matrix_extractElement_BOOL(&edge1, G1, node1,
                                                     state->order[k]);
                        GrB_Matrix_extractElement_BOOL(&edge2, G1, node2,
                                                     state->order[k]);
                        if (edge1) conn1++;
                        if (edge2) conn2++;
                    }

                    // Get degrees
                    GrB_Vector v;
                    GrB_Vector_new(&v, GrB_BOOL, n);
                    GrB_Index deg1, deg2;

                    GrB_Col_extract(v, NULL, NULL, G1, GrB_ALL, n, node1, NULL);
                    GrB_reduce(&deg1, NULL, GrB_PLUS_MONOID_INT64, v, NULL);

                    GrB_Col_extract(v, NULL, NULL, G1, GrB_ALL, n, node2, NULL);
                    GrB_reduce(&deg2, NULL, GrB_PLUS_MONOID_INT64, v, NULL);

                    GrB_free(&v);

                    // Compare by criteria
                    if (conn1 > conn2 ||
                        (conn1 == conn2 && deg1 > deg2) ||
                        (conn1 == conn2 && deg1 == deg2 &&
                         labels1->label_freqs[labels1->labels[node1]] <
                         labels1->label_freqs[labels1->labels[node2]])) {
                        break;
                    }

                    level_nodes[j + 1] = level_nodes[j];
                    j--;
                }
                level_nodes[j + 1] = key;
            }

            // Add level nodes to order
            for (GrB_Index i = 0; i < level_size; i++) {
                state->order[order_idx++] = level_nodes[i];
                visited[level_nodes[i]] = true;
            }

            // Compute next frontier
            GrB_mxv(frontier, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL,
                   G1, frontier, NULL);

            // Remove visited nodes
            for (GrB_Index i = 0; i < n; i++) {
                if (visited[i]) {
                    GrB_Vector_setElement(frontier, false, i);
                }
            }

            free(level_nodes);
        }

        GrB_free(&frontier);
    }

    free(visited);
}

// Terminal set computation
static void update_terminal_sets(const GrB_Matrix G1, const GrB_Matrix G2,
                               const VF2State* state,
                               const GrB_Index node1, const GrB_Index node2) {
    // Update in-terminal sets
    GrB_Vector v1, v2;
    GrB_Vector_new(&v1, GrB_BOOL, state->size1);
    GrB_Vector_new(&v2, GrB_BOOL, state->size2);

    // Get neighbors
    GrB_Col_extract(v1, NULL, NULL, G1, GrB_ALL, state->size1, node1, NULL);
    GrB_Col_extract(v2, NULL, NULL, G2, GrB_ALL, state->size2, node2, NULL);

    for (GrB_Index i = 0; i < state->size1; i++) {
        bool is_neighbor;
        GrB_Vector_extractElement_BOOL(&is_neighbor, v1, i);
        if (is_neighbor && !state->in_1[i] && state->core_1[i] == state->size2) {
            state->in_1[i] = true;
        }
    }

    for (GrB_Index i = 0; i < state->size2; i++) {
        bool is_neighbor;
        GrB_Vector_extractElement_BOOL(&is_neighbor, v2, i);
        if (is_neighbor && !state->in_2[i] && state->core_2[i] == state->size1) {
            state->in_2[i] = true;
        }
    }

    // Update out-terminal sets similarly
    GrB_free(&v1);
    GrB_free(&v2);
}

// Feasibility checking functions
static bool check_syntactic_feasibility(const GrB_Matrix G1, const GrB_Matrix G2,
                                      const VF2State* state,
                                      const GrB_Index node1, const GrB_Index node2,
                                      const ProblemType type) {
    // Check edges to/from already matched nodes
    for (GrB_Index i = 0; i < state->depth; i++) {
        const GrB_Index m1 = state->order[i];
        const GrB_Index m2 = state->core_1[m1];

        bool edge1_in, edge1_out, edge2_in, edge2_out;

        GrB_Matrix_extractElement_BOOL(&edge1_in, G1, node1, m1);
        GrB_Matrix_extractElement_BOOL(&edge1_out, G1, m1, node1);
        GrB_Matrix_extractElement_BOOL(&edge2_in, G2, node2, m2);
        GrB_Matrix_extractElement_BOOL(&edge2_out, G2, m2, node2);

        if (type == SUB) {
            if ((edge1_in && !edge2_in) || (edge1_out && !edge2_out)) {
                return false;
            }
        } else { // ISO or IND
            if (edge1_in != edge2_in || edge1_out != edge2_out) {
                return false;
            }
        }
    }

    return true;
}

static bool check_semantic_feasibility(const NodeLabels* labels1, const NodeLabels* labels2,
                                     const GrB_Index node1, const GrB_Index node2) {
    return labels1->labels[node1] == labels2->labels[node2];
}

static bool check_label_compatibility(const GrB_Matrix G1, const GrB_Matrix G2,
                                    const VF2State* state,
                                    const NodeLabels* labels1, const NodeLabels* labels2,
                                    const GrB_Index node1, const GrB_Index node2,
                                    const ProblemType type) {
    // Count label frequencies in neighborhoods
    int* label_counts1 = (int*)calloc(labels1->num_labels, sizeof(int));
    int* label_counts2 = (int*)calloc(labels2->num_labels, sizeof(int));

    // Count in terminal sets
    for (GrB_Index i = 0; i < state->size1; i++) {
        if (state->in_1[i]) {
            bool edge;
            GrB_Matrix_extractElement_BOOL(&edge, G1, node1, i);
            if (edge) {
                label_counts1[labels1->labels[i]]++;
            }
        }
    }

    for (GrB_Index i = 0; i < state->size2; i++) {
        if (state->in_2[i]) {
            bool edge;
            GrB_Matrix_extractElement_BOOL(&edge, G2, node2, i);
            if (edge) {
                label_counts2[labels2->labels[i]]++;
            }
        }
    }

    // Check counts based on problem type
    for (int l = 0; l < labels1->num_labels; l++) {
        if (type == SUB) {
            if (label_counts1[l] > label_counts2[l]) {
                free(label_counts1);
                free(label_counts2);
                return false;
            }
        } else { // ISO or IND
            if (label_counts1[l] != label_counts2[l]) {
                free(label_counts1);
                free(label_counts2);
                return false;
            }
        }
    }

    free(label_counts1);
    free(label_counts2);
    return true;
}

// Main matching function
static bool vf2pp_match_internal(const GrB_Matrix G1, const GrB_Matrix G2,
                                NodeLabels* labels1, NodeLabels* labels2,
                                VF2State* state, const ProblemType type) {
    if (state->depth == state->size1) {
        return true;
    }

    const GrB_Index node1 = state->order[state->depth];

    // Get candidates with matching label
    for (GrB_Index node2 = 0; node2 < state->size2; node2++) {
        // Skip if already matched
        if (state->core_2[node2] != state->size1) {
            continue;
        }

        // Check semantic feasibility (label matching)
        if (!check_semantic_feasibility(labels1, labels2, node1, node2)) {
            continue;
        }

        // Check syntactic feasibility
        if (!check_syntactic_feasibility(G1, G2, state, node1, node2, type)) {
            continue;
        }

        // Check label compatibility (cutting rules)
        if (!check_label_compatibility(G1, G2, state, labels1, labels2,
                                     node1, node2, type)) {
            continue;
        }

        // Add to mapping
        state->core_1[node1] = node2;
        state->core_2[node2] = node1;
        state->depth++;

        // Update terminal sets
        update_terminal_sets(G1, G2, state, node1, node2);

        // Recurse
        if (vf2pp_match_internal(G1, G2, labels1, labels2, state, type)) {
            return true;
        }

        // Backtrack
        state->core_1[node1] = state->size2;
        state->core_2[node2] = state->size1;
        state->depth--;
    }

    return false;
}

// API implementation
NodeLabels* init_node_labels(int* labels, GrB_Index size) {
    NodeLabels* nl = (NodeLabels*)malloc(sizeof(NodeLabels));
    nl->labels = (int*)malloc(size * sizeof(int));
    memcpy(nl->labels, labels, size * sizeof(int));
    nl->size = size;

    // Count unique labels and frequencies
    nl->num_labels = 0;
    nl->label_freqs = (int*)calloc(size, sizeof(int));

    for (GrB_Index i = 0; i < size; i++) {
        nl->label_freqs[labels[i]]++;
        if (nl->label_freqs[labels[i]] == 1) {
            nl->num_labels++;
        }
    }

    return nl;
}

void free_node_labels(NodeLabels* labels) {
    free(labels->labels);
    free(labels->label_freqs);
    free(labels);
}

VF2State* init_state(GrB_Matrix G1, GrB_Matrix G2) {
    GrB_Index n1, n2;
    GrB_Matrix_nrows(&n1, G1);
    GrB_Matrix_nrows(&n2, G2);

    VF2State* state = (VF2State*)malloc(sizeof(VF2State));

    state->core_1 = (GrB_Index*)malloc(n1 * sizeof(GrB_Index));
    state->core_2 = (GrB_Index*)malloc(n2 * sizeof(GrB_Index));
    state->order = (GrB_Index*)malloc(n1 * sizeof(GrB_Index));

    state->in_1 = (bool*)calloc(n1, sizeof(bool));
    state->in_2 = (bool*)calloc(n2, sizeof(bool));
    state->out_1 = (bool*)calloc(n1, sizeof(bool));
    state->out_2 = (bool*)calloc(n2, sizeof(bool));

    state->size1 = n1;
    state->size2 = n2;
    state->depth = 0;

    // Initialize mappings to invalid values
    for (GrB_Index i = 0; i < n1; i++) {
        state->core_1[i] = n2;
    }
    for (GrB_Index i = 0; i < n2; i++) {
        state->core_2[i] = n1;
    }

    return state;
}

void free_state(VF2State* state) {
    free(state->core_1);
    free(state->core_2);
    free(state->order);
    free(state->in_1);
    free(state->in_2);
    free(state->out_1);
    free(state->out_2);
    free(state);
}

int compare_int64(const void *a, const void *b) {
    int64_t va = *(const int64_t*)a;
    int64_t vb = *(const int64_t*)b;
    // 1, 0 or -1 depending on whether va is greater, equal or less than vb
    return (va > vb) - (va < vb);
}

int compare_sizes(bool *is_isomorphic,
                 const LAGraph_Graph G1,
                 const LAGraph_Graph G2,
                 GrB_Index *size1,
                 GrB_Index *size2,
                 char *msg) {
    // We only need rows because Ajd Matrices are squared by definition
    GRB_TRY(GrB_Matrix_nrows(size1, G1->A));
    GRB_TRY(GrB_Matrix_nrows(size2, G2->A));
    if (size1 != size2) {
        (*is_isomorphic) = false;
    }
    return GrB_SUCCESS;
}

int compute_degrees(bool *is_isomorphic,
                    const LAGraph_Graph G1,
                    const LAGraph_Graph G2,
                    GrB_Vector degrees1,
                    GrB_Vector degrees2,
                    GrB_Index *nvals1,
                    GrB_Index *nvals2,
                    char *msg) {
    // Calculate degrees nvals
    LAGraph_Cached_OutDegree(G1, msg);
    LAGraph_Cached_OutDegree(G2, msg);
    if (G1->kind != LAGraph_ADJACENCY_DIRECTED) {
        LAGraph_Cached_InDegree(G1, msg);
        LAGraph_Cached_InDegree(G2, msg);
        GRB_TRY(GrB_eWiseAdd(degrees1, NULL, NULL, GrB_PLUS_INT64, G1->in_degree, G1->out_degree, NULL));
        GRB_TRY(GrB_eWiseAdd(degrees2, NULL, NULL, GrB_PLUS_INT64, G2->in_degree, G2->out_degree, NULL));
    } else {
        degrees1 = G1->in_degree;
        degrees2 = G2->in_degree;
    }

    GrB_Vector_nvals(nvals1, degrees1);
    GrB_Vector_nvals(nvals2, degrees2);
    if(nvals1 != nvals2) {
        (*is_isomorphic) = false;
    }

    return (GrB_SUCCESS);
}

int allocate_and_extract_tuples(int64_t *values1,
                                int64_t *values2,
                                GrB_Index nvals1,
                                GrB_Index nvals2,
                                const GrB_Vector degrees1,
                                const GrB_Vector degrees2,
                                char *msg) {
    if (values1 == NULL || values2 == NULL) {
        if (values1) free(values1);
        if (values2) free(values2);
        return GrB_OUT_OF_MEMORY;
    }
    GrB_Index actual_nvals;
    GRB_TRY(GrB_Vector_extractTuples(NULL, values1, &actual_nvals, degrees1));
    GRB_TRY(GrB_Vector_extractTuples(NULL, values2, &actual_nvals, degrees2));
    return GrB_SUCCESS;
}

int compute_triangles(bool *is_isomorphic,
                    const LAGraph_Graph G,
                    GrB_Vector triangles,
                    GrB_Index *nvals,
                    char *msg) {
    GrB_Matrix L = NULL, U = NULL;
    GrB_Matrix C = NULL, D = NULL;
    GrB_Vector temp = NULL;
    GrB_Index n;

    // Get matrix size
    GRB_TRY(GrB_Matrix_nrows(&n, G->A));

    // Extract L and U triangular parts
    GRB_TRY(GrB_Matrix_new(&L, GrB_INT64, n, n));
    GRB_TRY(GrB_Matrix_new(&U, GrB_INT64, n, n));
    GRB_TRY(GrB_select(L, NULL, NULL, GrB_TRIL, G->A, 0, NULL));
    GRB_TRY(GrB_select(U, NULL, NULL, GrB_TRIU, G->A, 0, NULL));

    // Create output vector and temporary matrices
    GRB_TRY(GrB_Vector_new(&triangles, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&temp, GrB_INT64, n));
    GRB_TRY(GrB_Matrix_new(&C, GrB_INT64, n, n));
    GRB_TRY(GrB_Matrix_new(&D, GrB_INT64, n, n));

    // Use plus_pair semiring like in Python version
    GrB_Semiring semiring = GxB_PLUS_PAIR_INT64;

    // C = (L * L') .* L
    GRB_TRY(GrB_mxm(C, L, NULL, semiring, L, L, GrB_DESC_ST0));

    // D = (U * L') .* U
    GRB_TRY(GrB_mxm(D, U, NULL, semiring, U, L, GrB_DESC_ST0));

    // Add up the three components like in Python:
    // C.reduce_rowwise + C.reduce_columnwise + D.reduce_rowwise
     GRB_TRY(GrB_reduce(triangles, NULL, GrB_PLUS_MONOID_INT64, C, GrB_DESC_T0));  // rowwise C
     GRB_TRY(GrB_reduce(temp, NULL, GrB_PLUS_MONOID_INT64, C, GrB_DESC_T1));        // columnwise C
     GRB_TRY(GrB_eWiseAdd(triangles, NULL, NULL, GrB_PLUS_INT64, triangles, temp, NULL));
     GRB_TRY(GrB_reduce(temp, NULL, GrB_PLUS_MONOID_INT64, D, GrB_DESC_T0));        // rowwise D
     GRB_TRY(GrB_eWiseAdd(triangles, NULL, NULL, GrB_PLUS_INT64, triangles, temp, NULL));

    // Free temporary objects
    GrB_Matrix_free(&L);
    GrB_Matrix_free(&U);
    GrB_Matrix_free(&C);
    GrB_Matrix_free(&D);
    GrB_Vector_free(&temp);

    return GrB_SUCCESS;
    return (GrB_SUCCESS);
}

int compare_nvals(bool *is_isomorphic,
                  GrB_Vector triangles1,
                  GrB_Vector triangles2,
                  GrB_Index *triangles_nvals1,
                  GrB_Index *triangles_nvals2,
                  GrB_Index *degrees_nvals1,
                  GrB_Index *degrees_nvals2,
                  GrB_Index size1,
                  GrB_Index size2,
                  char *msg) {
    // First check if triangle counts match
    if (triangles_nvals1 != triangles_nvals2) {
        *is_isomorphic = false;
        return GrB_SUCCESS;
    }

    // Create mask for entries present in degrees but not in triangles
    GrB_Vector mask1, mask2;
    GrB_Vector_new(&mask1, GrB_BOOL, size1);  // size1 is size of G1
    GrB_Vector_new(&mask2, GrB_BOOL, size2);  // size2 is size of G2

    // Get structure of triangles vectors
    GrB_Vector struct1, struct2;
    GrB_Vector_new(&struct1, GrB_BOOL, size1);
    GrB_Vector_new(&struct2, GrB_BOOL, size2);
    GrB_Vector_apply(struct1, NULL, NULL, GrB_IDENTITY_BOOL, triangles1, NULL);
    GrB_Vector_apply(struct2, NULL, NULL, GrB_IDENTITY_BOOL, triangles2, NULL);

    // Complement to get entries in degrees but not in triangles
     GrB_Vector_apply(mask1, NULL, NULL, GrB_LNOT, struct1, NULL);
     GrB_Vector_apply(mask2, NULL, NULL, GrB_LNOT, struct2, NULL);

    // Create vectors to hold the filled triangle counts
    GrB_Vector filled_triangles1, filled_triangles2;
    GrB_Vector_dup(&filled_triangles1, triangles1);
    GrB_Vector_dup(&filled_triangles2, triangles2);

    // Fill missing entries with 0 where degrees exist but triangles don't
    GrB_Vector zeros1, zeros2;
    GrB_Vector_new(&zeros1, GrB_INT64, size1);
    GrB_Vector_new(&zeros2, GrB_INT64, size2);

    // Use mask to set zeros where needed
    GrB_Vector_assign(filled_triangles1, mask1, NULL, 0, GrB_ALL, size1, NULL);
    GrB_Vector_assign(filled_triangles2, mask2, NULL, 0, GrB_ALL, size2, NULL);
}

int sort_and_compare(bool *is_isomorphic,
                     int64_t *values1,
                     int64_t *values2,
                     GrB_Index nvals1,
                     GrB_Index nvals2) {
    qsort(values1, nvals1, sizeof(int64_t), compare_int64);
    qsort(values2, nvals2, sizeof(int64_t), compare_int64);
    for (GrB_Index i = 0; i < nvals1; i++) {
        if (values1[i] != values2[i]) {
            (*is_isomorphic) = false;
            break;
        }
    }
    return GrB_SUCCESS;
}

int do_fastest(bool *is_isomorphic,
                 const LAGraph_Graph G1,
                 const LAGraph_Graph G2,
                 char *msg) {
    GrB_Index size1, size2;
    GRB_TRY(compare_sizes(is_isomorphic, G1, G2, &size1, &size2, msg));
    if(*is_isomorphic == false) {
      return GrB_SUCCESS;
    }

    GrB_Vector degrees1, degrees2;
    GrB_Index nvals1, nvals2;
    GRB_TRY(compute_degrees(is_isomorphic, G1, G2, degrees1, degrees2, &nvals1, &nvals2, msg));
    if(*is_isomorphic == false) {
      return GrB_SUCCESS;
    }

    // Compare the degrees
    int64_t *values1 = malloc(nvals1 * sizeof(int64_t));
    int64_t *values2 = malloc(nvals2 * sizeof(int64_t));
    allocate_and_extract_tuples(values1, values2, nvals1, nvals2, degrees1, degrees2, msg);
    sort_and_compare(is_isomorphic, values1, values2, nvals1, nvals2);

    free(values1);
    free(values2);
    return (GrB_SUCCESS);
}

int do_fast(bool *is_isomorphic,
                 const LAGraph_Graph G1,
                 const LAGraph_Graph G2,
                 char *msg) {
    GrB_Index size1, size2;
    GRB_TRY(compare_sizes(is_isomorphic, G1, G2, &size1, &size2, msg));
    if(*is_isomorphic == false) {
        return GrB_SUCCESS;
    }

    GrB_Vector degrees1, degrees2;
    GrB_Index degrees_nvals1, degrees_nvals2;
    GRB_TRY(compute_degrees(is_isomorphic, G1, G2, degrees1, degrees2, &degrees_nvals1, &degrees_nvals2, msg));
    if(*is_isomorphic == false) {
        return GrB_SUCCESS;
    }

    GrB_Vector triangles1 = malloc(sizeof(GrB_Vector)), triangles2 = malloc(sizeof(GrB_Vector));
    GrB_Index triangles_nvals1, triangles_nvals2;
    GRB_TRY(compute_triangles(is_isomorphic, G1, triangles1, &triangles_nvals1, msg));
    GRB_TRY(compute_triangles(is_isomorphic, G2, triangles2, &triangles_nvals2, msg));
    if(*is_isomorphic == false) {
        return GrB_SUCCESS;
    }

    GRB_TRY(compare_nvals(is_isomorphic, triangles1, triangles2, &triangles_nvals1, &triangles_nvals2, &degrees_nvals1, &degrees_nvals2, size1, size2, msg));
}