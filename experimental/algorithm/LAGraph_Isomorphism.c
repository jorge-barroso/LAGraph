//
// Created by jorge on 01/25/25.
//


#include <GraphBLAS.h>
#include <stdbool.h>

enum Level { FASTEST, FAST, SLOW, SLOWEST };

//#define LG_FREE_WORK                \
//{                                   \
//    GrB_free (&degrees1) ;              \
//    GrB_free (&degrees2) ;              \
//}
#define LG_FREE_WORK {}

#include "LG_internal.h"


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

    return (GrB_SUCCESS);
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