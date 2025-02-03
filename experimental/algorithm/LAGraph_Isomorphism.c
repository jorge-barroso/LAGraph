//
// Created by jorge on 01/25/25.
//


#include <GraphBLAS.h>
#include <stdbool.h>

enum Level { FASTEST, FAST, SLOW, SLOWEST };

#define LG_FREE_WORK                \
{                                   \
    GrB_free (&degrees1) ;              \
    GrB_free (&degrees2) ;              \
}
#include "LG_internal.h"



int fastest_check(bool *is_isomorphic, const LAGraph_Graph G1, const LAGraph_Graph G2, GrB_Vector degrees1, GrB_Vector degrees2, char *msg);

int LAGraph_IsIsomorphic(
    // output
    bool *is_isomorphic,
    // inputs
    LAGraph_Graph G1,
    LAGraph_Graph G2,
    enum Level level,
    char *msg) {
    LG_CLEAR_MSG;

    (*is_isomorphic) = true;

   // TODO have to make them not null!!
   GrB_Vector degrees1, degrees2;
   GRB_TRY(fastest_check(is_isomorphic, G1, G2, degrees1, degrees2, msg));
   if(!(*is_isomorphic) || level == FASTEST) {
     return (GrB_SUCCESS);
   }

//    GrB_Vector triangles1 = NULL, triangles2 = NULL;
//    GRB_TRY(LAGraph_Tricount(&triangles1, G1, msg));
//    GRB_TRY(LAGraph_Tricount(&triangles2, G2, msg));
//
//    GrB_Index tri_nvals1, tri_nvals2;
//    GRB_TRY(GrB_Vector_nvals(&tri_nvals1, triangles1));
//    GRB_TRY(GrB_Vector_nvals(&tri_nvals2, triangles2));
//    if (tri_nvals1 != tri_nvals2) {
//        (*is_isomorphic) = false;
//        return GrB_SUCCESS;
//    }

    return (GrB_SUCCESS);
}

int compare_int64(const void *a, const void *b) {
    int64_t va = *(const int64_t*)a;
    int64_t vb = *(const int64_t*)b;
    // 1, 0 or -1 depending on whether va is greater, equal or less than vb
    return (va > vb) - (va < vb);
}

int fastest_check(bool *is_isomorphic,
                 const LAGraph_Graph G1,
                 const LAGraph_Graph G2,
                 GrB_Vector degrees1,
                 GrB_Vector degrees2,
                 char *msg) {
    if (G1->kind != G2->kind) {
        return false;
    }

    // We only need rows because Ajd Matrices are squared by definition
    GrB_Index size1, size2;
    GRB_TRY(GrB_Matrix_nrows(&size1, G1->A));
    GRB_TRY(GrB_Matrix_nrows(&size2, G2->A));
    if (size1 != size2) {
      (*is_isomorphic) = false;
      return GrB_SUCCESS;
    }

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
    GrB_Index nvals1, nvals2;
    GrB_Vector_nvals(&nvals1, degrees1);
    GrB_Vector_nvals(&nvals2, degrees2);
    if(nvals1 != nvals2) {
      return false;
    }

    // Compare the degrees
    int64_t *values1 = malloc(nvals1 * sizeof(int64_t));
    int64_t *values2 = malloc(nvals2 * sizeof(int64_t));
    if (values1 == NULL || values2 == NULL) {
        if (values1) free(values1);
        if (values2) free(values2);
        return GrB_OUT_OF_MEMORY;
    }
    GrB_Index actual_nvals;
    GRB_TRY(GrB_Vector_extractTuples(NULL, values1, &actual_nvals, degrees1));
    GRB_TRY(GrB_Vector_extractTuples(NULL, values2, &actual_nvals, degrees2));
    qsort(values1, nvals1, sizeof(int64_t), compare_int64);
    qsort(values2, nvals2, sizeof(int64_t), compare_int64);
    bool equal = true;
    for (GrB_Index i = 0; i < nvals1; i++) {
        if (values1[i] != values2[i]) {
            (*is_isomorphic) = false;
            break;
        }
    }

    free(values1);
    free(values2);
    return (GrB_SUCCESS);
}