//
// Created by jorge on 9/15/23.
//

#define LG_FREE_WORK                \
{                                   \
    GrB_free (&L) ;              \
}

#include "LG_internal.h"

int LAGraph_IsTournament(
        // output
        bool *is_tournament,

        // input
        LAGraph_Graph G,
        char *msg) {
    LG_CLEAR_MSG;

    // Declare items
    GrB_Index nvals, nrows, ncols;
    GrB_Matrix L = NULL;

    // Get matrix dimensions
    GRB_TRY(GrB_Matrix_nvals(&nvals, G->A));
    GRB_TRY(GrB_Matrix_nrows(&nrows, G->A));
    GRB_TRY(GrB_Matrix_ncols(&ncols, G->A));
    LAGRAPH_TRY(LAGraph_Cached_NSelfEdges(G, msg));

    // Find whether we can quickly assert that this is not a is_tournament graph
    if (nvals != nrows * (ncols - 1) / 2.0 || G->nself_edges > 0) {
        (*is_tournament) = false;
        return (GrB_SUCCESS);
    }

    GRB_TRY(GrB_Matrix_new(&L, GrB_FP64, nrows, ncols));
    GRB_TRY(LAGraph_Cached_AT(G, msg));
//    Grb_DESC_ST0
    GRB_TRY(GrB_select(L, G->A, NULL, GrB_TRIL, G->AT, -1, GrB_DESC_S));
    GxB_print(L, 3);
    GRB_TRY(GrB_Matrix_nvals(&nvals, L));
    (*is_tournament) = nvals == 0;

// Transpose matrix
//    GrB_Matrix A_transposed;
//    GRB_TRY(GrB_Matrix_new(&A_transposed, GrB_UINT64, ncols, nrows));
//    GRB_TRY(GrB_transpose(A_transposed, GrB_NULL, GrB_NULL, G->A, GrB_NULL));
//    GxB_print(A_transposed, 3);
//
//    // Create selector for upper triangle
//    GrB_Matrix triu;
//    GRB_TRY(GrB_Matrix_new(&triu, GrB_UINT64, ncols, nrows));
//    GRB_TRY(GrB_select(triu, NULL, NULL, GrB_TRIL, A_transposed, -1, NULL));
//    GxB_print(triu, 3);
//
//    // Mask the matrix
//    GrB_Matrix C;
//    GRB_TRY(GrB_Matrix_new(&C, GrB_UINT64, ncols, nrows));
//    GRB_TRY(GrB_mxm(C, G->A, GrB_NULL, GrB_PLUS_TIMES_SEMIRING_UINT64, triu, A_transposed, GrB_NULL));
//
//    // Count non-zero values
//    GrB_Index nvals_after;
//    GRB_TRY(GrB_Matrix_nvals(&nvals_after, C));
//    printf("nvals: %lu\n", nvals_after);
//    printf("\nis tournament 2: %s\n", nvals_after == 0 ? "true" : "false");
//    (*is_tournament) = nvals_after == 0;

    LG_FREE_WORK;
    return (GrB_SUCCESS);
}