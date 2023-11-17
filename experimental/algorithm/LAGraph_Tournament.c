//
// Created by jorge on 9/15/23.
//

#define LG_FREE_WORK                \
{                                   \
    GrB_free (&L) ;              \
    GrB_free (&diag) ;              \
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
    GrB_Matrix diag = NULL;

    // Get matrix dimensions
    GRB_TRY(GrB_Matrix_nvals(&nvals, G->A));
    GRB_TRY(GrB_Matrix_nrows(&nrows, G->A));
    GRB_TRY(GrB_Matrix_ncols(&ncols, G->A));

    // Transpose and cache
    GRB_TRY(LAGraph_Cached_AT(G, msg));

//        LAGRAPH_TRY(LAGraph_Cached_NSelfEdges(G, msg));
    GRB_TRY(GrB_Matrix_new(&diag, GrB_FP64, nrows, ncols));
    GRB_TRY(GrB_select(diag, G->A, NULL, GrB_DIAG, G->AT, 0, GrB_DESC_S));
    GRB_TRY(GrB_Matrix_nvals(&G->nself_edges, diag));

    // Find whether we can quickly assert that this is not a is_tournament graph
    if (nvals != nrows * (ncols - 1) / 2.0 || G->nself_edges > 0) {
        (*is_tournament) = false;
        return (GrB_SUCCESS);
    }

    GRB_TRY(GrB_Matrix_new(&L, GrB_FP64, nrows, ncols));
    GRB_TRY(GrB_select(L, G->A, NULL, GrB_TRIL, G->AT, -1, GrB_DESC_S));
    GxB_print(L, 3);
    GRB_TRY(GrB_Matrix_nvals(&nvals, L));
    (*is_tournament) = nvals == 0;

    LG_FREE_WORK;
    return (GrB_SUCCESS);
}