//
// Created by jorge on 9/15/23.
//

#define LG_FREE_WORK                \
{                                   \
    GrB_free (&U) ;              \
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
    GrB_Matrix A = G->A;
    GrB_Index nvals, nrows, ncols;
    GrB_Matrix U = NULL;

    // Get matrix dimensions
    GRB_TRY(GrB_Matrix_nvals(&nvals, A));
    GRB_TRY(GrB_Matrix_nrows(&nrows, A));
    GRB_TRY(GrB_Matrix_ncols(&ncols, A));

    // Find whether we can quickly assert that this is not a is_tournament graph
    printf("self edges? %lb", G->nself_edges);
    if (nvals != nrows * (ncols - 1) / 2.0 || G->nself_edges) {
        printf("\nis tournament 1: %b\n", false);
        (*is_tournament) = false;
        return (GrB_SUCCESS) ;
    }

    GRB_TRY (GrB_select(U, NULL, NULL, GrB_TRIU, A, (int64_t) 1, NULL));
    GRB_TRY(GrB_Matrix_nvals(&nvals, U));
    printf("\nis tournament 2: %b\n", nvals == 0);
    (*is_tournament) = nvals == 0;

    LG_FREE_WORK;
    return (GrB_SUCCESS) ;
}