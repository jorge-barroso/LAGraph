//
// Created by jorge on 9/15/23.
//

#define LG_FREE_WORK                \
{                                   \
    GrB_free (&U) ;              \
}

#include "LG_internal.h"

bool is_tournament(LAGraph_Graph G,
                   char *msg) {
    GrB_Matrix A = G->A;
    GrB_Index nvals, nrows, ncols;
    GrB_Matrix U = NULL;
    GRB_TRY(GrB_Matrix_nvals(&nvals, A));
    GRB_TRY(GrB_Matrix_nrows(&nrows, A));
    GRB_TRY(GrB_Matrix_ncols(&ncols, A));
    if (nvals != nrows * floor((ncols - 1) / 2.0) || G->nself_edges) {
        return false;
    }

    GRB_TRY (GrB_select(U, NULL, NULL, GrB_TRIU, A, (int64_t) 1, NULL));
    GRB_TRY(GrB_Matrix_nvals(&nvals, U));

    return nvals == 0;
}