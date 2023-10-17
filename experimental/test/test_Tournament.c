#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>

bool tournament = false;
char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;

#define LEN 512
char filename[LEN + 1];

void test_IsTournament(void) {

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    LAGraph_Init(msg);
    GrB_Matrix Y = NULL, A = NULL;

    //--------------------------------------------------------------------------
    // test with the west0067 matrix
    //--------------------------------------------------------------------------

    // create the graph
    snprintf(filename, LEN, LG_DATA_DIR "%s", "square.mtx");
    FILE *f = fopen(filename, "r");
    TEST_CHECK (f != NULL);
    OK (LAGraph_MMRead(&A, f, msg));
    OK (fclose(f));
    OK (LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
    TEST_CHECK (A == NULL);    // A has been moved into G->A

    OK(LAGraph_IsTournament(&tournament, G, msg));

    TEST_CHECK(tournament == true);

    OK (LAGraph_Delete(&G, msg));
    LAGraph_Finalize(msg);
}