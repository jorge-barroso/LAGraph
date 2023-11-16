#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>

#define LEN 512

char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;
GrB_Matrix Y = NULL, A = NULL;
bool tournament = false;
char filename[LEN + 1];

typedef struct {
    const bool is_tournament;
    const char *name;
} matrix_info;

const matrix_info files[] =
        {
                {true,  "square_tournament.mtx"},
                {false, "square_diagonal_vals.mtx"}, // when i == j (diagonal), cell(i, j) should be 0
        };

void test_IsTournament(void) {

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------
    LAGraph_Init(msg);

    for (int i = 0; i < sizeof(files) / sizeof(files[0]); ++i) {


        const char *aname = files[i].name;
        const bool expectedReturn = files[i].is_tournament;
        if (strlen(aname) == 0) break;
        printf("\n================================== %s: ==================================\n", aname);
        TEST_CASE (aname);
        // create the graph
        snprintf(filename, LEN, LG_DATA_DIR "%s", aname);
        FILE *f = fopen(filename, "r");
        TEST_CHECK (f != NULL);
        OK (LAGraph_MMRead(&A, f, msg));
        OK (fclose(f));
        OK (LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
        TEST_CHECK (A == NULL);    // A has been moved into G->A

        OK(LAGraph_IsTournament(&tournament, G, msg));

        printf("\nis tournament: %b\n", tournament);
        TEST_CHECK(tournament == expectedReturn);

        OK (LAGraph_Delete(&G, msg));
    }

    LAGraph_Finalize(msg);
}

//----------------------------------------------------------------------------
// the make program is created by acutest, and it runs a list of tests:
//----------------------------------------------------------------------------

TEST_LIST =
        {
                {"IsTournament", test_IsTournament},
                {NULL, NULL}
        };