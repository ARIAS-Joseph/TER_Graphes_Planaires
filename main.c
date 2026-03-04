#include <string.h>
#include <time.h>
#include "graph.h"
#include "planar_graph_creator.h"
#include "permutations.h"

void test_multiple_horton(Graph* g, const char* filename) {

    const int n = g->nb_edges;
    int *inv = calloc(n, sizeof(int));

    multiple_horton(g, inv, 1000);

    printf("\n=== RÉSULTATS ===\n");
    printf("Nombre de bases minimales distinctes trouvées: %d\n\n", g->nb_minimal_bases);

    for (int b = 0; b < g->nb_minimal_bases; b++) {
        printf("Minimal cycle basis %d:\n", b + 1);
        for (int i = 0; i < g->basis_dimension; i++) {
            printf("  Cycle %d: length = %d, edges = [", i + 1,
                   g->minimals_basis[b].cycles[i].length);
            for (int e = 0; e < g->nb_edges; e++) {
                if (g->minimals_basis[b].cycles[i].edges_ids[e] == 1) {
                    printf(" %d", e);
                }
            }
            printf(" ]\n");
        }
    }
    save_graph(g, filename);
    free(inv);
}

void test_planar_graph() {
    const int nb_vertex = 50;
    Graph* planar = create_planar_graph(nb_vertex, 10000);
    test_multiple_horton(planar, "planar_test.txt");
    free_graph(planar);
    Graph* outer_planar = create_outer_planar_graph(nb_vertex, 70);
    test_multiple_horton(outer_planar, "outer_planar_test.txt");
    free_graph(outer_planar);

}

void test_planar_circle() {
    const int nb_vertex = 50;
    Graph *g = test_circle(nb_vertex);

    save_graph(g, "circle_graph.txt");
    free_graph(g);
}

void test_tree() {
    const int nb_vertex = 50;
    Graph* g = create_graph();
    double created[nb_vertex][2];
    int i = 0;
    while (i < nb_vertex) {
        const double x = rand() % (nb_vertex*100);
        const double y = rand() % (nb_vertex*100);
        int already_exists = 0;
        for (int j = 0; j < i; j++) {
            if (created[j][0] == x && created[j][1] == y) {
                already_exists = 1;
                break;
            }
        }
        if (!already_exists) {
            created[i][0] = x;
            created[i][1] = y;
            create_vertex(g, x, y);
            i++;
        }
    }

    int *perm = malloc(g->nb_vertex * (g->nb_vertex - 1) * sizeof(int));
    create_all_edges(g->nb_vertex, perm);
    fisher_yates_shuffle(perm, g->nb_vertex * (g->nb_vertex - 1));

    create_tree(g, perm);

    save_graph(g, "tree_graph.txt");
    free_graph(g);
    free(perm);
}

void test_permutations() {
    const int n = 4; // Number of elements in the permutation
    int *inv = calloc(n, sizeof(int));
    int *perm = calloc(n, sizeof(int));
    inversion_to_permutation(inv, perm, n);
    printf("Inversion table: [");
    for (int i = 0; i < n; i++) {
        printf("%d ", inv[i]);
    }
    printf("] -> Permutation: [");
    for (int i = 0; i < n; i++) {
        printf("%d ", perm[i]);
    }
    printf("]\n");

    while(next_inversion_table(inv, n))
    {
        memset(perm, 0, n * sizeof(int));
        inversion_to_permutation(inv, perm, n);
        printf("Inversion table: [");
        for (int i = 0; i < n; i++) {
            printf("%d ", inv[i]);
        }
        printf("] -> Permutation: [");
        for (int i = 0; i < n; i++) {
            printf("%d ", perm[i]);
        }
        printf("]\n");
    }

    free(inv);
    free(perm);
}

void test_all_edges_creation() {
    const int nb_vertex = 5;
    int *perm = malloc(nb_vertex * (nb_vertex - 1) * sizeof(int));
    create_all_edges(nb_vertex, perm);
    printf("All edges for %d vertices:\n", nb_vertex);
    for (int i = 0; i < nb_vertex * (nb_vertex - 1); i += 2) {
        printf("Edge: (%d, %d)\n", perm[i], perm[i + 1]);
    }
    free(perm);
}

void test_fisher_yates() {
    const int n = 10;
    int *array = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        array[i] = i;
    }
    printf("Original array: [");
    for (int i = 0; i < n; i++) {
        printf("%d ", array[i]);
    }
    printf("]\n");

    fisher_yates_shuffle(array, n);
    printf("Shuffled array: [");
    for (int i = 0; i < n; i++) {
        printf("%d ", array[i]);
    }
    printf("]\n");

    free(array);
}

int main() {
    srand(time(NULL));
    test_permutations();
    // test_tree();
    // test_planar_circle();
    // test_planar_graph();
    //test_all_edges_creation();
    //test_fisher_yates();

    return 0;
}