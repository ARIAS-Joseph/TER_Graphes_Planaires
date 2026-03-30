#include <stdio.h>
#include <stdlib.h>
#include <direct.h>
#include <time.h>
#include <errno.h>

#include "graph.h"
#include "planar_graph_creator.h"

void process_and_save(char* type, Graph *g, int id) {
    FILE *csv = fopen("result.csv", "a");
    if (!csv) {
        printf("Error opening csv file\n");
        return;
    }

    int *inv = calloc(g->nb_edges, sizeof(int));
    multiple_horton(g, inv, 1000);

    const int dir_all = _mkdir("all_graphs");
    if (dir_all == -1 && errno != EEXIST) {
        printf("Error creating directory all_graphs\n");
        free(inv);
        return;
    }

    char dir_buffer[100];
    snprintf(dir_buffer, 100, "./all_graphs/%s_%dbasis", type, g->nb_minimal_bases);
    const int dir = _mkdir(dir_buffer);
    if (dir == -1 && errno != EEXIST) {
        printf("Error creating directory %s\n", dir_buffer);
        free(inv);
        return;
    }

    char buffer[100];
    snprintf(buffer, 100, "./all_graphs/%s_%dbasis/%s_%dvertices_%dedges_%dbasis_id%d.txt", type, g->nb_minimal_bases, type, g->nb_vertex, g->nb_edges, g->nb_minimal_bases, id);
    save_graph(g, buffer);

    fprintf(csv,"%s, %d, %d, %d, %d, %d\n", type, g->nb_vertex, g->nb_edges, g->nb_minimal_bases, g->face_basis != -1, g->face_basis_outer != -1);
    fclose(csv);
    delete_graph(g);
    free(inv);
}

void study_graph(const int nb_vertex, const int nb_test) {

    for (int i = 0; i < nb_test; i++) {
        int nb_edges = (rand() % ((3*nb_vertex - 6) - nb_vertex + 1)) + nb_vertex;
        Graph* planar = create_planar_graph(nb_vertex, nb_edges);
        process_and_save("PLANAR", planar, i);
    }
}

int main(int argc, char* argv[]) {
    freopen("NUL", "w", stdout);
    FILE *csv = fopen("result.csv", "w+");
    if (!csv) {
        printf("Error opening csv file\n");
        return 0;
    }
    fprintf(csv,"Type, nb_vertices, nb_edges, nb_minimal_basis, has_faces, has_faces_outer\n");
    fclose(csv);
    srand(time(NULL));
    for (int i = 10; i < 21; i++) {
        study_graph(i, 1000);
    }
}