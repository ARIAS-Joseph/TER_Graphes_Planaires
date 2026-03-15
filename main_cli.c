#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "graph.h"
#include "permutations.h"
#include "planar_graph_creator.h"

static void reset_bases(Graph *g) {
    if (!g->minimals_basis) return;
    for (int b = 0; b < g->nb_minimal_bases; b++) {
        for (int i = 0; i < g->basis_dimension; i++) {
            free(g->minimals_basis[b].cycles[i].edges_ids);
            free(g->minimals_basis[b].cycles[i].edges_labels);
        }
        free(g->minimals_basis[b].cycles);
    }
    free(g->minimals_basis);
    g->minimals_basis = NULL;
    g->nb_minimal_bases = 0;
    g->basis_dimension = 0;
}

static void save_clean_horton_save(Graph *g, const char *file, int n) {

    save_graph(g, file);
    delete_graph(g);

    g = create_graph();
    load_graph(g, file);

    if (g->nb_edges >= 1) {
        int *inv = calloc(g->nb_edges, sizeof(int));
        if (!inv) { perror("calloc inv"); exit(1); }
        multiple_horton(g, inv, n);
        free(inv);
    }

    save_graph(g, file);
    delete_graph(g);
}

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage:\n"
        "  %s run_horton     <file> <n>\n"
        "  %s add_vertex     <file> <x> <y> [n_horton]\n"
        "  %s delete_vertex  <file> <vid> [n_horton]\n"
        "  %s add_edge       <file> <u> <v> [n_horton]\n"
        "  %s delete_edge    <file> <eid> [n_horton]\n"
        "  %s split_edges    <file> <k> [n_horton] <eid1> [eid2 ...]\n"
        "  %s generate_graph <file> <planar|outer_planar> <nb_vertex> <nb_edges> [n_horton]\n",
        prog, prog, prog, prog, prog, prog, prog);
    exit(1);
}

static int opt_int(int argc, char *argv[], int idx, int def) {
    if (idx < argc) return atoi(argv[idx]);
    return def;
}

int main(int argc, char *argv[]) {
    srand(time(NULL));
    if (argc < 3) usage(argv[0]);

    const char *action = argv[1];
    const char *file   = argv[2];

    if (strcmp(action, "generate_graph") == 0) {
        if (argc < 6) usage(argv[0]);
        const char *type    = argv[3];
        const int nb_vertex = atoi(argv[4]);
        const int nb_edges  = atoi(argv[5]);
        const int n_horton  = opt_int(argc, argv, 6, 10);

        Graph *g;
        if (strcmp(type, "outer_planar") == 0)
            g = create_outer_planar_graph(nb_vertex, nb_edges);
        else
            g = create_planar_graph(nb_vertex, nb_edges);

        if (g->nb_edges >= 1) {
            int *inv = calloc(g->nb_edges, sizeof(int));
            if (!inv) { perror("calloc"); exit(1); }
            multiple_horton(g, inv, n_horton);
            free(inv);
        }
        save_graph(g, file);
        delete_graph(g);
        return 0;
    }

    Graph *g = create_graph();
    load_graph(g, file);

    if (strcmp(action, "run_horton") == 0) {
        const int n = opt_int(argc, argv, 3, 10);
        if (g->nb_edges >= 1) {
            int *inv = calloc(g->nb_edges, sizeof(int));
            if (!inv) { perror("calloc"); exit(1); }
            multiple_horton(g, inv, n);
            free(inv);
        }
        save_graph(g, file);
        delete_graph(g);

    } else if (strcmp(action, "add_vertex") == 0) {
        if (argc < 5) usage(argv[0]);
        create_vertex(g, atof(argv[3]), atof(argv[4]));
        save_clean_horton_save(g, file, opt_int(argc, argv, 5, 10));

    } else if (strcmp(action, "delete_vertex") == 0) {
        if (argc < 4) usage(argv[0]);
        const int vid = atoi(argv[3]);
        for (int e = 0; e < g->nb_edges; e++) {
            if (g->edges[e].deleted) continue;
            if (g->edges[e].u == vid || g->edges[e].v == vid)
                delete_edge(g, e);
        }
        g->vertices[vid].deleted = 1;
        reset_bases(g);
        save_clean_horton_save(g, file, opt_int(argc, argv, 4, 10));

    } else if (strcmp(action, "add_edge") == 0) {
        if (argc < 5) usage(argv[0]);
        create_edge(g, atoi(argv[3]), atoi(argv[4]));
        save_clean_horton_save(g, file, opt_int(argc, argv, 5, 10));

    } else if (strcmp(action, "delete_edge") == 0) {
        if (argc < 4) usage(argv[0]);
        delete_edge(g, atoi(argv[3]));
        reset_bases(g);
        save_clean_horton_save(g, file, opt_int(argc, argv, 4, 10));

    } else if (strcmp(action, "split_edges") == 0) {
        /* split_edges <file> <k> [n_horton] <eid1> [eid2 ...] */
        if (argc < 5) usage(argv[0]);
        const int k = atoi(argv[3]);
        const int n_horton = opt_int(argc, argv, 4, 10);
        /* Les edge IDs commencent à argv[5] si n_horton fourni, argv[4] sinon.
         * Pour simplifier : on cherche tous les args >= argv[4] qui sont des IDs
         * valides. On suppose que n_horton est toujours fourni par le JS. */
        const int nb_eids  = argc - 5;
        if (nb_eids <= 0) { fprintf(stderr, "split_edges: aucun edge ID fourni\n"); exit(1); }
        int *eids = malloc(nb_eids * sizeof(int));
        if (!eids) { perror("malloc"); exit(1); }
        for (int i = 0; i < nb_eids; i++) eids[i] = atoi(argv[5 + i]);
        for (int i = 0; i < nb_eids; i++) split_edge(g, eids[i], k);
        free(eids);
        save_clean_horton_save(g, file, n_horton);

    } else {
        fprintf(stderr, "Action inconnue : %s\n", action);
        delete_graph(g);
        usage(argv[0]);
    }

    return 0;
}