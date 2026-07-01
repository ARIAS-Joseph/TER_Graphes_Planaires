#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "dfs_graph.h"

static void print_dfs_graph(const DfsGraph *dfs) {
    printf("=== Sommets ===\n\n");
    printf("v\tparent\tdfs_num\n");

    for (int v = 0; v < dfs->vertices_count; ++v) {
        printf("%d\t%d\t%d\n", v, dfs->vertices[v].parent, dfs->vertices[v].dfs_num);
    }

    printf("\n=== Arcs DFS ===\n\n");
    printf("id\tfrom\tto\ttype\n");

    for (int i = 0; i < dfs->edge_count; ++i) {
        const DfsEdge *e = &dfs->edges[i];
        printf("%d\t%d\t%d\t%s\n", e->id, dfs->vertices[e->from].dfs_num,
            dfs->vertices[e->to].dfs_num,
            e->type == TREE_EDGE ? "TREE" : "BACK");
    }

    printf("\n=== Successeurs arcs ===\n\n");
    for (int e = 0; e < dfs->edge_count; ++e) {
        for (int i = 0; i < dfs->edges[e].nb_successors; i++) {
            const int succ_id = dfs->edges[e].successors[i];
            printf("Arc %d->%d (id=%d) successeur: %d->%d (id=%d)\n",
                dfs->vertices[dfs->edges[e].from].dfs_num,
                dfs->vertices[dfs->edges[e].to].dfs_num,
                e,
                dfs->vertices[dfs->edges[succ_id].from].dfs_num,
                dfs->vertices[dfs->edges[succ_id].to].dfs_num,
                succ_id);
        }
    }

    printf("\n=== Successeurs sommets ===\n\n");
    for (int v = 0; v < dfs->vertices_count; ++v) {
        printf("\nSommet %d (id=%d) successeurs: ", dfs->vertices[v].dfs_num, v);
        for (int i = 0; i < dfs->vertices[v].nb_children; i++) {
            const int child_id = dfs->vertices[v].children[i];
            printf("%d ", dfs->vertices[child_id].dfs_num);
        }
    }
}

static void print_low_values_phi(const DfsGraph *dfs) {
    printf("\n\n=== Low values et phi arêtes ===\n\n");
    printf("e\tlow1\tlow2\tphi\tsingulière\n");

    for (int i = 0; i < dfs->edge_count; ++i) {
        const DfsEdge *e = &dfs->edges[i];
        printf("%d->%d\t%d\t%d\t%d\t%s\n", dfs->vertices[e->from].dfs_num,
            dfs->vertices[e->to].dfs_num, e->low1, e->low2, e->phi, e->singular? "Oui" : "Non");
    }

    printf("\n\n=== Phi sommets ===\n\n");
    printf("v (dfs num)\tphi\n");

    for (int i = 0; i < dfs->vertices_count; ++i) {
        printf("%d\t", dfs->vertices[i].dfs_num);
        for (int j = 0; j < dfs->vertices[i].nb_out_edges; ++j) {
            const int edge = dfs->vertices[i].out_edges[j];
            printf("%d->%d (phi=%d)\t", dfs->vertices[dfs->edges[edge].from].dfs_num, dfs->vertices[dfs->edges[edge].to].dfs_num, dfs->edges[edge].phi);
        }
        printf("\n");
    }
}

static void print_singular_sets(const DfsGraph *dfs) {

    printf("\n=== E0 ===\n\n");

    for (int i = 0; i < dfs->E0_size; ++i) {

        int edge_id = dfs->E0[i];

        const DfsEdge *e = &dfs->edges[edge_id];

        printf("[%d->%d] ",
               dfs->vertices[e->from].dfs_num,
               dfs->vertices[e->to].dfs_num);
    }

    printf("\n");

    printf("\n=== Singular sets ===\n\n");

    for (int i = 0; i < dfs->nb_singular_sets; ++i) {

        const SingularSet *set = &dfs->singular_sets[i];

        printf("S(v=%d, low1=%d): ",
               dfs->vertices[set->vertex].dfs_num,
               set->low1);

        for (int j = 0; j < set->size; ++j) {

            int edge_id = set->edges[j];

            const DfsEdge *e = &dfs->edges[edge_id];

            printf("[%d->%d] ",
                   dfs->vertices[e->from].dfs_num,
                   dfs->vertices[e->to].dfs_num);
        }

        printf("\n");
    }
}

static void print_same_diff_and_diff(const DfsGraph *dfs) {
    printf("=== SAME ===\n");
    for (int i = 0; i < dfs->SAME.nb_groups; i++) {
        printf("groupe %d : { ", i);
        for (int j = 0; j < dfs->SAME.groups[i].size; j++) {
            int e = dfs->SAME.groups[i].edges[j];
            printf("%d->%d ", dfs->vertices[dfs->edges[e].from].dfs_num,
                               dfs->vertices[dfs->edges[e].to].dfs_num);
        }
        printf("}\n");
    }

    printf("=== DIFF ===\n");
    for (int i = 0; i < dfs->DIFF.nb_groups; i++) {
        int p = dfs->DIFF.partner[i];
        if (p == -1) printf("groupe %d : libre\n", i);
        else if (p > i) printf("(groupe %d, groupe %d)\n", i, p);
    }

    printf("=== SAME' ===\n");
    for (int i = 0; i < dfs->SAME_prime.nb_groups; i++) {
        printf("groupe %d : { ", i);
        for (int j = 0; j < dfs->SAME_prime.groups[i].size; j++) {
            int e = dfs->SAME_prime.groups[i].edges[j];
            printf("%d->%d ", dfs->vertices[dfs->edges[e].from].dfs_num,
                               dfs->vertices[dfs->edges[e].to].dfs_num);
        }
        printf("}\n");
    }

    printf("=== DIFF' ===\n");
    for (int i = 0; i < dfs->DIFF_prime.nb_groups; i++) {
        int p = dfs->DIFF_prime.partner[i];
        if (p == -1) printf("groupe %d : libre\n", i);
        else if (p > i) printf("(groupe %d, groupe %d)\n", i, p);
    }
}

int main(void) {
    Graph *g = create_graph();

    for (int i = 0; i < 6; ++i) {
        create_vertex(g, 0.0, 0.0);
    }

    create_edge(g, 0, 1);
    create_edge(g, 1, 2);
    create_edge(g, 2, 3);
    create_edge(g, 3, 4);
    create_edge(g, 3, 5);
    create_edge(g, 4, 0);
    create_edge(g, 4, 1);
    create_edge(g, 2, 5);
    create_edge(g, 2, 0);

    printf("Graph created with %d vertices and %d edges.\n", g->nb_vertices, g->nb_edges);

    DfsGraph *dfs = build_dfs_graph(g);
    // print_dfs_graph(dfs);

    compute_low_values(dfs);
    build_phi_lists(dfs);
    //print_low_values_phi(dfs);

    build_singular_sets(dfs);
    //print_singular_sets(dfs);

    compute_same_diff(dfs);
    build_SAME_DIFF_prime(dfs);
    //print_same_diff_and_diff(dfs);

    const EmbeddingSet result = enumerate_embeddings(dfs);

    for (int i=0; i<result.count; i++) {
        printf("Embedding %d:\n", i+1);
        for (int v=0; v<g->nb_vertices; v++) {
            printf("Vertex %d: ", v);
            for (int j=0; j<result.embeddings[i].rotation_lengths[v]; j++) {
                printf("%d ", result.embeddings[i].rotation_system[v][j]);
            }
            printf("\n");
        }
    }

    free_dfs_graph(dfs);
    delete_graph(g);

    return 1;
}