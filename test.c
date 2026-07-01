#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "dfs_graph.h"

static void dump_pipeline(const char *label, Graph *g) {
    printf("\n########## %s ##########\n", label);

    printf("Adjacency:\n");
    for (int v = 0; v < g->nb_vertices; v++) {
        printf("  %d:", v);
        for (int i = 0; i < g->neighbors[v].count; i++)
            printf(" %d", g->neighbors[v].neighbors[i]);
        printf("\n");
    }

    DfsGraph *dfs = build_dfs_graph(g);
    compute_low_values(dfs);
    build_phi_lists(dfs);
    build_singular_sets(dfs);
    compute_same_diff(dfs);
    build_SAME_DIFF_prime(dfs);

    printf("\nDFS tree (parent_edge):\n");
    for (int v = 0; v < dfs->vertices_count; v++) {
        int pe = dfs->vertices[v].parent_edge;
        if (pe == -1) printf("  v=%d root (dfs_num=%d)\n", v, dfs->vertices[v].dfs_num);
        else printf("  v=%d dfs_num=%d  parent_edge=%d (%d->%d)\n",
                    v, dfs->vertices[v].dfs_num, pe,
                    dfs->edges[pe].from, dfs->edges[pe].to);
    }

    printf("\nSAME' groups:\n");
    for (int i = 0; i < dfs->SAME_prime.nb_groups; i++) {
        printf("  group %d: {", i);
        for (int j = 0; j < dfs->SAME_prime.groups[i].size; j++) {
            int e = dfs->SAME_prime.groups[i].edges[j];
            printf(" %d->%d", dfs->edges[e].from, dfs->edges[e].to);
        }
        printf(" }  partner=%d\n", dfs->DIFF_prime.partner[i]);
    }

    int d = 0, s = 0;
    int *done = calloc(dfs->SAME_prime.nb_groups, sizeof(int));
    for (int i = 0; i < dfs->SAME_prime.nb_groups; i++) {
        if (done[i]) continue;
        int p = dfs->DIFF_prime.partner[i];
        if (p == -1) { s++; done[i] = 1; }
        else { d++; done[i] = done[p] = 1; }
    }
    free(done);
    printf("\nd=%d  s=%d\n", d, s);

    printf("\nSingular sets (RS) and h(x):\n");
    for (int si = 0; si < dfs->nb_singular_sets; si++) {
        SingularSet *set = &dfs->singular_sets[si];
        printf("  set %d (vertex dfs_num=%d, low1=%d, g(x)=%d): {",
               si, dfs->vertices[set->vertex].dfs_num, set->low1, set->size);
        for (int j = 0; j < set->size; j++) {
            int e = set->edges[j];
            printf(" %d->%d", dfs->edges[e].from, dfs->edges[e].to);
        }
        printf(" }\n");
    }

    EmbeddingSet result = enumerate_embeddings(dfs);
    printf("\n>>> enumerate_embeddings count = %d\n", result.count);

    free_embedding_set(&result);
    free_dfs_graph(dfs);
}

int main(void) {

    /* Version A : ordre observé via create_planar_graph (4 embeddings) */
    Graph *gA = create_graph();
    for (int i = 0; i < 6; i++) create_vertex(gA, 0, 0);
    create_edge(gA, 0, 3);
    create_edge(gA, 0, 5);
    create_edge(gA, 0, 1);
    create_edge(gA, 1, 2);
    create_edge(gA, 1, 4);
    create_edge(gA, 1, 3);
    create_edge(gA, 2, 3);
    create_edge(gA, 2, 4);
    create_edge(gA, 3, 5);
    for (int i = 0; i < 4; i++) delete_vertex(gA, 0);
    for (int i = 0; i < 4; i++) create_vertex(gA, 0, 0);
    create_edge(gA, 0, 3);
    create_edge(gA, 0, 5);
    create_edge(gA, 0, 1);
    create_edge(gA, 1, 2);
    create_edge(gA, 1, 4);
    create_edge(gA, 1, 3);
    create_edge(gA, 2, 3);
    create_edge(gA, 2, 4);
    create_edge(gA, 3, 5);
    dump_pipeline("VERSION A (create_planar_graph order)", gA);
    delete_graph(gA);

    /* Version B : ordre manuel (8 embeddings) */
    Graph *gB = create_graph();
    for (int i = 0; i < 6; i++) create_vertex(gB, 0, 0);
    create_edge(gB, 0, 1);
    create_edge(gB, 0, 3);
    create_edge(gB, 0, 5);
    create_edge(gB, 1, 2);
    create_edge(gB, 1, 3);
    create_edge(gB, 1, 4);
    create_edge(gB, 2, 3);
    create_edge(gB, 2, 4);
    create_edge(gB, 3, 5);
    dump_pipeline("VERSION B (manual order)", gB);
    delete_graph(gB);

    return 0;
}