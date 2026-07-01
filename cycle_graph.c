#include "cycle_graph.h"

#include <stdlib.h>

/**
 * @brief Compute the centroid of cycle c in basis mb of graph g.
 *
 * Each edge of the cycle contributes both its endpoints. In a simple cycle every
 * vertex is incident to exactly 2 edges, so each vertex is counted twice —
 * the average is still correct.
 *
 * @param g    Original graph.
 * @param mb   Minimal basis.
 * @param c    Index of the cycle within mb (0 .. dim-1).
 * @param[out] cx  X coordinate of the centroid.
 * @param[out] cy  Y coordinate of the centroid.
 */
static void cycle_centroid(const Graph *g, const Minimal_basis *mb, int c,
                           double *cx, double *cy)
{
    *cx = 0.0;
    *cy = 0.0;
    int count = 0;

    for (int e = 0; e < g->nb_edges; e++) {
        if (g->edges[e].deleted)          continue;
        if (!mb->cycles[c].edges_ids[e])  continue;

        int u = g->edges[e].u;
        int v = g->edges[e].v;
        *cx += g->vertices[u].x + g->vertices[v].x;
        *cy += g->vertices[u].y + g->vertices[v].y;
        count += 2;   /* each vertex counted once per incident cycle-edge */
    }

    if (count > 0) {
        *cx /= count;
        *cy /= count;
    }
}

/**
 * @brief Build the cycle graph Γ(mb) for a given minimal basis.
 *
 * Nodes  : one per cycle in mb, placed at the cycle's centroid in g.
 * Edges  : an edge (i, j) is added when cycles i and j share at least one
 *          edge of the original graph g.
 *
 * @param g    Original graph.
 * @param mb   Minimal basis of g.
 * @param dim  Number of cycles in the basis (g->basis_dimension).
 * @return     Freshly allocated cycle graph (caller must delete_graph() it).
 */
static Graph *build_cycle_graph(const Graph *g, const Minimal_basis *mb, int dim)
{
    Graph *cg = create_graph();

    for (int c = 0; c < dim; c++) {
        double cx, cy;
        cycle_centroid(g, mb, c, &cx, &cy);
        create_vertex(cg, cx, cy);
    }

    for (int i = 0; i < dim; i++) {
        for (int j = i + 1; j < dim; j++) {
            int shared = 0;
            for (int e = 0; e < g->nb_edges && !shared; e++) {
                if (mb->cycles[i].edges_ids[e] && mb->cycles[j].edges_ids[e])
                    shared = 1;
            }
            if (shared)
                create_edge(cg, i, j);
        }
    }

    return cg;
}

/**
 * @brief Return 1 if any minimal basis of cg has every edge appearing in
 *        at most 2 cycles, 0 otherwise.
 *
 * This is a necessary condition for a face basis in a planar graph: each edge
 * borders exactly 2 faces, so it can appear in at most 2 basis cycles.
 */
static int any_basis_satisfies_planarity_condition(const Graph *cg)
{
    if (cg->nb_minimal_bases == 0) return 0;

    for (int b = 0; b < cg->nb_minimal_bases; b++) {
        const Minimal_basis *mb = &cg->minimals_basis[b];
        int ok = 1;

        for (int e = 0; e < cg->nb_edges && ok; e++) {
            if (cg->edges[e].deleted) continue;

            int cnt = 0;
            for (int c = 0; c < cg->basis_dimension; c++)
                cnt += mb->cycles[c].edges_ids[e];

            if (cnt == 0) continue;
            if (cnt > 2) ok = 0;
        }

        if (ok) return 1;
    }
    return 0;
}

/**
 * @brief Return 1 if the cycle graph cg appears to be planar, 0 otherwise.
 *
 * If cg has no edges it is trivially planar (independent set of nodes).
 * Otherwise Horton is run with nb_permutations random numberings and the
 * planarity condition (every edge in ≤ 2 basis cycles) is checked.
 */
static int cycle_graph_is_planar(Graph *cg, int nb_permutations)
{
    if (cg->nb_edges == 0)
        return 1;

    int *inv_e = malloc(cg->nb_edges * sizeof(int));
    int *inv_v = malloc(cg->nb_vertices * sizeof(int));
    if (!inv_e || !inv_v) {
        free(inv_e);
        free(inv_v);
        perror("cycle_graph_is_planar: malloc");
        return 0;
    }

    multiple_horton(cg, inv_e, inv_v, nb_permutations);

    free(inv_e);
    free(inv_v);

    if (cg->basis_dimension == 0)
        return 1;

    return any_basis_satisfies_planarity_condition(cg);
}

int check_cycle_graphs_planarity(Graph *g, int nb_permutations)
{
    if (!g || g->nb_minimal_bases == 0 || g->basis_dimension == 0)
        return 0;

    int nb_planar = 0;
    int dim = g->basis_dimension;

    for (int b = 0; b < g->nb_minimal_bases; b++) {
        Graph *cg = build_cycle_graph(g, &g->minimals_basis[b], dim);

        if (cycle_graph_is_planar(cg, nb_permutations))
            nb_planar++;

        delete_graph(cg);
    }

    if (nb_planar == 0) return 0;
    if (nb_planar == g->nb_minimal_bases) return 2;
    return 1;
}