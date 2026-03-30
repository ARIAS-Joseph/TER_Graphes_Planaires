/**
 * @file graph.c
 * @brief Implementation of the planar graph data structure and associated algorithms.
 *
 * This file provides the full implementation of the Graph structure, including
 * vertex and edge management, face detection, shortest-path computation, and
 * Horton's minimal cycle basis algorithm.
 *
 * ### Self-contained invariant: topology mutations auto-invalidate stale data
 *
 * Any function that changes the graph topology (create_edge, delete_edge,
 * create_vertex, delete_vertex, split_edge) automatically calls the internal
 * helper invalidate_bases() before modifying the structure.  This guarantees:
 *
 *  - cycle basis vectors (edges_ids) are never read with a stale nb_edges size,
 *  - face data is not compared against a topology it was not computed for,
 *  - callers external to this file never need to manually reset derived state.
 *
 * Multiple consecutive calls to invalidate_bases() are safe: after the first
 * call nb_faces == 0 and minimals_basis == NULL, so subsequent calls are no-ops.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include "graph.h"
#include "permutations.h"

/**
 * @brief Invalidate all derived data that depends on the graph topology.
 *
 * Must be called whenever the topology changes (vertices or edges added/removed).
 * Frees:
 *  - all minimal cycle bases (edges_ids vectors sized to the old nb_edges),
 *  - all per-face path data (edges_ids, edges_labels, vertices_ids).
 *
 * Resets the corresponding counters so that algorithms recompute from scratch
 * on the new topology.
 *
 * This function is @b idempotent: after the first call nb_faces == 0 and
 * minimals_basis == NULL, so every subsequent call traverses empty loops.
 *
 * @note This function is static.  It is an implementation detail of graph.c and
 *       is called automatically by every topology-mutating public function.
 *
 * @param g Graph whose derived data must be invalidated.
 */
static void invalidate_bases(Graph *g) {
    if (!g) return;

    /* Free all minimal cycle bases */
    if (g->minimals_basis) {
        for (int b = 0; b < g->nb_minimal_bases; b++) {
            for (int i = 0; i < g->basis_dimension; i++) {
                free(g->minimals_basis[b].cycles[i].edges_ids);
                free(g->minimals_basis[b].cycles[i].edges_labels);
                free(g->minimals_basis[b].cycles[i].vertices_ids);
            }
            free(g->minimals_basis[b].cycles);
        }
        free(g->minimals_basis);
        g->minimals_basis = NULL;
        g->nb_minimal_bases = 0;
        g->basis_dimension = 0;
    }

    for (int i = 0; i < g->nb_faces; i++) {
        free(g->faces[i].edges_ids);
        free(g->faces[i].edges_labels);
        free(g->faces[i].vertices_ids);
        g->faces[i].edges_ids = NULL;
        g->faces[i].edges_labels = NULL;
        g->faces[i].vertices_ids = NULL;
    }
    g->nb_faces = 0;
    g->face_basis = -1;
}

/**
 * @brief Allocate and initialise an empty graph.
 *
 * Creates a new Graph with default dynamic-array capacities for vertices, edges,
 * neighbour lists, and Horton cycles.  All algorithm-result fields
 * (edge_indices, predecessors, distances, minimals_basis, face data) are
 * initialised to NULL or zero.
 *
 * @return Pointer to the newly created graph.  Exits on allocation failure.
 */
Graph *create_graph() {
    Graph *graph = calloc(1, sizeof(Graph));
    if (!graph) { perror("create_graph: graph"); exit(1); }

    graph->capacity_vertices = 4;
    graph->vertices = malloc(graph->capacity_vertices * sizeof(Vertex));
    if (!graph->vertices) { free(graph); perror("create_graph: vertices"); exit(1); }

    graph->neighbors = malloc(graph->capacity_vertices * sizeof(Neighbor_list));
    if (!graph->neighbors) {
        free(graph->vertices); free(graph);
        perror("create_graph: neighbors"); exit(1);
    }
    for (int i = 0; i < graph->capacity_vertices; i++) {
        graph->neighbors[i].neighbors = NULL;
        graph->neighbors[i].angles = NULL;
        graph->neighbors[i].count = 0;
        graph->neighbors[i].capacity = 0;
    }

    graph->capacity_edges = 4;
    graph->edges = malloc(graph->capacity_edges * sizeof(Edge));
    if (!graph->edges) {
        free(graph->neighbors); free(graph->vertices); free(graph);
        perror("create_graph: edges"); exit(1);
    }

    /* faces capacity: upper bound from Euler's formula for planar graphs */
    graph->faces = malloc(
        (2 + graph->capacity_edges + graph->capacity_vertices) * sizeof(Path));
    if (!graph->faces) {
        free(graph->edges); free(graph->neighbors);
        free(graph->vertices); free(graph);
        perror("create_graph: faces"); exit(1);
    }

    graph->capacity_horton_cycles = 4;
    graph->horton_cycles = malloc(graph->capacity_horton_cycles * sizeof(Path));
    if (!graph->horton_cycles) {
        free(graph->faces); free(graph->edges);
        free(graph->neighbors); free(graph->vertices); free(graph);
        perror("create_graph: horton_cycles"); exit(1);
    }

    graph->face_basis = -1;
    graph->face_basis_outer = -1;
    graph->edge_indices = NULL;
    graph->predecessors = NULL;
    graph->distances = NULL;
    graph->nb_minimal_bases = 0;
    graph->minimals_basis = NULL;
    graph->basis_dimension = 0;

    return graph;
}

/**
 * @brief Free all memory associated with a graph.
 *
 * Releases vertices, edges, neighbour lists (including the angles arrays),
 * algorithm matrices (edge_indices, predecessors, distances), face path arrays,
 * Horton cycles, minimal cycle bases, and the graph structure itself.
 *
 * @param g Graph to destroy.  No-op if NULL.
 */
void delete_graph(Graph *g) {
    if (!g) return;

    free(g->vertices);
    free(g->edges);

    if (g->neighbors) {
        for (int i = 0; i < g->nb_vertex; i++) {
            free(g->neighbors[i].neighbors);
            free(g->neighbors[i].angles);
        }
        free(g->neighbors);
    }

    free(g->edge_indices);
    free(g->predecessors);
    free(g->distances);

    if (g->faces) {
        for (int i = 0; i < g->nb_faces; i++) {
            free(g->faces[i].edges_ids);
            free(g->faces[i].edges_labels);
            free(g->faces[i].vertices_ids);
        }
        free(g->faces);
    }

    if (g->horton_cycles) {
        for (int i = 0; i < g->nb_horton_cycles; i++) {
            free(g->horton_cycles[i].edges_ids);
            free(g->horton_cycles[i].edges_labels);
            free(g->horton_cycles[i].vertices_ids);
        }
        free(g->horton_cycles);
    }

    if (g->minimals_basis) {
        for (int b = 0; b < g->nb_minimal_bases; b++) {
            for (int i = 0; i < g->basis_dimension; i++) {
                free(g->minimals_basis[b].cycles[i].edges_ids);
                free(g->minimals_basis[b].cycles[i].edges_labels);
                free(g->minimals_basis[b].cycles[i].vertices_ids);
            }
            free(g->minimals_basis[b].cycles);
        }
        free(g->minimals_basis);
    }

    free(g);
}

/**
 * @brief Save the graph to a text file.
 *
 * File format (all IDs are re-numbered sequentially, deleted elements skipped):
 * @code
 * V E M D face_basis
 * seq_id  x  y          (V vertex lines)
 * seq_id  u  v          (E edge lines)
 * eid ...               (M * D cycle lines: space-separated edge IDs per cycle)
 * @endcode
 *
 * @param g        Graph to save.
 * @param filename Output file path.
 */
void save_graph(const Graph *g, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) { perror("save_graph: fopen"); exit(1); }
    if (!g || !g->vertices || !g->edges) {
        fclose(file); perror("save_graph: invalid graph"); exit(1);
    }

    /* Build remapping tables: old index → sequential ID, -1 if deleted */
    int *vid_map = calloc(g->nb_vertex, sizeof(int));
    int *eid_map = calloc(g->nb_edges,  sizeof(int));
    if (!vid_map || !eid_map) { perror("save_graph: maps"); exit(1); }

    int nb_v = 0, nb_e = 0;
    for (int i = 0; i < g->nb_vertex; i++)
        vid_map[i] = g->vertices[i].deleted ? -1 : nb_v++;
    for (int i = 0; i < g->nb_edges; i++)
        eid_map[i] = g->edges[i].deleted    ? -1 : nb_e++;

    if (g->minimals_basis)
        fprintf(file, "%d %d %d %d %d\n",
                nb_v, nb_e, g->nb_minimal_bases, g->basis_dimension, g->face_basis);
    else
        fprintf(file, "%d %d 0 0 -1\n", nb_v, nb_e);

    for (int i = 0; i < g->nb_vertex; i++) {
        if (g->vertices[i].deleted) continue;
        fprintf(file, "%d %lf %lf\n", vid_map[i], g->vertices[i].x, g->vertices[i].y);
    }

    for (int i = 0; i < g->nb_edges; i++) {
        if (g->edges[i].deleted) continue;
        const int nu = vid_map[g->edges[i].u];
        const int nv = vid_map[g->edges[i].v];
        if (nu < 0 || nv < 0) continue;
        fprintf(file, "%d %d %d\n", eid_map[i], nu, nv);
    }

    if (g->minimals_basis) {
        for (int b = 0; b < g->nb_minimal_bases; b++) {
            for (int i = 0; i < g->basis_dimension; i++) {
                for (int e = 0; e < g->nb_edges; e++) {
                    if (g->edges[e].deleted) continue;
                    if (g->minimals_basis[b].cycles[i].edges_ids[e] == 1)
                        fprintf(file, "%d ", eid_map[e]);
                }
                fprintf(file, "\n");
            }
        }
    }

    free(vid_map);
    free(eid_map);
    fclose(file);
}

/**
 * @brief Load a graph from a file produced by save_graph().
 *
 * @p g must be a freshly created, empty graph.  After loading, the is_faces
 * flag of each basis is restored from the face_basis field in the file header.
 *
 * @note create_edge() and create_vertex() are called internally during loading
 *       and will invoke invalidate_bases() on each call — this is harmless
 *       since no bases have been added yet at that point.
 *
 * @param g        Destination graph (must be empty).
 * @param filename File to read.
 */
void load_graph(Graph *g, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) { perror("load_graph: fopen"); exit(1); }

    int nb_vertex, nb_edges, nb_minimal_bases, basis_dimension, face_basis;
    if (fscanf(file, "%d %d %d %d %d\n",
               &nb_vertex, &nb_edges, &nb_minimal_bases,
               &basis_dimension, &face_basis) != 5) {
        fclose(file); perror("load_graph: header"); exit(1);
    }

    fprintf(stdout, "Loading vertices...\n"); fflush(stdout);
    for (int i = 0; i < nb_vertex; i++) {
        int id; double x, y;
        if (fscanf(file, "%d %lf %lf\n", &id, &x, &y) != 3) {
            fclose(file); perror("load_graph: vertex"); exit(1);
        }
        create_vertex(g, x, y);
    }

    fprintf(stdout, "Loading edges...\n"); fflush(stdout);
    for (int i = 0; i < nb_edges; i++) {
        int id, u, v;
        if (fscanf(file, "%d %d %d", &id, &u, &v) != 3) {
            fclose(file); perror("load_graph: edge"); exit(1);
        }
        create_edge(g, u, v);
    }
    { char tmp[256]; fgets(tmp, sizeof(tmp), file); } /* consume trailing newline */

    if (nb_minimal_bases > 0 && basis_dimension > 0) {
        /*
         * Directly populate the bases without going through the public API,
         * since bases are being restored (not computed), so invalidate_bases()
         * must NOT be triggered here.
         */
        g->minimals_basis = calloc(nb_minimal_bases, sizeof(Minimal_basis));
        if (!g->minimals_basis) { perror("load_graph: minimals_basis"); exit(1); }
        g->nb_minimal_bases = nb_minimal_bases;
        g->basis_dimension  = basis_dimension;
        g->face_basis       = face_basis;

        fprintf(stdout, "Loading bases (nb=%d dim=%d)...\n",
                nb_minimal_bases, basis_dimension); fflush(stdout);

        for (int b = 0; b < nb_minimal_bases; b++) {
            g->minimals_basis[b].cycles = calloc(basis_dimension, sizeof(Path));
            if (!g->minimals_basis[b].cycles) { perror("load_graph: cycles"); exit(1); }

            /* Restore the is_faces flag: face_basis is 1-based, b is 0-based */
            g->minimals_basis[b].is_faces =
                (face_basis > 0 && face_basis == b + 1) ? 1 : 0;

            for (int i = 0; i < basis_dimension; i++) {
                g->minimals_basis[b].cycles[i].edges_ids =
                    calloc(g->nb_edges, sizeof(uint32_t));
                if (!g->minimals_basis[b].cycles[i].edges_ids) {
                    perror("load_graph: edges_ids"); exit(1);
                }

                char line[2048];
                if (!fgets(line, sizeof(line), file)) break;

                int count = 0;
                char *tok = strtok(line, " \t\r\n");
                while (tok) {
                    int eid = atoi(tok);
                    if (eid >= 0 && eid < g->nb_edges) {
                        g->minimals_basis[b].cycles[i].edges_ids[eid] = 1;
                        count++;
                    }
                    tok = strtok(NULL, " \t\r\n");
                }
                g->minimals_basis[b].cycles[i].length = count;
            }
        }
    }

    fprintf(stdout, "Load done.\n"); fflush(stdout);
    fclose(file);
}

/**
 * @brief Add a vertex at the given coordinates.
 *
 * Appends a new Vertex and initialises its (empty) neighbour list.  Doubles
 * array capacities when needed.
 *
 * @note Calls invalidate_bases() because adding a vertex changes V, which
 *       changes the cycle-space dimension D = E − V + 1, rendering all
 *       existing cycle basis vectors invalid.
 *
 * @param g  Graph to modify.
 * @param x  X coordinate (used for visualisation and face detection).
 * @param y  Y coordinate.
 */
void create_vertex(Graph *g, const double x, const double y) {
    invalidate_bases(g);

    if (g->nb_vertex == g->capacity_vertices) {
        const int old_cap = g->capacity_vertices;
        g->capacity_vertices *= 2;
        const int new_cap = g->capacity_vertices;

        Vertex *tv = realloc(g->vertices, new_cap * sizeof(Vertex));
        if (!tv) { perror("create_vertex: vertices"); exit(1); }
        g->vertices = tv;

        Neighbor_list *tn = realloc(g->neighbors, new_cap * sizeof(Neighbor_list));
        if (!tn) { perror("create_vertex: neighbors"); exit(1); }
        g->neighbors = tn;

        for (int i = old_cap; i < new_cap; i++) {
            g->neighbors[i].neighbors = NULL;
            g->neighbors[i].angles    = NULL;
            g->neighbors[i].count     = 0;
            g->neighbors[i].capacity  = 0;
        }

        Path *tf = realloc(g->faces,
            (2 + g->capacity_edges + new_cap) * sizeof(Path));
        if (!tf) { perror("create_vertex: faces"); exit(1); }
        g->faces = tf;
    }

    const int id = g->nb_vertex;
    g->vertices[id] = (Vertex){x, y, id, 0, 0};
    g->nb_vertex++;
}

/**
 * @brief Add an undirected edge between two existing vertices.
 *
 * Inserts the edge into the edges array and updates both endpoints' neighbour
 * lists, maintaining sorted polar-angle order (required by find_faces()).
 * Doubles array capacities when needed.
 *
 * @note Calls invalidate_bases() because adding an edge increases nb_edges,
 *       making all existing edges_ids vectors (sized to the old nb_edges) invalid.
 *
 * @param g     Graph to modify.
 * @param v1_id ID of the first endpoint (must be valid and not deleted).
 * @param v2_id ID of the second endpoint (must be valid and not deleted).
 */
void create_edge(Graph *g, const int v1_id, const int v2_id) {
    invalidate_bases(g);

    if (g->nb_edges == g->capacity_edges) {
        g->capacity_edges *= 2;

        Edge *te = realloc(g->edges, g->capacity_edges * sizeof(Edge));
        if (!te) { perror("create_edge: edges"); exit(1); }
        g->edges = te;

        Path *tf = realloc(g->faces,
            (2 + g->capacity_edges + g->capacity_vertices) * sizeof(Path));
        if (!tf) { perror("create_edge: faces"); exit(1); }
        g->faces = tf;
    }

    /* Update sorted neighbour lists for both endpoints */
    const int ids[2] = {v1_id, v2_id};
    for (int i = 0; i < 2; i++) {
        const double angle = atan2l(
            g->vertices[ids[(i + 1) % 2]].y - g->vertices[ids[i]].y,
            g->vertices[ids[(i + 1) % 2]].x - g->vertices[ids[i]].x);

        Neighbor_list *nl = &g->neighbors[ids[i]];
        if (nl->count == nl->capacity) {
            nl->capacity = nl->capacity == 0 ? 4 : nl->capacity * 2;
            int *tn = realloc(nl->neighbors, nl->capacity * sizeof(int));
            if (!tn) { perror("create_edge: nl->neighbors"); exit(1); }
            nl->neighbors = tn;
            double *ta = realloc(nl->angles, nl->capacity * sizeof(double));
            if (!ta) { perror("create_edge: nl->angles"); exit(1); }
            nl->angles = ta;
        }

        /* Insertion sort to maintain ascending angle order */
        int pos = nl->count;
        while (pos > 0 && nl->angles[pos - 1] > angle) {
            nl->neighbors[pos] = nl->neighbors[pos - 1];
            nl->angles[pos]    = nl->angles[pos - 1];
            pos--;
        }
        nl->neighbors[pos] = ids[(i + 1) % 2];
        nl->angles[pos] = angle;
        nl->count++;
    }

    g->edges[g->nb_edges] =
        (Edge){v1_id, v2_id, g->nb_edges, g->nb_edges, 0, 0, {-1, -1}};
    g->nb_edges++;
}

/**
 * @brief Mark an edge as deleted and remove it from both endpoints' neighbour lists.
 *
 * The edge slot is kept in the array (its index is preserved) but its deleted
 * flag is set to 1 so all algorithms skip it.
 *
 * @note Calls invalidate_bases() because existing cycle basis vectors may
 *       reference this edge ID, making the stored bases semantically invalid.
 *
 * @param g    Graph to modify.
 * @param e_id Index of the edge to delete.
 */
void delete_edge(Graph *g, const int e_id) {
    invalidate_bases(g);

    g->edges[e_id].deleted = 1;

    const int u = g->edges[e_id].u;
    const int v = g->edges[e_id].v;

    /* Remove v from u's neighbour list (swap-with-last, O(1)) */
    for (int i = 0; i < g->neighbors[u].count; i++) {
        if (g->neighbors[u].neighbors[i] == v) {
            const int last = --g->neighbors[u].count;
            g->neighbors[u].neighbors[i] = g->neighbors[u].neighbors[last];
            g->neighbors[u].angles[i]    = g->neighbors[u].angles[last];
            break;
        }
    }
    /* Remove u from v's neighbour list */
    for (int i = 0; i < g->neighbors[v].count; i++) {
        if (g->neighbors[v].neighbors[i] == u) {
            const int last = --g->neighbors[v].count;
            g->neighbors[v].neighbors[i] = g->neighbors[v].neighbors[last];
            g->neighbors[v].angles[i]    = g->neighbors[v].angles[last];
            break;
        }
    }
}

/**
 * @brief Mark a vertex as deleted together with all its incident edges.
 *
 * invalidate_bases() is called implicitly through each delete_edge() call.
 *
 * @param g    Graph to modify.
 * @param v_id ID of the vertex to delete.
 */
void delete_vertex(Graph *g, const int v_id) {
    for (int e = 0; e < g->nb_edges; e++) {
        if (g->edges[e].deleted) continue;
        if (g->edges[e].u == v_id || g->edges[e].v == v_id)
            delete_edge(g, e);
    }
    g->vertices[v_id].deleted = 1;
}

/**
 * @brief Subdivide an edge by inserting @p number_vertex_to_add intermediate vertices.
 *
 * Deletes edge (u, v) and replaces it with the path
 * u − w₁ − … − wₖ − v, where each wᵢ lies at equal intervals along (u, v).
 *
 * invalidate_bases() is called implicitly through delete_edge() and
 * create_edge() / create_vertex(); all calls are idempotent.
 *
 * @param g                    Graph to modify.
 * @param edge_id              Index of the edge to split.
 * @param number_vertex_to_add Number of intermediate vertices to insert (k ≥ 1).
 */
void split_edge(Graph *g, const int edge_id, const int number_vertex_to_add) {
    const int    u   = g->edges[edge_id].u;
    const int    v   = g->edges[edge_id].v;
    const double ux  = g->vertices[u].x, uy = g->vertices[u].y;
    const double vx  = g->vertices[v].x, vy = g->vertices[v].y;

    delete_edge(g, edge_id);

    int prev = u;
    for (int i = 1; i <= number_vertex_to_add; i++) {
        const double t = (double)i / (number_vertex_to_add + 1);
        create_vertex(g, ux + t * (vx - ux), uy + t * (vy - uy));
        const int nv = g->nb_vertex - 1;
        create_edge(g, prev, nv);
        prev = nv;
    }
    create_edge(g, prev, v);
}

/**
 * @brief Subdivide multiple edges, distributing intermediate vertices evenly.
 *
 * The @p number_vertex_to_add vertices are distributed across all edges as
 * evenly as possible; the first (number_vertex_to_add % number_edge_to_split)
 * edges each receive one extra vertex.
 *
 * @param g                    Graph to modify.
 * @param edge_ids             Array of edge indices to split.
 * @param number_edge_to_split Number of edges to split.
 * @param number_vertex_to_add Total number of intermediate vertices to insert.
 */
void split_edges(Graph *g, const int *edge_ids,
                 const int number_edge_to_split, const int number_vertex_to_add) {
    if (number_edge_to_split <= 0) return;
    const int base  = number_vertex_to_add / number_edge_to_split;
    const int extra = number_vertex_to_add % number_edge_to_split;
    for (int i = 0; i < number_edge_to_split; i++)
        split_edge(g, edge_ids[i], base + (i < extra ? 1 : 0));
}

/**
 * @brief Sort a vertex's neighbour list by polar angle in ascending order.
 *
 * Uses bubble sort on the cached angles array.  O(deg²) but degrees are small
 * for planar graphs.
 *
 * @param g    Graph.
 * @param v_id Vertex whose neighbour list is to be sorted.
 */
void sort_neighbors_by_atan(const Graph *g, const int v_id) {
    const Neighbor_list *nl = &g->neighbors[v_id];
    for (int i = 0; i < nl->count - 1; i++) {
        for (int j = 0; j < nl->count - i - 1; j++) {
            if (nl->angles[j] > nl->angles[j + 1]) {
                const int tn = nl->neighbors[j];
                const double ta = nl->angles[j];
                nl->neighbors[j] = nl->neighbors[j + 1];
                nl->angles[j] = nl->angles[j + 1];
                nl->neighbors[j + 1] = tn;
                nl->angles[j + 1] = ta;
            }
        }
    }
}

/**
 * @brief Move a vertex to new coordinates and update all geometry-derived data.
 *
 * Recalculates the stored polar angles for the moved vertex and all its
 * neighbours, re-sorts neighbour lists, then recomputes the face decomposition.
 *
 * @note Moving a vertex does NOT change the graph topology (V, E, and adjacency
 *       are unchanged), so cycle bases remain valid and invalidate_bases() is
 *       NOT called.  Only face data changes because the planar embedding may
 *       change when vertex positions change.
 *
 * @param g     Graph to modify.
 * @param v     ID of the vertex to move.
 * @param new_x New X coordinate.
 * @param new_y New Y coordinate.
 */
void move_vertex(Graph *g, const int v, const double new_x, const double new_y) {
    g->vertices[v].x = new_x;
    g->vertices[v].y = new_y;

    /*
     * Recompute stored angles BEFORE sorting.
     * sort_neighbors_by_atan() sorts by cached values; calling it with stale
     * angles would produce an incorrect order and corrupt face detection.
     */
    Neighbor_list *nlv = &g->neighbors[v];
    for (int i = 0; i < nlv->count; i++) {
        const int nb = nlv->neighbors[i];
        nlv->angles[i] = atan2(g->vertices[nb].y - new_y,
                               g->vertices[nb].x - new_x);

        Neighbor_list *nl_nb = &g->neighbors[nb];
        for (int j = 0; j < nl_nb->count; j++) {
            if (nl_nb->neighbors[j] == v) {
                nl_nb->angles[j] = atan2(new_y - g->vertices[nb].y,
                                         new_x - g->vertices[nb].x);
                break;
            }
        }
    }

    sort_neighbors_by_atan(g, v);
    for (int i = 0; i < g->neighbors[v].count; i++)
        sort_neighbors_by_atan(g, g->neighbors[v].neighbors[i]);

    find_faces(g);
}

/**
 * @brief Allocate and initialise the n×n matrices for Horton's algorithm.
 *
 * Allocates (or reallocates) three n×n integer matrices:
 *  - edge_indices[u*n+v]: index of edge (u,v) in g->edges, or -1,
 *  - predecessors[u*n+v]: predecessor of v on the shortest path from u,
 *  - distances[u*n+v]:    BFS distance from u to v.
 *
 * All entries initialised to -1; edge_indices filled for non-deleted edges.
 *
 * @param graph Graph to prepare.
 */
void prepare_graph_matrices(Graph *graph) {
    if (!graph->edges || graph->nb_vertex == 0) {
        printf("prepare_graph_matrices: empty graph, skipping.\n");
        return;
    }

    const int n    = graph->nb_vertex;
    const int size = n * n;

    free(graph->edge_indices);
    free(graph->predecessors);
    free(graph->distances);

    graph->edge_indices = malloc(size * sizeof(int));
    graph->predecessors = malloc(size * sizeof(int));
    graph->distances    = malloc(size * sizeof(int));
    if (!graph->edge_indices || !graph->predecessors || !graph->distances) {
        perror("prepare_graph_matrices: malloc"); exit(1);
    }

    for (int i = 0; i < size; i++) {
        graph->edge_indices[i] = -1;
        graph->predecessors[i] = -1;
        graph->distances[i]    = -1;
    }

    for (int i = 0; i < graph->nb_edges; i++) {
        if (graph->edges[i].deleted) continue;
        const int u = (int)graph->edges[i].u;
        const int v = (int)graph->edges[i].v;
        graph->edge_indices[u * n + v] = i;
        graph->edge_indices[v * n + u] = i;
    }
}

/**
 * @brief Compute all-pairs BFS shortest paths and predecessor tables.
 *
 * Runs one BFS per non-deleted source vertex.  Then reconstructs the
 * predecessor[u*n+v] entry for every reachable pair (u, v) by walking
 * backward along the distance array.  When several predecessors yield the
 * same distance, the one with the smallest vertex label is chosen
 * (deterministic tie-breaking across permutations).
 *
 * @param g Graph whose matrices (allocated by prepare_graph_matrices) are filled.
 */
void compute_all_shortest_paths(const Graph *g) {
    int *queue = malloc(g->nb_vertex * sizeof(int));
    if (!queue) { perror("compute_all_shortest_paths: queue"); exit(1); }

    for (int s = 0; s < g->nb_vertex; s++) {
        if (g->vertices[s].deleted) continue;

        g->distances[s * g->nb_vertex + s] = 0;
        int front = 0, rear = 1;
        queue[front] = s;

        while (front < rear) {
            const int u = queue[front++];
            for (int i = 0; i < g->neighbors[u].count; i++) {
                const int v = g->neighbors[u].neighbors[i];
                if (g->vertices[v].deleted) continue;
                if (g->distances[s * g->nb_vertex + v] == -1) {
                    g->distances[s * g->nb_vertex + v] =
                        g->distances[s * g->nb_vertex + u] + 1;
                    queue[rear++] = v;
                }
            }
        }
    }
    free(queue);

    for (int u = 0; u < g->nb_vertex; u++) {
        if (g->vertices[u].deleted) continue;
        for (int v = 0; v < g->nb_vertex; v++) {
            if (g->vertices[v].deleted) continue;
            int dist = g->distances[u * g->nb_vertex + v];
            if (dist <= 0) continue;

            int current = v;
            while (dist > 0) {
                int best = -1;
                for (int i = 0; i < g->neighbors[current].count; i++) {
                    const int nb = g->neighbors[current].neighbors[i];
                    if (g->vertices[nb].deleted) continue;
                    if (g->distances[u * g->nb_vertex + nb] == dist - 1) {
                        if (best == -1 ||
                            g->vertices[nb].label < g->vertices[best].label)
                            best = nb;
                    }
                }
                if (best == -1) {
                    fprintf(stderr,
                        "compute_all_shortest_paths: no predecessor for vertex %d "
                        "on path %d→%d\n", current, u, v);
                    exit(1);
                }
                g->predecessors[u * g->nb_vertex + current] = best;
                current = best;
                dist--;
            }
        }
    }
}

/**
 * @brief Build the shortest path from vertex u to vertex v as a Path object.
 *
 * Follows the predecessor chain stored in g->predecessors.  The returned Path
 * owns its edges_ids, vertices_ids, and edges_labels arrays; the caller must
 * free them.
 *
 * @param g  Graph with precomputed shortest paths.
 * @param u  Source vertex.
 * @param v  Destination vertex.
 * @return   A Path from u to v.  Exits on error.
 */
Path create_path(const Graph *g, const int u, const int v) {
    if (!g->predecessors || !g->edges) {
        perror("create_path: shortest paths not computed"); exit(1);
    }

    Path p = {0};
    p.edges_ids = calloc(g->nb_edges, sizeof(uint32_t));
    if (!p.edges_ids) { perror("create_path: edges_ids"); exit(1); }

    p.vertices_ids = calloc(g->nb_vertex, sizeof(uint32_t));
    if (!p.vertices_ids) {
        free(p.edges_ids); perror("create_path: vertices_ids"); exit(1);
    }

    p.edges_labels = calloc(g->nb_edges, sizeof(uint32_t));
    if (!p.edges_labels) {
        free(p.edges_ids); free(p.vertices_ids);
        perror("create_path: edges_labels"); exit(1);
    }

    p.vertices_ids[u] = 1;
    p.vertices_ids[v] = 1;

    int current = v;
    while (current != u) {
        const int pred = g->predecessors[u * g->nb_vertex + current];
        if (pred == -1) {
            fprintf(stderr, "create_path: missing predecessor for %d on %d→%d\n",
                    current, u, v);
            free(p.edges_ids); free(p.vertices_ids); free(p.edges_labels);
            exit(1);
        }
        const int eidx = g->edge_indices[current * g->nb_vertex + pred];
        if (eidx == -1) {
            fprintf(stderr, "create_path: no edge between %d and %d\n", current, pred);
            free(p.edges_ids); free(p.vertices_ids); free(p.edges_labels);
            exit(1);
        }
        p.edges_ids[eidx] = 1;
        p.vertices_ids[g->edges[eidx].u] = 1;
        p.vertices_ids[g->edges[eidx].v] = 1;
        p.edges_labels[g->edges[eidx].label] = 1;
        p.length++;
        current = pred;
    }
    return p;
}

void find_faces(Graph *g) {
    // 1. Nettoyage de la mémoire si des faces existent déjà
    if (g->faces) {
        for (int i = 0; i < g->nb_faces; i++) {
            free(g->faces[i].edges_ids);
            free(g->faces[i].edges_labels);
            free(g->faces[i].vertices_ids);
        }
    }

    // 2. Réinitialisation des propriétés des arêtes
    for (int e = 0; e < g->nb_edges; e++) {
        g->edges[e].face_id[0] = -1;
        g->edges[e].face_id[1] = -1;
        g->edges[e].is_outer   = 0;
    }

    const int n = g->nb_vertex;
    const int m = g->nb_edges;
    if (!g->edge_indices) prepare_graph_matrices(g);

    // Tableau pour marquer les demi-arêtes (u -> p) visitées
    int *visited = calloc(n * n, sizeof(int));
    if (!visited) { perror("find_faces: visited"); exit(1); }

    g->nb_faces = 0;
    g->outer_face = -1;
    double max_area = -1e18; // On cherche l'aire positive maximale pour la face externe

    for (int v = 0; v < n; v++) {
        if (g->vertices[v].deleted) continue;

        for (int i = 0; i < g->neighbors[v].count; i++) {
            const int nb = g->neighbors[v].neighbors[i];
            if (visited[v * n + nb]) continue;

            Path *face = &g->faces[g->nb_faces];
            face->edges_ids    = calloc(m, sizeof(uint32_t));
            face->vertices_ids = calloc(n, sizeof(uint32_t));
            face->edges_labels = calloc(m, sizeof(uint32_t));
            face->length       = 0;

            double area = 0.0;
            int u = v, p = nb;

            // Parcours du cycle de la face
            while (!visited[u * n + p]) {
                visited[u * n + p] = 1;
                face->length++;

                // Calcul de l'aire signée (Shoelace formula)
                // area = sum (x_u * y_p - x_p * y_u)
                area += (double)g->vertices[u].x * g->vertices[p].y;
                area -= (double)g->vertices[p].x * g->vertices[u].y;

                const int eid = g->edge_indices[u * n + p];
                if (eid >= 0) {
                    face->edges_ids[eid] = 1;
                    face->edges_labels[g->edges[eid].label] = 1;
                    if (g->edges[eid].face_id[0] == -1)
                        g->edges[eid].face_id[0] = g->nb_faces;
                    else
                        g->edges[eid].face_id[1] = g->nb_faces;
                }
                face->vertices_ids[u] = 1;

                // Sélection du prochain sommet en tournant à droite
                // On cherche l'indice de u dans la liste de p, puis on prend l'élément précédent
                const Neighbor_list *nl = &g->neighbors[p];
                int idx = -1;
                for (int j = 0; j < nl->count; j++) {
                    if (nl->neighbors[j] == u) { idx = j; break; }
                }

                // Le virage à droite (idx - 1) assure que les faces internes sont CW
                // et que la bordure externe est CCW.
                const int w = nl->neighbors[(idx - 1 + nl->count) % nl->count];
                u = p; p = w;
            }

            // Identification de la face externe
            // Dans ce système de virage à droite, seule l'aire de la face externe est positive
            if (area > max_area) {
                max_area = area;
                g->outer_face = g->nb_faces;
            }

            g->nb_faces++;
        }
    }

    free(visited);
}

/**
 * @brief Append a cycle to the graph's Horton cycle list.
 *
 * Doubles the horton_cycles capacity when the array is full.  Ownership of
 * the cycle's inner arrays is transferred to the graph.
 *
 * @param g  Graph to which the cycle is appended.
 * @param c  Cycle to add.
 */
void add_horton_cycle(Graph *g, const Path c) {
    if (g->nb_horton_cycles == g->capacity_horton_cycles) {
        const int nc = g->capacity_horton_cycles * 2;
        Path *tmp = realloc(g->horton_cycles, nc * sizeof(Path));
        if (!tmp) { perror("add_horton_cycle: realloc"); exit(1); }
        g->horton_cycles = tmp;
        g->capacity_horton_cycles = nc;
    }
    g->horton_cycles[g->nb_horton_cycles++] = c;
}

/**
 * @brief Generate all candidate cycles for Horton's algorithm.
 *
 * For each non-deleted root vertex r and each non-deleted edge (x, y) not
 * incident to r, computes the two shortest paths r→x and r→y.  If the paths
 * share only the root r as a common vertex, their symmetric difference plus
 * edge (x, y) forms a valid Horton cycle, which is added to g->horton_cycles.
 *
 * @param g Graph with precomputed shortest paths.
 */
void find_horton_cycles(Graph *g) {
    if (!g->predecessors || !g->edges) {
        printf("find_horton_cycles: shortest paths not computed, skipping.\n");
        return;
    }

    for (int v = 0; v < g->nb_vertex; v++) {
        if (g->vertices[v].deleted) continue;

        for (int e = 0; e < g->nb_edges; e++) {
            if (g->edges[e].deleted) continue;
            const int x = g->edges[e].u;
            const int y = g->edges[e].v;
            if (x == v || y == v) continue;

            const int px = g->predecessors[v * g->nb_vertex + x];
            const int py = g->predecessors[v * g->nb_vertex + y];
            if (px == -1 || py == -1) {
                fprintf(stderr,
                    "find_horton_cycles: no path from %d to %d or %d, "
                    "skipping edge (%d,%d)\n", v, x, y, x, y);
                continue;
            }

            const Path pv_x = create_path(g, v, x);
            const Path pv_y = create_path(g, v, y);

            uint32_t *expected = calloc(g->nb_vertex, sizeof(uint32_t));
            uint32_t *actual   = calloc(g->nb_vertex, sizeof(uint32_t));
            if (!expected || !actual) {
                free(expected); free(actual); perror("find_horton_cycles"); exit(1);
            }
            expected[v] = 1;
            for (int i = 0; i < g->nb_vertex; i++)
                actual[i] = pv_x.vertices_ids[i] & pv_y.vertices_ids[i];

            if (memcmp(expected, actual, g->nb_vertex * sizeof(uint32_t)) == 0) {
                uint32_t *ec = calloc(g->nb_edges,  sizeof(uint32_t));
                uint32_t *vc = calloc(g->nb_vertex, sizeof(uint32_t));
                uint32_t *el = calloc(g->nb_edges,  sizeof(uint32_t));
                if (!ec || !vc || !el) {
                    free(ec); free(vc); free(el);
                    free(expected); free(actual);
                    perror("find_horton_cycles: cycle alloc"); exit(1);
                }

                for (int i = 0; i < g->nb_edges; i++) {
                    ec[i] = pv_x.edges_ids[i] ^ pv_y.edges_ids[i];
                    if (ec[i]) el[g->edges[i].label] = 1;
                }
                for (int i = 0; i < g->nb_vertex; i++)
                    vc[i] = pv_x.vertices_ids[i] ^ pv_y.vertices_ids[i];

                ec[e] = 1;
                el[g->edges[e].label] = 1;

                add_horton_cycle(g, (Path){
                    .edges_ids    = ec,
                    .vertices_ids = vc,
                    .edges_labels = el,
                    .length       = pv_x.length + pv_y.length + 1
                });
            }

            free(pv_x.edges_ids); free(pv_x.vertices_ids); free(pv_x.edges_labels);
            free(pv_y.edges_ids); free(pv_y.vertices_ids); free(pv_y.edges_labels);
            free(expected); free(actual);
        }
    }
}

/** @brief Global used by qsort callbacks to communicate the edge count. */
static int nb_edges_for_qsort = 0;

/**
 * @brief Compare two cycles: first by length, then by edges_ids lexicographic order.
 *
 * Used by qsort() to sort Horton cycles before the greedy cover step.
 */
int compare_cycles(const void *a, const void *b) {
    const Path *c1 = a, *c2 = b;
    if (c1->length != c2->length) return c1->length - c2->length;
    return memcmp(c1->edges_ids, c2->edges_ids,
                  nb_edges_for_qsort * sizeof(uint32_t));
}

/**
 * @brief Sort the Horton cycle list by length (shortest first).
 * @param g Graph whose horton_cycles are to be sorted.
 */
void order_horton_cycles(const Graph *g) {
    nb_edges_for_qsort = g->nb_edges;
    if (g->nb_horton_cycles > 1)
        qsort(g->horton_cycles, g->nb_horton_cycles, sizeof(Path), compare_cycles);
}

/**
 * @brief Compare two cycles: first by length, then by edges_labels order.
 *
 * Used to explore different bases across edge-label permutations in
 * multiple_horton().
 */
int compare_cycles_permutations(const void *a, const void *b) {
    const Path *c1 = a, *c2 = b;
    if (c1->length != c2->length) return c1->length - c2->length;
    for (int i = 0; i < nb_edges_for_qsort; i++) {
        if (c1->edges_labels[i] != c2->edges_labels[i])
            return (int)c1->edges_labels[i] - (int)c2->edges_labels[i];
    }
    return 0;
}

/**
 * @brief Sort the Horton cycle list using the label-based comparator.
 * @param g Graph whose horton_cycles are to be sorted.
 */
void order_orton_cycles_permutations(const Graph *g) {
    nb_edges_for_qsort = g->nb_edges;
    if (g->nb_horton_cycles > 1)
        qsort(g->horton_cycles, g->nb_horton_cycles,
              sizeof(Path), compare_cycles_permutations);
}

/**
 * @brief Check whether two cycles have identical edge sets and length.
 *
 * @param c1       First cycle.
 * @param c2       Second cycle.
 * @param nb_edges Size of the edges_ids vectors.
 * @return 1 if identical, 0 otherwise.
 */
int are_same_cycles(const Path *c1, const Path *c2, const int nb_edges) {
    if (c1->length != c2->length) return 0;
    return memcmp(c1->edges_ids, c2->edges_ids, nb_edges * sizeof(uint32_t)) == 0;
}

/** @brief qsort wrapper for are_same_cycles using nb_edges_for_qsort. */
int are_same_cycles_qsort(const void *c1, const void *c2) {
    return are_same_cycles(c1, c2, nb_edges_for_qsort);
}

/**
 * @brief Check whether two minimal cycle bases are identical (as unordered sets).
 *
 * Sorts copies of both bases by cycle weight, then compares them element-wise.
 *
 * @param a               First basis.
 * @param b               Second basis.
 * @param nb_edges        Size of the edges_ids vectors.
 * @param basis_dimension Number of cycles in each basis.
 * @return 1 if the bases are identical, 0 otherwise.
 */
int are_same_minimal_bases(const void *a, const void *b,
                           const int nb_edges, const int basis_dimension) {
    const Minimal_basis *mb1 = a, *mb2 = b;

    Path *s1 = malloc(basis_dimension * sizeof(Path));
    Path *s2 = malloc(basis_dimension * sizeof(Path));
    if (!s1 || !s2) { free(s1); free(s2); perror("are_same_minimal_bases"); exit(1); }

    for (int i = 0; i < basis_dimension; i++) {
        s1[i].length    = mb1->cycles[i].length;
        s1[i].edges_ids = malloc(nb_edges * sizeof(uint32_t));
        s2[i].length    = mb2->cycles[i].length;
        s2[i].edges_ids = malloc(nb_edges * sizeof(uint32_t));
        if (!s1[i].edges_ids || !s2[i].edges_ids) {
            perror("are_same_minimal_bases: malloc"); exit(1);
        }
        memcpy(s1[i].edges_ids, mb1->cycles[i].edges_ids, nb_edges * sizeof(uint32_t));
        memcpy(s2[i].edges_ids, mb2->cycles[i].edges_ids, nb_edges * sizeof(uint32_t));
    }

    nb_edges_for_qsort = nb_edges;
    qsort(s1, basis_dimension, sizeof(Path), compare_cycles);
    qsort(s2, basis_dimension, sizeof(Path), compare_cycles);

    int result = 1;
    for (int i = 0; i < basis_dimension; i++) {
        if (!are_same_cycles(&s1[i], &s2[i], nb_edges)) { result = 0; break; }
    }

    for (int i = 0; i < basis_dimension; i++) {
        free(s1[i].edges_ids);
        free(s2[i].edges_ids);
    }
    free(s1);
    free(s2);
    return result;
}

/**
 * @brief Check whether a candidate basis already appears in the graph's list.
 *
 * @param g         Graph containing the known bases.
 * @param candidate Candidate basis to search for.
 * @return 1 if the candidate duplicates an existing basis, 0 otherwise.
 */
int already_exists_basis(const Graph *g, const Minimal_basis candidate) {
    if (!g->minimals_basis || g->nb_minimal_bases == 0) return 0;
    for (int b = 0; b < g->nb_minimal_bases; b++) {
        if (are_same_minimal_bases(&g->minimals_basis[b], &candidate,
                                   g->nb_edges, g->basis_dimension))
            return 1;
    }
    return 0;
}

/**
 * @brief Test whether a minimal basis equals the set of interior faces.
 *
 * Compares the candidate basis @p mb (sorted by edge set) to the interior
 * faces of the graph (all faces except the outer face, also sorted).
 *
 * @param g  Graph with a precomputed face decomposition.
 * @param mb Candidate minimal basis to test.
 * @return 1 if @p mb matches the face basis, 0 otherwise.
 */
int mb_is_faces(Graph *g, const Minimal_basis *mb) {

    if (!g->faces || g->nb_faces == 0) find_faces(g);

    Path *s1 = malloc(g->basis_dimension * sizeof(Path));
    Path *s2 = malloc(g->nb_faces * sizeof(Path));

    if (!s1 || !s2) {
        free(s1); free(s2);
        perror("mb_is_faces: allocation failure");
        exit(1);
    }

    for (int i = 0; i < g->basis_dimension; i++) {
        s1[i].length = mb->cycles[i].length;
        s1[i].edges_ids = malloc(g->nb_edges * sizeof(uint32_t));
        if (!s1[i].edges_ids) { perror("malloc s1 edges"); exit(1); }
        memcpy(s1[i].edges_ids, mb->cycles[i].edges_ids, g->nb_edges * sizeof(uint32_t));
    }

    for (int i = 0; i < g->nb_faces; i++) {
        s2[i].length = g->faces[i].length;
        s2[i].edges_ids = malloc(g->nb_edges * sizeof(uint32_t));
        if (!s2[i].edges_ids) { perror("malloc s2 edges"); exit(1); }
        memcpy(s2[i].edges_ids, g->faces[i].edges_ids, g->nb_edges * sizeof(uint32_t));
    }

    nb_edges_for_qsort = g->nb_edges;
    qsort(s1, g->basis_dimension, sizeof(Path), compare_cycles);
    qsort(s2, g->nb_faces, sizeof(Path), compare_cycles);

    int result = 1;
    int missed_idx = -1;
    int i_s1 = 0;
    int i_s2 = 0;

    while (i_s1 < g->basis_dimension && i_s2 < g->nb_faces) {
        if (are_same_cycles(&s1[i_s1], &s2[i_s2], g->nb_edges)) {
            i_s1++;
            i_s2++;
        } else {
            if (missed_idx != -1) {
                result = 0;
                break;
            }
            missed_idx = i_s2;
            i_s2++;
        }
    }

    if (result && missed_idx == -1 && i_s2 < g->nb_faces) {
        missed_idx = i_s2;
    }

    int final_return = 0;
    if (result && i_s1 == g->basis_dimension) {
        if (are_same_cycles(&s2[missed_idx], &g->faces[g->outer_face], g->nb_edges)) {
            final_return = 1;
        } else {
            final_return = 2;
        }
    }

    for (int i = 0; i < g->nb_faces; i++) {
        if (i < g->basis_dimension) {
            free(s1[i].edges_ids);
        }
        free(s2[i].edges_ids);
    }
    free(s1);
    free(s2);

    return final_return;
}

/**
 * @brief Extract one minimal cycle basis from the sorted Horton cycle list.
 *
 * Iterates over g->horton_cycles (sorted by length) and greedily selects each
 * cycle that is linearly independent of the already-chosen ones, using Gaussian
 * elimination over GF(2) on the binary edge vectors.
 *
 * If the dimension of the new basis differs from the stored one (which should
 * not happen if invalidate_bases() was called before every topology mutation),
 * all stale bases are cleared and the new basis replaces them.
 *
 * The is_faces flag is set when the computed basis matches the interior face set
 * (Euler: F = E − V + 2, so nb_faces = basis_dimension + 1).
 *
 * @param g Graph.  g->horton_cycles must be non-empty and sorted.
 */
void greedy_cycle_cover(Graph *g) {
    if (!g->horton_cycles || g->nb_horton_cycles == 0) {
        g->basis_dimension = 0;
        printf("greedy_cycle_cover: no Horton cycles, skipping.\n");
        return;
    }

    uint32_t **basis = calloc(g->nb_edges, sizeof(uint32_t *));
    Path      *mcb   = calloc(g->nb_horton_cycles, sizeof(Path));
    if (!basis || !mcb) { perror("greedy_cycle_cover: alloc"); exit(1); }

    int mcb_count = 0;

    for (int i = 0; i < g->nb_horton_cycles; i++) {
        const Path cand = g->horton_cycles[i];

        uint32_t *vec = malloc(g->nb_edges * sizeof(uint32_t));
        if (!vec) { perror("greedy_cycle_cover: vec"); exit(1); }
        memcpy(vec, cand.edges_ids, g->nb_edges * sizeof(uint32_t));

        int independent = 0;
        for (int e = 0; e < g->nb_edges; e++) {
            if (!vec[e]) continue;
            if (!basis[e]) { basis[e] = vec; independent = 1; break; }
            for (int k = 0; k < g->nb_edges; k++) vec[k] ^= basis[e][k];
        }

        if (independent) {
            mcb[mcb_count].length = cand.length;
            mcb[mcb_count].edges_ids = malloc(g->nb_edges * sizeof(uint32_t));
            mcb[mcb_count].edges_labels = malloc(g->nb_edges * sizeof(uint32_t));
            if (!mcb[mcb_count].edges_ids || !mcb[mcb_count].edges_labels) {
                perror("greedy_cycle_cover: mcb entry"); exit(1);
            }
            memcpy(mcb[mcb_count].edges_ids,    cand.edges_ids,    g->nb_edges * sizeof(uint32_t));
            memcpy(mcb[mcb_count].edges_labels, cand.edges_labels, g->nb_edges * sizeof(uint32_t));
            mcb_count++;
        } else {
            free(vec);
        }
    }

    for (int e = 0; e < g->nb_edges; e++) free(basis[e]);
    free(basis);

    if (mcb_count > 0) {
        Path *ex = realloc(mcb, mcb_count * sizeof(Path));
        if (ex) mcb = ex;
        else { perror("greedy_cycle_cover: realloc mcb"); exit(1); }
    }

    Minimal_basis mb = {.cycles = mcb};

    if (g->nb_minimal_bases != 0) {
        if (mcb_count != g->basis_dimension) {
            /*
             * Safety net: dimension mismatch means topology changed without
             * triggering invalidate_bases().  Discard all stale bases.
             * This branch should never be reached in normal usage.
             */
            fprintf(stderr,
                "greedy_cycle_cover: dimension changed (%d → %d), "
                "clearing stale bases.\n", g->basis_dimension, mcb_count);
            for (int b = 0; b < g->nb_minimal_bases; b++) {
                for (int i = 0; i < g->basis_dimension; i++) {
                    free(g->minimals_basis[b].cycles[i].edges_ids);
                    free(g->minimals_basis[b].cycles[i].edges_labels);
                    free(g->minimals_basis[b].cycles[i].vertices_ids);
                }
                free(g->minimals_basis[b].cycles);
            }
            free(g->minimals_basis);
            g->minimals_basis   = NULL;
            g->nb_minimal_bases = 0;
            g->face_basis       = -1;
        } else if (already_exists_basis(g, mb)) {
            for (int k = 0; k < mcb_count; k++) {
                free(mcb[k].edges_ids);
                free(mcb[k].edges_labels);
            }
            free(mcb);
            return;
        }
    }

    g->basis_dimension = mcb_count;

    /*
     * Check if this basis equals the set of interior faces.
     * nb_faces includes the outer face, so the interior count is nb_faces - 1,
     * which equals the cycle-space dimension D for a connected planar graph.
     */
    mb.is_faces = 0;
    if (g->faces && g->nb_faces == mcb_count + 1) {
        const Minimal_basis face_mb = {.cycles = g->faces};
        const int is_face = mb_is_faces(g, &face_mb);
        if (is_face == 1) {
            mb.is_faces = 1;
            g->face_basis = g->nb_minimal_bases + 1; /* 1-based */
        } else if (is_face == 2) {
            mb.is_faces_outer = 1;
            g->face_basis_outer = g->nb_minimal_bases + 1;
        }
    }

    g->nb_minimal_bases++;
    Minimal_basis *tmp = realloc(g->minimals_basis,
                                 g->nb_minimal_bases * sizeof(Minimal_basis));
    if (!tmp) { perror("greedy_cycle_cover: realloc minimals_basis"); exit(1); }
    g->minimals_basis = tmp;
    g->minimals_basis[g->nb_minimal_bases - 1] = mb;
}

/**
 * @brief Free the Horton cycle list without touching any other graph state.
 *
 * Called at the start of each Horton run so stale cycles from the previous
 * invocation do not pollute the new enumeration.
 *
 * @param g Graph whose Horton cycle list is to be cleared.
 */
void reset_graph_results(Graph *g) {
    if (g->horton_cycles) {
        for (int i = 0; i < g->nb_horton_cycles; i++) {
            free(g->horton_cycles[i].edges_ids);
            free(g->horton_cycles[i].vertices_ids);
            free(g->horton_cycles[i].edges_labels);
        }
        g->nb_horton_cycles = 0;
    }
}

/**
 * @brief Run one iteration of Horton's minimal cycle basis algorithm.
 *
 * Pipeline:
 *  1. reset_graph_results()         — clear previous Horton cycle list,
 *  2. prepare_graph_matrices()      — (re-)allocate edge_indices / BFS matrices,
 *  3. compute_all_shortest_paths()  — BFS from every non-deleted vertex,
 *  4. find_faces() if nb_faces == 0 — face decomposition for face-basis detection,
 *  5. find_horton_cycles()          — enumerate candidate cycles,
 *  6. order_orton_cycles_permutations() — sort by length + label order,
 *  7. greedy_cycle_cover()          — extract one minimal basis.
 *
 * @param g Graph on which to run the algorithm.
 */
void horton(Graph *g) {
    reset_graph_results(g);
    prepare_graph_matrices(g);
    compute_all_shortest_paths(g);
    find_faces(g);
    find_horton_cycles(g);
    order_orton_cycles_permutations(g);
    greedy_cycle_cover(g);
}

/**
 * @brief Run Horton's algorithm multiple times with different edge-label permutations.
 *
 * By permuting edge labels between runs, the greedy cover selects different
 * cycles among those of equal length, potentially discovering new minimal bases.
 * Duplicate bases are detected and skipped by greedy_cycle_cover().
 *
 * @param g               Graph on which to run.
 * @param inv             Caller-allocated inversion table buffer (size: g->nb_edges).
 *                        Contents are overwritten on every iteration.
 * @param max_permutations Number of permuted runs to execute.
 */
void multiple_horton(Graph *g, int *inv, const int max_permutations) {
    if (!inv || !g->edges) {
        printf("multiple_horton: missing inversion table or edges, skipping.\n");
        return;
    }

    int *perm = calloc(g->nb_edges, sizeof(int));
    if (!perm) { perror("multiple_horton: perm"); exit(1); }

    generate_random_inversion_table(inv, g->nb_edges);
    inversion_to_permutation(inv, perm, g->nb_edges);
    for (int i = 0; i < g->nb_edges; i++) g->edges[i].label = perm[i];
    horton(g);

    for (int cnt = 1; cnt < max_permutations; cnt++) {
        memset(perm, 0, g->nb_edges * sizeof(int));
        generate_random_inversion_table(inv, g->nb_edges);
        inversion_to_permutation(inv, perm, g->nb_edges);
        for (int i = 0; i < g->nb_edges; i++) g->edges[i].label = perm[i];
        printf("%d / %d\n", cnt, max_permutations);
        horton(g);
    }

    free(perm);
}