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
 * Must be called whenever the topology changes (vertices or edges added/removed). Resets the
 * corresponding counters so that algorithms recompute from scratch on the new topology.
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

    /* Free all face-related structures */
    for (int i = 0; i < g->nb_faces; i++) {
        free(g->faces[i].edges_ids);
        free(g->faces[i].edges_labels);
        free(g->faces[i].vertices_ids);
        g->faces[i].edges_ids = NULL;
        g->faces[i].edges_labels = NULL;
        g->faces[i].vertices_ids = NULL;
    }
    g->nb_faces = 0;
    g->outer_face = -1;
    g->face_basis = -1;
    g->no_face_basis_possible = 0;
    g->nb_face_basis_outer = 0;
    free(g->face_basis_outer);
    g->face_basis_outer = NULL;
}

/**
 * @brief Allocate and initialize an empty graph.
 *
 * @return Pointer to the newly created graph. Exits on allocation failure.
 */
Graph *create_graph() {
    Graph *graph = calloc(1, sizeof(Graph));
    if (!graph) { perror("create_graph: graph"); exit(1); }

    graph->capacity_vertices = 4;
    graph->vertices = malloc(graph->capacity_vertices * sizeof(Vertex));
    if (!graph->vertices) { free(graph); perror("create_graph: vertices"); exit(1); }

    graph->neighbors = malloc(graph->capacity_vertices * sizeof(Neighbor_list));
    if (!graph->neighbors) {
        free(graph->vertices);
        free(graph);
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

    graph->face_basis_outer = NULL;

    graph->face_basis = -1;
    graph->outer_face = -1;
    graph->nb_face_basis_outer = 0;
    graph->no_face_basis_possible = 0;
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
 * @param g Graph to destroy.
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
    free(g->face_basis_outer);

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
 * V E B D face_basis outer_face_basis_1 ... outer_face_basis_n
 * vertex_id x_coordinate y_coordinate (V vertex lines)
 * edge_id id_vertex_u id_vertex_v (E edge lines)
 * edge_id ... (B * D cycle lines: space-separated edge IDs per cycle)
 *
 * V is the number of vertices, E the number of edges, B the number of different minimal basis,
 * D the dimension of the basis, face_basis the index of the basis representing the internal faces
 * and outer_face_basis_i the index of the basis representing the outer face and all the internal
 * ones except one.
 *
 * @param g Graph to save.
 * @param filename Output file path.
 */
void save_graph(const Graph *g, const char *filename) {

    FILE *file = fopen(filename, "w");
    if (!file) { perror("save_graph: fopen"); exit(1); }
    if (!g || !g->vertices || !g->edges) {
        fclose(file); perror("save_graph: invalid graph"); exit(1);
    }

    /* Build remapping tables */
    int *vertices_id_map = calloc(g->nb_vertex, sizeof(int));
    int *edges_id_map = calloc(g->nb_edges,  sizeof(int));
    if (!vertices_id_map || !edges_id_map) { perror("save_graph: maps"); exit(1); }

    int nb_v = 0, nb_e = 0;
    for (int i = 0; i < g->nb_vertex; i++)
        vertices_id_map[i] = g->vertices[i].deleted ? -1 : nb_v++;
    for (int i = 0; i < g->nb_edges; i++)
        edges_id_map[i] = g->edges[i].deleted ? -1 : nb_e++;

    if (g->minimals_basis) {
        fprintf(file, "%d %d %d %d %d %d",
                nb_v, nb_e, g->nb_minimal_bases, g->basis_dimension,
                g->face_basis, g->nb_face_basis_outer);
        for (int i = 0; i < g->nb_face_basis_outer; i++)
            fprintf(file, " %d", g->face_basis_outer[i]);
        fprintf(file, "\n");
    } else {
        fprintf(file, "%d %d 0 0 -1 0\n", nb_v, nb_e);
    }

    for (int i = 0; i < g->nb_vertex; i++) {
        if (g->vertices[i].deleted) continue;
        fprintf(file, "%d %lf %lf\n", vertices_id_map[i], g->vertices[i].x, g->vertices[i].y);
    }

    for (int i = 0; i < g->nb_edges; i++) {
        if (g->edges[i].deleted) continue;
        const int nu = vertices_id_map[g->edges[i].u];
        const int nv = vertices_id_map[g->edges[i].v];
        if (nu < 0 || nv < 0) continue;
        fprintf(file, "%d %d %d\n", edges_id_map[i], nu, nv);
    }

    if (g->minimals_basis) {
        for (int b = 0; b < g->nb_minimal_bases; b++) {
            for (int i = 0; i < g->basis_dimension; i++) {
                for (int e = 0; e < g->nb_edges; e++) {
                    if (g->edges[e].deleted) continue;
                    if (g->minimals_basis[b].cycles[i].edges_ids[e] == 1)
                        fprintf(file, "%d ", edges_id_map[e]);
                }
                fprintf(file, "\n");
            }
        }
    }

    free(vertices_id_map);
    free(edges_id_map);
    fclose(file);
}

/**
 * @brief Load a graph from a file produced by save_graph().
 *
 * g must be a newly created and empty graph. The format should be the same that presented in
 * save_graph().
 *
 * @param g Destination graph (must be empty).
 * @param filename File to read.
 */
void load_graph(Graph *g, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) { perror("load_graph: fopen"); exit(1); }

    int nb_vertex, nb_edges, nb_minimal_bases, basis_dimension, face_basis;
    int nb_face_basis_outer = 0;

    /* Read fields: V E B D face_basis nb_outer_bases */
    const int n_read = fscanf(file, "%d %d %d %d %d %d",
                        &nb_vertex, &nb_edges, &nb_minimal_bases,
                        &basis_dimension, &face_basis, &nb_face_basis_outer);
    if (n_read < 5) {
        fclose(file); perror("load_graph: header"); exit(1);
    }

    if (nb_face_basis_outer < 0) nb_face_basis_outer = 0;

    /* Read the (possibly empty) list of outer-basis indices that follows on the same line */
    int *face_basis_outer_arr = NULL;
    if (nb_face_basis_outer > 0) {
        face_basis_outer_arr = malloc(nb_face_basis_outer * sizeof(int));
        if (!face_basis_outer_arr) { perror("load_graph: face_basis_outer_arr"); exit(1); }
        for (int i = 0; i < nb_face_basis_outer; i++) {
            if (fscanf(file, "%d", &face_basis_outer_arr[i]) != 1) {
                fclose(file); perror("load_graph: face_basis_outer index"); exit(1);
            }
        }
    }
    /* Consume the rest of the header line */
    { char tmp[512]; fgets(tmp, sizeof(tmp), file); }

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
        g->minimals_basis = calloc(nb_minimal_bases, sizeof(Minimal_basis));
        if (!g->minimals_basis) { perror("load_graph: minimals_basis"); exit(1); }
        g->nb_minimal_bases = nb_minimal_bases;
        g->basis_dimension = basis_dimension;
        g->face_basis = face_basis;
        g->nb_face_basis_outer = nb_face_basis_outer;
        g->face_basis_outer = face_basis_outer_arr;

        fprintf(stdout, "Loading bases (nb=%d dim=%d)...\n",
                nb_minimal_bases, basis_dimension); fflush(stdout);

        for (int b = 0; b < nb_minimal_bases; b++) {
            g->minimals_basis[b].cycles = calloc(basis_dimension, sizeof(Path));
            if (!g->minimals_basis[b].cycles) { perror("load_graph: cycles"); exit(1); }

            /* is_faces: check if this basis is the face basis */
            g->minimals_basis[b].is_faces =
                (face_basis >= 0 && face_basis == b) ? 1 : 0;

            /* is_faces_outer: check if b appears in the outer-basis array */
            g->minimals_basis[b].is_faces_outer = 0;
            for (int k = 0; k < nb_face_basis_outer; k++) {
                if (face_basis_outer_arr && face_basis_outer_arr[k] == b) {
                    g->minimals_basis[b].is_faces_outer = 1;
                    break;
                }
            }

            /* Line buffer on the heap, not the stack. */
            const int LINE_BUF = 1 << 20; /* 1 MB in order to handle very long cycle lines */
            char *line = malloc(LINE_BUF);
            if (!line) { perror("load_graph: line buffer"); exit(1); }

            for (int i = 0; i < basis_dimension; i++) {
                g->minimals_basis[b].cycles[i].edges_ids =
                    calloc(g->nb_edges, sizeof(uint32_t));
                if (!g->minimals_basis[b].cycles[i].edges_ids) {
                    free(line); perror("load_graph: edges_ids"); exit(1);
                }

                if (!fgets(line, LINE_BUF, file)) break;

                int count = 0;
                const char *tok = strtok(line, " \t\r\n");
                while (tok) {
                    const int edge_id = atoi(tok);
                    if (edge_id >= 0 && edge_id < g->nb_edges) {
                        g->minimals_basis[b].cycles[i].edges_ids[edge_id] = 1;
                        count++;
                    }
                    tok = strtok(NULL, " \t\r\n");
                }
                g->minimals_basis[b].cycles[i].length = count;
            }
            free(line);
        }
    }

    fprintf(stdout, "Load done.\n"); fflush(stdout);
    fclose(file);
}

/**
 * @brief Add a vertex at the given coordinates.
 *
 * Appends a new Vertex and initializes its (empty) neighbor list. Doubles array capacities when
 * needed.
 *
 * @param g Graph to modify.
 * @param x X coordinate.
 * @param y Y coordinate.
 */
void create_vertex(Graph *g, const double x, const double y) {
    invalidate_bases(g);

    if (g->nb_vertex == g->capacity_vertices) {
        const int old_cap = g->capacity_vertices;
        g->capacity_vertices *= 2;
        const int new_cap = g->capacity_vertices;

        Vertex *temp_vertex = realloc(g->vertices, new_cap * sizeof(Vertex));
        if (!temp_vertex) { perror("create_vertex: vertices"); exit(1); }
        g->vertices = temp_vertex;

        Neighbor_list *temp_neighbors = realloc(g->neighbors, new_cap * sizeof(Neighbor_list));
        if (!temp_neighbors) { perror("create_vertex: neighbors"); exit(1); }
        g->neighbors = temp_neighbors;

        for (int i = old_cap; i < new_cap; i++) {
            g->neighbors[i].neighbors = NULL;
            g->neighbors[i].angles = NULL;
            g->neighbors[i].count = 0;
            g->neighbors[i].capacity = 0;
        }

        Path *temp_faces = realloc(g->faces,
            (2 + g->capacity_edges + new_cap) * sizeof(Path));
        if (!temp_faces) { perror("create_vertex: faces"); exit(1); }
        g->faces = temp_faces;
    }

    const int id = g->nb_vertex;
    g->vertices[id] = (Vertex){x, y, id, 0, 0};
    g->nb_vertex++;
}

/**
 * @brief Add an edge between two existing vertices.
 *
 * Inserts the edge into the edges array and updates both neighbor lists, maintaining sorted
 * polar-angle order. Doubles array capacities when needed.
 *
 * @param g Graph to modify.
 * @param v1_id ID of the first vertex (must be valid and not deleted).
 * @param v2_id ID of the second vertex (must be valid and not deleted).
 */
void create_edge(Graph *g, const int v1_id, const int v2_id) {
    invalidate_bases(g);

    if (g->nb_edges == g->capacity_edges) {
        g->capacity_edges *= 2;

        Edge *temp_edges = realloc(g->edges, g->capacity_edges * sizeof(Edge));
        if (!temp_edges) { perror("create_edge: edges"); exit(1); }
        g->edges = temp_edges;

        Path *temp_faces = realloc(g->faces,
            (2 + g->capacity_edges + g->capacity_vertices) * sizeof(Path));
        if (!temp_faces) { perror("create_edge: faces"); exit(1); }
        g->faces = temp_faces;
    }

    /* Update sorted neighbor lists for both endpoints */
    const int ids[2] = {v1_id, v2_id};
    for (int i = 0; i < 2; i++) {
        const double angle = atan2l(
            g->vertices[ids[(i + 1) % 2]].y - g->vertices[ids[i]].y,
            g->vertices[ids[(i + 1) % 2]].x - g->vertices[ids[i]].x);

        Neighbor_list *neigh_list = &g->neighbors[ids[i]];
        if (neigh_list->count == neigh_list->capacity) {
            neigh_list->capacity = neigh_list->capacity == 0 ? 4 : neigh_list->capacity * 2;
            int *temp_neigh = realloc(neigh_list->neighbors, neigh_list->capacity * sizeof(int));
            if (!temp_neigh) { perror("create_edge: nl->neighbors"); exit(1); }
            neigh_list->neighbors = temp_neigh;
            double *ta = realloc(neigh_list->angles, neigh_list->capacity * sizeof(double));
            if (!ta) { perror("create_edge: nl->angles"); exit(1); }
            neigh_list->angles = ta;
        }

        /* Insertion sort to maintain ascending angle order */
        int pos = neigh_list->count;
        while (pos > 0 && neigh_list->angles[pos - 1] > angle) {
            neigh_list->neighbors[pos] = neigh_list->neighbors[pos - 1];
            neigh_list->angles[pos] = neigh_list->angles[pos - 1];
            pos--;
        }
        neigh_list->neighbors[pos] = ids[(i + 1) % 2];
        neigh_list->angles[pos] = angle;
        neigh_list->count++;
    }

    g->edges[g->nb_edges] =
        (Edge){v1_id, v2_id, g->nb_edges, g->nb_edges, 0, 0, {-1, -1}};
    g->nb_edges++;
}

/**
 * @brief Internal edge deletion — marks the slot deleted and updates neighbor lists, but does NOT
 * compact the arrays.
 *
 * Must be used instead of delete_edge() when several edges are removed in a single pass
 * (like delete_vertex) so that IDs remain stable throughout the loop. The caller is responsible for
 * calling compact_graph() afterward.
 */
static void delete_edge_impl(Graph *g, const int e_id) {
    invalidate_bases(g);

    g->edges[e_id].deleted = 1;

    const int u = g->edges[e_id].u;
    const int v = g->edges[e_id].v;

    /* Remove v from u's neighbour list */
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
            g->neighbors[v].angles[i] = g->neighbors[v].angles[last];
            break;
        }
    }
}

/**
 * @brief Compact the graph arrays by removing all deleted vertex/edge slots.
 *
 * Does nothing if no vertex or edge has been deleted since the last compaction.
 * invalidate_bases() must have been called before compact_graph() so that bases
 * and faces (which embed old IDs) are already freed.
 *
 * @param g Graph to compact.
 */
void compact_graph(Graph *g) {
    if (!g) return;

    /* Build remapping tables */
    int *vmap = malloc(g->nb_vertex * sizeof(int));
    int *emap = malloc(g->nb_edges  * sizeof(int));
    if (!vmap || !emap) { free(vmap); free(emap); perror("compact_graph: maps"); exit(1); }

    int new_nv = 0;
    for (int i = 0; i < g->nb_vertex; i++)
        vmap[i] = g->vertices[i].deleted ? -1 : new_nv++;

    int new_ne = 0;
    for (int i = 0; i < g->nb_edges; i++)
        emap[i] = g->edges[i].deleted ? -1 : new_ne++;

    /* Nothing to compact */
    if (new_nv == g->nb_vertex && new_ne == g->nb_edges) {
        free(vmap); free(emap);
        return;
    }

    /* Compact vertices, update id and label */
    for (int i = 0; i < g->nb_vertex; i++) {
        if (vmap[i] < 0) continue;
        g->vertices[vmap[i]] = g->vertices[i];
        g->vertices[vmap[i]].id = vmap[i];
        g->vertices[vmap[i]].label = vmap[i];
    }
    g->nb_vertex = new_nv;

    /* Compact edges, remap vertices and id */
    for (int i = 0; i < g->nb_edges; i++) {
        if (emap[i] < 0) continue;
        g->edges[emap[i]] = g->edges[i];
        g->edges[emap[i]].u = vmap[g->edges[i].u];
        g->edges[emap[i]].v = vmap[g->edges[i].v];
        g->edges[emap[i]].id = emap[i];
    }
    g->nb_edges = new_ne;

    /* Rebuild neighbor lists */
    for (int i = 0; i < g->capacity_vertices; i++) {
        free(g->neighbors[i].neighbors);
        free(g->neighbors[i].angles);
        g->neighbors[i].neighbors = NULL;
        g->neighbors[i].angles = NULL;
        g->neighbors[i].count = 0;
        g->neighbors[i].capacity = 0;
    }

    /* Re-insert each edge using insertion sort by angle (same as create_edge) */
    for (int e = 0; e < new_ne; e++) {
        const int ep[2] = {g->edges[e].u, g->edges[e].v};
        for (int s = 0; s < 2; s++) {
            const int src = ep[s], dst = ep[s ^ 1];
            const double angle = atan2l(
                g->vertices[dst].y - g->vertices[src].y,
                g->vertices[dst].x - g->vertices[src].x);

            Neighbor_list *nl = &g->neighbors[src];
            if (nl->count == nl->capacity) {
                nl->capacity = nl->capacity == 0 ? 4 : nl->capacity * 2;
                int    *tn = realloc(nl->neighbors, nl->capacity * sizeof(int));
                double *ta = realloc(nl->angles,    nl->capacity * sizeof(double));
                if (!tn || !ta) { perror("compact_graph: neighbor realloc"); exit(1); }
                nl->neighbors = tn;
                nl->angles    = ta;
            }

            /* Insertion sort to maintain ascending angle order */
            int pos = nl->count;
            while (pos > 0 && nl->angles[pos - 1] > angle) {
                nl->neighbors[pos] = nl->neighbors[pos - 1];
                nl->angles[pos]    = nl->angles[pos - 1];
                pos--;
            }
            nl->neighbors[pos] = dst;
            nl->angles[pos]    = angle;
            nl->count++;
        }
    }

    /* Free BFS matrices (reallocated by prepare_graph_matrices) */
    free(g->edge_indices); g->edge_indices = NULL;
    free(g->predecessors); g->predecessors = NULL;
    free(g->distances); g->distances = NULL;

    free(vmap);
    free(emap);
}

/**
 * @brief Delete an edge and compact the graph arrays.
 *
 * @param g Graph to modify.
 * @param e_id Index of the edge to delete.
 */
void delete_edge(Graph *g, const int e_id) {
    delete_edge_impl(g, e_id);
    invalidate_bases(g);
    compact_graph(g);
}

/**
 * @brief Delete a vertex and all its incident edges, then compact the graph.
 *
 * @param g Graph to modify.
 * @param v_id ID of the vertex to delete.
 */
void delete_vertex(Graph *g, const int v_id) {
    for (int e = 0; e < g->nb_edges; e++) {
        if (g->edges[e].deleted) continue;
        if (g->edges[e].u == v_id || g->edges[e].v == v_id)
            delete_edge_impl(g, e);
    }
    g->vertices[v_id].deleted = 1;
    invalidate_bases(g);
    compact_graph(g);
}

/**
 * @brief Subdivide an edge by inserting number_vertex_to_add intermediate vertices.
 *
 * Deletes edge (u, v) and replaces it with the path
 * u − w1 − … − wk − v, where each wi lies at equal intervals along (u, v).
 *
 * @param g Graph to modify.
 * @param edge_id Index of the edge to split.
 * @param number_vertex_to_add Number of intermediate vertices to insert (k ≥ 1).
 */
void split_edge(Graph *g, const int edge_id, const int number_vertex_to_add) {
    const int u = g->edges[edge_id].u;
    const int v = g->edges[edge_id].v;
    const double ux = g->vertices[u].x, uy = g->vertices[u].y;
    const double vx = g->vertices[v].x, vy = g->vertices[v].y;

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
 * The number_vertex_to_add vertices are distributed across all edges as
 * evenly as possible; the first (number_vertex_to_add % number_edge_to_split)
 * edges each receive one extra vertex.
 *
 * @param g Graph to modify.
 * @param edge_ids Array of edge indices to split.
 * @param number_edge_to_split Number of edges to split.
 * @param number_vertex_to_add Total number of intermediate vertices to insert.
 */
void split_edges(Graph *g, const int *edge_ids,
                 const int number_edge_to_split, const int number_vertex_to_add) {
    if (number_edge_to_split <= 0) return;
    const int base = number_vertex_to_add / number_edge_to_split;
    const int extra = number_vertex_to_add % number_edge_to_split;
    for (int i = 0; i < number_edge_to_split; i++)
        split_edge(g, edge_ids[i], base + (i < extra ? 1 : 0));
}

/**
 * @brief Sort a vertex's neighbor list by polar angle using atan2 in ascending order.
 *
 * @param g Graph.
 * @param v_id Vertex whose neighbor list is to be sorted.
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
 * neighbors, re-sorts neighbor lists, then recomputes the face decomposition.
 *
 * @param g Graph to modify.
 * @param v ID of the vertex to move.
 * @param new_x New X coordinate.
 * @param new_y New Y coordinate.
 */
void move_vertex(Graph *g, const int v, const double new_x, const double new_y) {
    g->vertices[v].x = new_x;
    g->vertices[v].y = new_y;

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
 * @brief Allocate and initialize the n×n matrices for Horton's algorithm.
 *
 * Allocates (or reallocates) three n×n integer matrices:
 *  - edge_indices[u*n+v]: index of edge (u,v) in g->edges, or -1,
 *  - predecessors[u*n+v]: predecessor of v on the shortest path from u,
 *  - distances[u*n+v]: BFS distance from u to v.
 *
 * All entries initialized to -1; edge_indices filled for non-deleted edges.
 *
 * @param graph Graph to prepare.
 */
void prepare_graph_matrices(Graph *graph) {
    if (!graph->edges || graph->nb_vertex == 0) {
        printf("prepare_graph_matrices: empty graph, skipping.\n");
        return;
    }

    const int n = graph->nb_vertex;
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
        const int u = graph->edges[i].u;
        const int v = graph->edges[i].v;
        graph->edge_indices[u * n + v] = i;
        graph->edge_indices[v * n + u] = i;
    }
}

/**
 * @brief Compute all-pairs BFS shortest paths and predecessor tables.
 *
 * Runs one BFS per non-deleted source vertex. Then reconstructs the predecessor[u*n+v] entry for
 * every reachable pair (u, v) by walking backward along the distance array. When several
 * predecessors yield the same distance, the one with the smallest vertex label is chosen.
 *
 * @param g Graph whose matrices (allocated by prepare_graph_matrices) are filled.
 */
void compute_all_shortest_paths(const Graph *g) {
    int *queue = malloc(g->nb_vertex * sizeof(int));
    if (!queue) { perror("compute_all_shortest_paths: queue"); exit(1); }

    /* BFS from each vertex */
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

    /* Predecessors finding in order to compute shortest paths */
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
 * Follows the predecessor chain stored in g->predecessors. The returned Path
 * owns its edges_ids, vertices_ids, and edges_labels arrays; the caller must
 * free them.
 *
 * @param g Graph with precomputed shortest paths.
 * @param u Source vertex.
 * @param v Destination vertex.
 * @return A Path from u to v.
 */
Path create_path(const Graph *g, const int u, const int v) {
    if (!g->predecessors || !g->edges) {
        perror("create_path: shortest paths not computed"); exit(1);
    }

    Path p;
    p.length = 0;
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

    /* Path computation */
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

/*
 * @brief Find the outer face among all faces and mark its edges as outer edges. Update the
 * g->outer_face field.
 *
 * @param g Graph with computed faces.
 */
static void find_outer_face(Graph *g) {

    int leftmost = -1;
    for (int i = 0; i < g->nb_vertex; i++) {
        if (g->vertices[i].deleted || g->neighbors[i].count == 0) continue;
        if (leftmost == -1 || g->vertices[i].x < g->vertices[leftmost].x)
            leftmost = i;
    }

    if (leftmost == -1) { g->outer_face = 0; return; }

    const int neighbor = g->neighbors[leftmost].neighbors[0];
    int u = leftmost, v = neighbor, nb_outer_edge = 0;
    do {
        const int eid = g->edge_indices[u * g->nb_vertex + v];
        if (eid >= 0) g->edges[eid].is_outer = 1;
        nb_outer_edge++;
        const Neighbor_list *nl = &g->neighbors[v];
        int idx = -1;
        for (int j = 0; j < nl->count; j++) {
            if (nl->neighbors[j] == u) { idx = j; break; }
        }

        if (idx < 0) break;
        const int w = nl->neighbors[(idx + 1) % nl->count];
        u = v; v = w;

    } while (u != leftmost || v != neighbor);

    for (int f=0; f < g->nb_faces; f++) {
        const Path face = g->faces[f];
        if (face.length == nb_outer_edge) {
            int match = 1;
            for (int e = 0; e < g->nb_edges; e++) {
                if (face.edges_ids[e] != 0 && !g->edges[e].is_outer) {
                    match = 0;
                    break;
                }
            }
            if (match) {
                g->outer_face = f;
                return;
            }
        }
    }
}

/**
 * @brief Compute the faces of the planar graph and identify the outer face.
 *
 * @param g Graph to analyze. Must have its neighbor lists sorted in CCW order.
 */
void find_faces(Graph *g) {

    if (g->faces) {
        for (int i = 0; i < g->nb_faces; i++) {
            free(g->faces[i].edges_ids);
            free(g->faces[i].edges_labels);
            free(g->faces[i].vertices_ids);
            g->faces[i].edges_ids    = NULL;
            g->faces[i].edges_labels = NULL;
            g->faces[i].vertices_ids = NULL;
        }
    }
    g->nb_faces  = 0;
    g->outer_face = -1;

    for (int e = 0; e < g->nb_edges; e++) {
        g->edges[e].face_id[0] = -1;
        g->edges[e].face_id[1] = -1;
        g->edges[e].is_outer = 0;
    }

    const int n = g->nb_vertex;
    const int m = g->nb_edges;
    if (!g->edge_indices) prepare_graph_matrices(g);

    /* 'visited' tracks which directed half-edges (u->p) are processed */
    int *visited = calloc(n * n, sizeof(int));
    if (!visited) { perror("find_faces: visited"); exit(1); }

    /* Main Loop: iterate through every vertex and its neighbors */
    for (int v = 0; v < n; v++) {
        if (g->vertices[v].deleted) continue;

        for (int i = 0; i < g->neighbors[v].count; i++) {
            const int neigh = g->neighbors[v].neighbors[i];

            /* If the edge (v, neigh) hasn't been part of a face yet */
            if (visited[v * n + neigh]) continue;

            Path *face = &g->faces[g->nb_faces];
            face->edges_ids = calloc(m, sizeof(uint32_t));
            face->vertices_ids = calloc(n, sizeof(uint32_t));
            face->edges_labels = calloc(m, sizeof(uint32_t));
            face->length = 0;

            int u = v, p = neigh;

            /* walk around the face by following the "right-hand rule" */
            while (!visited[u * n + p]) {
                visited[u * n + p] = 1;
                face->length++;

                /* Assign the face ID to the current edge */
                const int eid = g->edge_indices[u * n + p];
                if (eid >= 0) {
                    face->edges_ids[eid] = 1;
                    /* Guard: label may be larger than current nb_edges
                     * if find_faces() is called outside a full Horton run. */
                    const int lbl = g->edges[eid].label;
                    if (lbl >= 0 && lbl < m)
                        face->edges_labels[lbl] = 1;
                    if (g->edges[eid].face_id[0] == -1)
                        g->edges[eid].face_id[0] = g->nb_faces;
                    else
                        g->edges[eid].face_id[1] = g->nb_faces;
                }
                face->vertices_ids[u] = 1;

                /*Select the next vertex by turning "right" by looking for the current vertex 'u'
                in the neighbor list of destination 'p' */
                const Neighbor_list *nl = &g->neighbors[p];
                int idx = -1;
                for (int j = 0; j < nl->count; j++) {
                    if (nl->neighbors[j] == u) { idx = j; break; }
                }

                /* Pick the next neighbor in the CCW sorted list (idx + 1). This corresponds to a
                 * "right turn" in the embedding */
                const int w = nl->neighbors[(idx + 1 + nl->count) % nl->count];
                u = p; p = w; // Move to the next edge (p, w)
            }
            g->nb_faces++;
        }
    }

    find_outer_face(g);
    free(visited);
}

/**
 * @brief Append a cycle to the graph's Horton cycle list.
 *
 * @param g Graph to which the cycle is appended.
 * @param c Cycle to add.
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
                    .edges_ids = ec,
                    .vertices_ids = vc,
                    .edges_labels = el,
                    .length = pv_x.length + pv_y.length + 1
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
 *
 * @param g Graph whose horton_cycles are to be sorted.
 */
void order_horton_cycles(const Graph *g) {
    nb_edges_for_qsort = g->nb_edges;
    if (g->nb_horton_cycles > 1)
        qsort(g->horton_cycles, g->nb_horton_cycles, sizeof(Path),
            compare_cycles);
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
 * @param c1 First cycle.
 * @param c2 Second cycle.
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
 * @param a First basis.
 * @param b Second basis.
 * @param nb_edges Size of the edges_ids vectors.
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
        s1[i].length = mb1->cycles[i].length;
        s1[i].edges_ids = malloc(nb_edges * sizeof(uint32_t));
        s2[i].length = mb2->cycles[i].length;
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
 * @param g Graph containing the known bases.
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
 * Compares the candidate basis mb (sorted by edge set) to the interior
 * faces of the graph (all faces except the outer face, also sorted).
 *
 * @param g  Graph with a precomputed face decomposition.
 * @param mb Candidate minimal basis to test.
 * @return 1 if mb matches the face basis, 0 otherwise.
 */
int mb_is_faces(Graph *g, const Minimal_basis *mb) {

    if (!g->faces || g->nb_faces == 0) find_faces(g);

    /* If find_outer_face() failed (because of disconnected graph), outer_face is -1.
     * In that case we cannot distinguish "interior faces only" from "outer face
     * included", so return 0 to avoid incorrectly tagging bases. */
    if (g->outer_face < 0) return 0;

    int D = g->basis_dimension;

    if (g->nb_faces != D + 1) return 0;

    /* Try all combination of faces, including the outer one */
    for (int skip = 0; skip < g->nb_faces; skip++) {
        int match_count = 0;
        for (int i = 0; i < D; i++) {
            int found = 0;
            for (int j = 0; j < g->nb_faces; j++) {
                if (j == skip) continue;
                if (are_same_cycles(&mb->cycles[i], &g->faces[j], g->nb_edges)) {
                    found = 1;
                    break;
                }
            }

            if (found) match_count++;
            else break;
        }

        if (match_count == D) {
            if (skip == g->outer_face) return 1; // matches interior faces
            return 2; // matches faces but includes the outer face
        }
    }

    return 0;
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
 * all old bases are cleared and the new basis replaces them.
 *
 * The is_faces flag is set when the computed basis matches the interior face set
 * (Euler: F = E − V + 2, so nb_faces = basis_dimension + 1).
 *
 * @param g Graph. g->horton_cycles must be non-empty and sorted.
 */
void greedy_cycle_cover(Graph *g) {
    if (!g->horton_cycles || g->nb_horton_cycles == 0) {
        g->basis_dimension = 0;
        printf("greedy_cycle_cover: no Horton cycles, skipping.\n");
        return;
    }

    uint32_t **basis = calloc(g->nb_edges, sizeof(uint32_t *));
    Path *mcb   = calloc(g->nb_horton_cycles, sizeof(Path));
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
             * triggering invalidate_bases(). Discard all old bases.
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
            g->minimals_basis = NULL;
            g->nb_minimal_bases = 0;
            g->face_basis = -1;
            g->nb_face_basis_outer = 0;
            free(g->face_basis_outer);
            g->face_basis_outer = NULL;
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
     * Check if this basis matches the face decomposition.
     *
     * We pass &mb (the newly computed MCB) to mb_is_faces, which compares each
     * cycle of the MCB against the graph's faces.
     *
     * mb_is_faces() returns:
     *   1  — MCB equals the D interior faces
     *   2  — MCB equals outer face + D-1 interior faces
     *   0  — no match
     */
    mb.is_faces = 0;
    mb.is_faces_outer = 0;
    if (g->faces && g->nb_faces == mcb_count + 1) {
        const int result = mb_is_faces(g, &mb);
        if (result == 1) {
            mb.is_faces   = 1;
            g->face_basis = g->nb_minimal_bases;
        } else if (result == 2) {
            mb.is_faces_outer = 1;
            int *tmp = realloc(g->face_basis_outer,
                               (g->nb_face_basis_outer + 1) * sizeof(int));
            if (!tmp) { perror("greedy_cycle_cover: realloc face_basis_outer"); exit(1); }
            g->face_basis_outer = tmp;
            g->face_basis_outer[g->nb_face_basis_outer] = g->nb_minimal_bases;
            g->nb_face_basis_outer++;
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
 * Called at the start of each Horton run so old cycles from the previous
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
 *  1. Clear previous Horton cycle list,
 *  2. Allocate edge_indices / BFS matrices,
 *  3. BFS from every non-deleted vertex to find shortest paths,
 *  4. Face decomposition for face-basis detection,
 *  5. Enumerate candidate cycles,
 *  6. Sort cycles by length and label order,
 *  7. Extract one minimal basis.
 *
 * @param g Graph on which to run the algorithm.
 */
void horton(Graph *g) {
    reset_graph_results(g);
    prepare_graph_matrices(g);
    compute_all_shortest_paths(g);
    if (g->faces == NULL || g->nb_faces == 0) find_faces(g);
    find_horton_cycles(g);
    order_orton_cycles_permutations(g);
    greedy_cycle_cover(g);
}

/**
 * @brief Try to add the interior face set as a minimal cycle basis.
 *
 * For a connected planar graph the interior faces always span the cycle space, so they always form
 * a basis — no independence check needed. However, we need to check if this  basis is *minimal*,
 * (if the sum of the face lengths equals the total weight of the minimal bases already found by
 * Horton). If so, the face basis is appended to the list (or tagged if Horton
 * already found it without recognizing it as the face basis).
 *
 * @param g Graph with a precomputed face decomposition.
 */
static void try_add_face_basis(Graph *g) {
    if (!g->faces || g->nb_faces == 0) return;
    if (g->outer_face < 0) return;
    if (g->basis_dimension <= 0) return;
    if (g->nb_faces != g->basis_dimension + 1) return;  /* Euler: F = E-V+2 */

    /* Sum of interior face lengths. */
    int face_weight = 0;
    for (int f = 0; f < g->nb_faces; f++) {
        if (f == g->outer_face) continue;
        face_weight += g->faces[f].length;
    }

    /* Sum of cycle lengths of any known minimal basis (all MCBs have the same
     * total weight — use the first one). */
    if (g->nb_minimal_bases == 0) return;
    int basis_weight = 0;
    for (int c = 0; c < g->basis_dimension; c++)
        basis_weight += g->minimals_basis[0].cycles[c].length;

    /* The face set is minimal only if weights match. */
    if (face_weight != basis_weight) return;

    /* Build the candidate and check for duplicates. */
    Path *mcb = calloc(g->basis_dimension, sizeof(Path));
    if (!mcb) { perror("try_add_face_basis: mcb"); exit(1); }

    int ci = 0;
    for (int f = 0; f < g->nb_faces; f++) {
        if (f == g->outer_face) continue;
        mcb[ci].length      = g->faces[f].length;
        mcb[ci].edges_ids   = malloc(g->nb_edges * sizeof(uint32_t));
        mcb[ci].edges_labels = calloc(g->nb_edges, sizeof(uint32_t));
        if (!mcb[ci].edges_ids || !mcb[ci].edges_labels) {
            perror("try_add_face_basis: mcb entry"); exit(1);
        }
        memcpy(mcb[ci].edges_ids, g->faces[f].edges_ids,
               g->nb_edges * sizeof(uint32_t));
        ci++;
    }

    Minimal_basis mb = { .cycles = mcb, .is_faces = 1, .is_faces_outer = 0 };

    if (already_exists_basis(g, mb)) {
        /* Horton found this basis but didn't tag it — fix the tag. */
        for (int b = 0; b < g->nb_minimal_bases; b++) {
            if (are_same_minimal_bases(&g->minimals_basis[b], &mb,
                                       g->nb_edges, g->basis_dimension)) {
                g->minimals_basis[b].is_faces = 1;
                g->face_basis = b;
                break;
            }
        }
        for (int c = 0; c < g->basis_dimension; c++) {
            free(mcb[c].edges_ids);
            free(mcb[c].edges_labels);
        }
        free(mcb);
        return;
    }

    /* New basis — append it. */
    g->face_basis = g->nb_minimal_bases;
    g->nb_minimal_bases++;
    Minimal_basis *tmp = realloc(g->minimals_basis,
                                 g->nb_minimal_bases * sizeof(Minimal_basis));
    if (!tmp) { perror("try_add_face_basis: realloc"); exit(1); }
    g->minimals_basis = tmp;
    g->minimals_basis[g->nb_minimal_bases - 1] = mb;
}

/**
 * @brief Return the maximum number of cycles any single edge appears in, across all cycles of a
 * given basis.
 */
static int max_edge_cycle_count(const Minimal_basis *mb,
                                const int nb_cycles, const int nb_edges) {
    int *count = calloc(nb_edges, sizeof(int));
    if (!count) { perror("max_edge_cycle_count"); exit(1); }
    for (int c = 0; c < nb_cycles; c++)
        for (int e = 0; e < nb_edges; e++)
            if (mb->cycles[c].edges_ids[e]) count[e]++;
    int mx = 0;
    for (int e = 0; e < nb_edges; e++) if (count[e] > mx) mx = count[e];
    free(count);
    return mx;
}

/**
 * @brief Update g->no_face_basis_possible.
 *
 * In any planar embedding, each interior edge borders exactly 2 faces, so each
 * edge can appear in at most 2 cycles of any face-basis. If every stored basis
 * has at least one edge that appears in 3 or more cycles, then no drawing of
 * the graph can yield a face decomposition equal to any known basis, and the
 * flag is set to 1. It is reset to 0 as soon as a face basis is found.
 */
static void update_no_face_basis_possible(Graph *g) {
    if (g->face_basis >= 0) { g->no_face_basis_possible = 0; return; }
    if (g->nb_minimal_bases == 0) return;

    for (int b = 0; b < g->nb_minimal_bases; b++) {
        const int mx = max_edge_cycle_count(&g->minimals_basis[b],
                                            g->basis_dimension, g->nb_edges);
        if (mx < 3) { g->no_face_basis_possible = 0; return; }
    }
    g->no_face_basis_possible = 1;
}

/**
 * @brief Run Horton's algorithm multiple times with different edge-label permutations.
 *
 * @param g Graph on which to run.
 * @param inv Caller-allocated inversion table buffer (size: g->nb_edges).
 * @param max_permutations Number of permuted runs to execute.
 */
void multiple_horton(Graph *g, int *inv, const int max_permutations) {
    if (!inv || !g->edges) {
        printf("multiple_horton: missing inversion table or edges, skipping.\n");
        return;
    }

    find_faces(g);

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

    /* If Horton never found the face basis, try adding it explicitly. */
    if (g->face_basis < 0)
        try_add_face_basis(g);

    /* Set the diagnostic flag (no face basis reachable regardless of drawing). */
    update_no_face_basis_possible(g);
}