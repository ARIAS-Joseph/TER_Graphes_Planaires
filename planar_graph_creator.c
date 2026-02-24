#include "planar_graph_creator.h"
#include <math.h>
#include "graph.h"

/**
 * @brief Check if two edges intersect in the plane.
 *
 * This function checks if the two edges a and b intersect in the plane. It first checks if the
 * edges share a vertex, in which case they do not intersect. Then it uses the cross product to check
 * if the edges intersect in the plane.
 *
 * @param a Pointer to the first edge.
 * @param b Pointer to the second edge.
 * @param g Pointer to the graph structure, used to access the vertices of the edges.
 * @return 1 if the edges intersect, 0 otherwise.
 */
int intersect(Edge *a, Edge *b, const Graph *g) {

    // Check if the edges share a vertex, in which case they do not intersect
    if (a->u == b->u || a->u == b->v || a->v == b->u || a->v == b->v) {
        return 0;
    }

    // Get the coordinates of the vertices of the edges
    const double e1_x1 = g->vertices[a->u].x;
    const double e1_y1 = g->vertices[a->u].y;
    const double e1_x2 = g->vertices[a->v].x;
    const double e1_y2 = g->vertices[a->v].y;
    const double e2_x1 = g->vertices[b->u].x;
    const double e2_y1 = g->vertices[b->u].y;
    const double e2_x2 = g->vertices[b->v].x;
    const double e2_y2 = g->vertices[b->v].y;

    // Compute the cross products to check if the edges intersect
    const double d1 = (e1_x2 - e1_x1) * (e2_y1 - e1_y1) - (e1_y2 - e1_y1) * (e2_x1 - e1_x1);
    const double d2 = (e1_x2 - e1_x1) * (e2_y2 - e1_y1) - (e1_y2 - e1_y1) * (e2_x2 - e1_x1);
    const double d3 = (e2_x2 - e2_x1) * (e1_y1 - e2_y1) - (e2_y2 - e2_y1) * (e1_x1 - e2_x1);
    const double d4 = (e2_x2 - e2_x1) * (e1_y2 - e2_y1) - (e2_y2 - e2_y1) * (e1_x2 - e2_x1);

    // The edges intersect if d1 and d2 have different signs and d3 and d4 have different signs
    return (d1 * d2 < 0) && (d3 * d4 < 0);

}

/**
 * @brief Try to add an edge between two random vertices.
 *
 * This function attempts to add a random edge between two distinct vertices in the graph. It first
 * checks if the edge already exists. Then it checks if the candidate
 * edge intersects with any existing edge in the graph using the intersect function. If it does not
 * intersect with any existing edge, it adds the edge to the graph and returns 1. If it intersects
 * with at least one existing edge, it does not add the edge and returns 0.
 *
 * @param g Pointer to the graph structure where the edge will be added if it does not create
 * intersections.
 * @return 1 if the edge was added successfully, 0 otherwise.
 */
int try_add_edge(Graph* g, const int v1, const int v2) {

    int already_exists = 0;
    for (int e = 0; e < g->nb_edges; e++) {
        if ((g->edges[e].u == v1 && g->edges[e].v == v2) ||
            (g->edges[e].u == v2 && g->edges[e].v == v1)) {
            already_exists = 1;
            break;
            }
    }
    if (already_exists) return 0;

    // Check if the candidate edge intersects with any existing edge
    Edge candidate_edge = (Edge){v1, v2, 1, -1, -1};
    int intersects = 0;
    for (int e = 0; e < g->nb_edges; e++) {
        if (intersect(&candidate_edge, &g->edges[e], g)) {
            intersects = 1;
            break;
        }
    }

    // Add the edge if it doesn't intersect
    if (!intersects) {
        create_edge(g, v1, v2, 1);
        return 1;
    }
    return 0;
}

void spread_label(Graph* g, const int edge_id) {
    const int v1_id = g->edges[edge_id].u;
    const int v2_id = g->edges[edge_id].v;

    const int label_v1 = g->vertices[v1_id].label;
    const int label_v2 = g->vertices[v2_id].label;

    const int min_label = label_v1 < label_v2 ? label_v1 : label_v2;
    g->vertices[v1_id].label = min_label;
    g->vertices[v2_id].label = min_label;

    for (int e = 0; e < g->nb_edges; e++) {
        if (g->edges[e].u == v1_id || g->edges[e].v == v1_id || g->edges[e].u == v2_id || g->edges[e].v == v2_id) {
            const int other_vertex_id = (g->edges[e].u == v1_id || g->edges[e].u == v2_id) ? g->edges[e].v : g->edges[e].u;
            if (g->vertices[other_vertex_id].label != min_label) {
                spread_label(g, e);
            }
        }
    }
}

void create_tree(Graph* g) {

    for (int i = 0; i < g->nb_vertex; i++) {
        g->vertices[i].label = i;
    }

    int unique_label = 0;
    int cnt = 0;
    while (!unique_label && cnt < 100000) {

        for (int i = 0; i < g->nb_vertex; i++) {
            if (g->vertices[i].label != g->vertices[(i+1)%g->nb_vertex].label) {
                int added = try_add_edge(g, i, (i+1)%g->nb_vertex);
                if (added) {
                    spread_label(g, g->nb_edges-1);
                } else {
                    int j = 0;
                    while (j < g->nb_vertex && !added) {
                        if (g->vertices[i].label != g->vertices[j].label) {
                            added = try_add_edge(g, i, j);
                            if (added) {
                                spread_label(g, g->nb_edges-1);
                                break;
                            }
                        }
                        j++;
                    }
                }
            }
        }

        unique_label = 1;
        for (int i=0; i < g->nb_vertex; i++) {
            if (g->vertices[i].label != 0) {
                unique_label = 0;
                break;
            }
        }
        cnt++;
    }
    if (!unique_label) {
        fprintf(stderr, "Failed to create a tree after %d attempts. The graph may be not connex.");
    }
}

/**
 * @brief Create an outerplanar graph with a specified number of vertices and a target number of
 * edges.
 *
 * This function creates an outerplanar graph with a specified number of vertices and a target
 * number of edges. It starts by creating a circular layout of vertices and connecting them in a
 * cycle to ensure the graph is outerplanar. Then, it attempts to add additional edges between
 * random pairs of vertices while ensuring that the new edge does not intersect with any existing
 * edge in the graph. The process continues until the target number of edges is reached or a maximum
 * number of attempts is exceeded to avoid infinite loops.
 *
 * @param nb_vertex The number of vertices in the outerplanar graph.
 * @param nb_edges_target The target number of edges in the outerplanar graph (excluding
 * the edges of the outer cycle).
 * @return A pointer to the created outerplanar graph structure.
 */
Graph* create_outer_planar_graph(const int nb_vertex, int nb_edges_target) {

    Graph *g = create_graph();

    // Create vertices in a circular layout
    for (int i = 0; i < nb_vertex; i++) {
        const double angle = 2.0 * M_PI * i / nb_vertex;
        const double x = cos(angle);
        const double y = sin(angle);
        create_vertex(g, x, y);
    }

    // Connect vertices in a cycle
    for (int i = 0; i < nb_vertex; i++) {
        create_edge(g, i, (i + 1) % nb_vertex, 1);
    }

    // The maximum number of edges in a planar graph with n vertices is 3n - 6, and we already have
    // n edges in the cycle, so we can only add at most 2n - 6 more edges
    nb_edges_target = nb_edges_target < 2 * nb_vertex - 6 ? nb_edges_target : 2 * nb_vertex - 6;

    // Try to add more edges until we reach the target or exhaust attempts
    const int MAX_ATTEMPTS = nb_edges_target * 10;
    int attempts = 0;

    while (g->nb_edges < nb_edges_target && attempts < MAX_ATTEMPTS) {
        attempts++;
        int v1 = rand() % g->nb_vertex;
        int v2 = rand() % g->nb_vertex;

        while (v1 == v2) {
            v1 = rand() % g->nb_vertex;
            v2 = rand() % g->nb_vertex;
        }
        const int added = try_add_edge(g, v1, v2);
        if (added) {
            attempts = 0; // Reset attempts after a successful addition
        }
    }

    return g;
}

Graph* create_planar_graph(const int nb_vertex, int nb_edges_target) {
    Graph *g = create_graph();

    // Create vertices in a nb_vertex*100 * nb_vertex*100 grid
    double created[nb_vertex][2];
    for (int i = 0; i < nb_vertex; i++) {
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
        } else {
            i--;
        }
    }

    create_tree(g);

    // The maximum number of edges in a planar graph with n vertices is 3n - 6, and we already have
    // n-1 edges in the cycle, so we can only add at most 2n - 5 more edges
    nb_edges_target = nb_edges_target < 2 * nb_vertex - 5 ? nb_edges_target : 2 * nb_vertex - 5;

    // Try to add more edges until we reach the target or exhaust attempts
    const int MAX_ATTEMPTS = nb_edges_target * 10;
    int attempts = 0;

    while (g->nb_edges < nb_edges_target && attempts < MAX_ATTEMPTS) {
        attempts++;
        int v1 = rand() % g->nb_vertex;
        int v2 = rand() % g->nb_vertex;

        while (v1 == v2) {
            v1 = rand() % g->nb_vertex;
            v2 = rand() % g->nb_vertex;
        }
        const int added = try_add_edge(g, v1, v2);
        if (added) {
            attempts = 0; // Reset attempts after a successful addition
        }
    }

    return g;
}

Graph* test_circle(const int nb_vertex) {

    Graph *g = create_graph();

    // Create vertices in a circular layout
    for (int i = 0; i < nb_vertex; i++) {
        const double angle = 2.0 * M_PI * i / nb_vertex;
        const double x = cos(angle);
        const double y = sin(angle);
        create_vertex(g, x, y);
    }

    // Connect vertices in a cycle
    for (int i = 0; i < nb_vertex; i++) {
        create_edge(g, i, (i + 1) % nb_vertex, 1);
    }
    return g;
}
