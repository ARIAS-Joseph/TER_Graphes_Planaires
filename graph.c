#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "permutations.h"
#include "graph.h"
#include <limits.h>

/**
 * @brief Create an edge between two vertices in the graph.
 *
 * This function adds an edge to the graph connecting the vertices with ids v1_id and v2_id, and
 * assigns it the specified weight. It also updates the neighbor lists of the involved vertices.
 *
 * @param g Pointer to the graph structure where the edge will be added.
 * @param v1_id ID of the first vertex.
 * @param v2_id ID of the second vertex.
 * @param weight Weight of the edge.
 *
 */
void create_edge(Graph *g, const int v1_id, const int v2_id, const int weight) {
    // Double the capacity of the edges array if the current capacity is reached and reallocate
    // memory
    if (g->nb_edges == g->capacity_edges) {
        g->capacity_edges *= 2;
        Edge *tmp = realloc(g->edges, g->capacity_edges * sizeof(Edge));
        if (tmp == NULL) {
            perror("Failed to reallocate memory for Graph.edges in create_edge");
            exit(1);
        }
        g->edges = tmp;
    }

    const int ids[2] = {v1_id, v2_id};

    // For each of the two vertices involved in the edge, we update their neighbor list to include
    // the other vertex as a neighbor
    for (int i = 0; i < 2; i++) {
        Neighbor_list *nl = &g->neighbors[ids[i]]; // Get the neighbor list of the vertex with id
        // ids[i]

        // Double the capacity of the neighbors array if the current capacity is reached and reallocate
        // memory
        if (nl->count == nl->capacity) {
            nl->capacity = nl->capacity == 0 ? 4 : nl->capacity * 2;
            int *tmp = realloc(nl->neighbors, nl->capacity * sizeof(int));
            if (tmp == NULL) {
                perror("Failed to reallocate memory for Graph.neighbors in create_edge");
                exit(1);
            }
            nl->neighbors = tmp;
        }
        nl->neighbors[nl->count++] = ids[(i + 1) % 2]; // Add the other vertex as a neighbor
    }

    // Add the new edge to the edges array and increment the number of edges
    g->edges[g->nb_edges] = (Edge){v1_id, v2_id, weight, g->nb_edges, g->nb_edges};
    g->nb_edges++;
}

/**
 * @brief Create a vertex in the graph.
 *
 * This function adds a vertex to the graph with the specified coordinates (x, y). It also
 * initializes the corresponding neighbor list for the new vertex.
 *
 * @param g Pointer to the graph structure where the vertex will be added.
 * @param x X coordinate of the vertex (optional, used for visualization).
 * @param y Y coordinate of the vertex (optional, used for visualization).
 */
void create_vertex(Graph *g, const double x, const double y) {
    // Double the capacity of the vertices array if the current capacity is reached and reallocate memory
    if (g->nb_vertex == g->capacity_vertices) {
        const int old_cap = g->capacity_vertices;
        g->capacity_vertices *= 2;
        const int new_cap = g->capacity_vertices;

        Vertex *tmp = realloc(g->vertices, new_cap * sizeof(Vertex));
        if (tmp == NULL) {
            perror("Failed to reallocate memory for Graph.vertices in create_vertex");
            exit(1);
        }
        g->vertices = tmp;

        Neighbor_list *tmp2 = realloc(g->neighbors, new_cap * sizeof(Neighbor_list));
        if (tmp2 == NULL) {
            perror("Failed to reallocate memory for Graph.neighbors in create_vertex");
            exit(1);
        }
        g->neighbors = tmp2;

        // Initialize the new neighbor lists for the future new vertices
        for (int i = old_cap; i < new_cap; i++) {
            g->neighbors[i].neighbors = NULL;
            g->neighbors[i].count = 0;
            g->neighbors[i].capacity = 0;
        }
    }

    // Get the id of the new vertex, add it to the vertices array, and increment the number of
    // vertices
    const int id = g->nb_vertex;
    g->vertices[id] = (Vertex){x, y, id};
    g->nb_vertex++;
}

/**
 * @brief Create an empty graph.
 *
 * This function initializes a new graph structure with default capacities for vertices, edges, and
 * Horton cycles. It allocates memory for the vertices array, edges array, neighbor lists, and
 * Horton cycles array. The edge_indices, predecessors, minimal_basis and distances matrices are initialized to
 * NULL.
 *
 * @return Pointer to the newly created graph structure.
 */
Graph *create_graph() {
    Graph *graph = calloc(1, sizeof(Graph));
    if (graph == NULL) {
        perror("Failed to allocate memory for Graph in create_graph");
        exit(1);
    }

    graph->capacity_vertices = 4;
    graph->vertices = malloc(graph->capacity_vertices * sizeof(Vertex));
    if (graph->vertices == NULL) {
        free(graph);
        perror("Failed to allocate memory for Graph.vertices in create_graph");
        exit(1);
    }
    graph->neighbors = malloc(graph->capacity_vertices * sizeof(Neighbor_list));
    if (graph->neighbors == NULL) {
        free(graph->vertices);
        free(graph);
        perror("Failed to allocate memory for Graph.neighbors in create_graph");
        exit(1);
    }

    for (int i = 0; i < graph->capacity_vertices; i++) {
        graph->neighbors[i].neighbors = NULL;
        graph->neighbors[i].count = 0;
        graph->neighbors[i].capacity = 0;
    }

    graph->capacity_edges = 4;
    graph->edges = malloc(graph->capacity_edges * sizeof(Edge));
    if (graph->edges == NULL) {
        free(graph->neighbors);
        free(graph->vertices);
        free(graph);
        perror("Failed to allocate memory for Graph.edges in create_graph");
        exit(1);
    }

    graph->capacity_horton_cycles = 4;
    graph->horton_cycles = malloc(graph->capacity_horton_cycles * sizeof(Path));
    if (graph->horton_cycles == NULL) {
        free(graph->edges);
        free(graph->neighbors);
        free(graph->vertices);
        free(graph);
        perror("Failed to allocate memory for Graph.horton_cycles in create_graph");
        exit(1);
    }

    graph->edge_indices = NULL;
    graph->predecessors = NULL;
    graph->distances = NULL;

    graph->nb_minimal_bases = 0;
    graph->minimals_basis = NULL;
    graph->basis_dimension = 0;

    return graph;
}

/**
 * @brief Compute all shortest paths in the graph using BFS.
 *
 * For each source vertex, a BFS is run to fill the distances and predecessors matrices. The
 * distances and predecessor matrix should be filled with -1's before calling this function.
 * The graph is considered unweighted: the primary criterion is the number of edges. When two
 * paths to a vertex have the same distance, the tiebreak is resolved by comparing the
 * cumulative sum of edge labels along each path. We keep the predecessor whose path yields the
 * smallest total label sum. This makes the BFS deterministic and sensitive to the current
 * edge-label permutation: different permutations (as produced by multiple_horton) yield
 * different BFS trees and therefore potentially different minimal cycle bases.
 *
 * @param g Pointer to the graph structure for which to compute the shortest paths.
 */
void compute_all_shortest_paths(const Graph *g) {

    // Creation of the array distance_label to store the cumulative sum of edge labels for the path
    // from the source to each vertex
    int *distance_label = malloc(g->nb_vertex * sizeof(int));
    if (!distance_label) {
        perror("Failed to allocate distance_label in compute_all_shortest_paths");
        exit(1);
    }

    // FIFO queue for the BFS
    int *queue = malloc(g->nb_vertex * sizeof(int));
    if (!queue) {
        free(distance_label);
        perror("Failed to allocate queue in compute_all_shortest_paths");
        exit(1);
    }

    // For each vertex s as the source of the BFS
    for (int s = 0; s < g->nb_vertex; s++) {

        // Reset the distance_label array for the new source vertex s, setting all distances to
        // infinity
        for (int i = 0; i < g->nb_vertex; i++) {
            distance_label[i] = INT_MAX;
        }
        distance_label[s] = 0; // The distance label of the source vertex to itself is 0

        int front = 0, rear = 1; // Pointers for the front and rear of the queue
        queue[front] = s;
        g->distances[(s * g->nb_vertex) + s] = 0; // Distance from s to itself is 0
        g->predecessors[s * g->nb_vertex + s] = s;

        // While the queue is not empty, we process the vertices in FIFO order
        while (front < rear) {
            const int u =  queue[front]; // Get the vertex at the front of the queue
            front++;

            // For each neighbor v of u, we check if we can improve the distance to v through u,
            // and if so, we update the distance, the distance_label, and the predecessor of v
            for (int i = 0; i < g->neighbors[u].count; i++) {
                const int v = g->neighbors[u].neighbors[i];
                const int edge_index = g->edge_indices[(u * g->nb_vertex) + v];

                // Get the label of the edge (u,v)
                const int current_edge_label = g->edges[edge_index].label;

                // If v has not been visited yet
                if (g->distances[s * g->nb_vertex + v] == -1) {
                    // Update the distance from s to v through u (which is the distance from s to u
                    // plus 1 edge)
                    g->distances[s * g->nb_vertex + v] = g->distances[s * g->nb_vertex + u] + 1;

                    // Update the distance_label of v to be the distance_label of u plus the label
                    // of the edge (u,v)
                    distance_label[v] = distance_label[u] + current_edge_label;
                    g->predecessors[s * g->nb_vertex + v] = u;

                    // Add v to the rear of the queue to be processed later
                    queue[rear] = v;
                    rear++;
                }
                if (g->distances[s * g->nb_vertex + v] == -1) {
                    g->distances[s * g->nb_vertex + v] = g->distances[s * g->nb_vertex + u] + 1;
                    distance_label[v] = distance_label[u] + current_edge_label;
                    g->predecessors[s * g->nb_vertex + v] = u;
                    queue[rear++] = v;
                }
            }
        }
    }

    // Free the memory allocated for the queue and distance_label arrays
    free(queue);
    free(distance_label);

}

Path create_path(const Graph *g, const int u, const int v) {
    if (g->predecessors == NULL || g->edges == NULL) {
        perror(
            "The function create_path was called but the shortest paths have not been "
            "computed");
        exit(1);
    }

    Path p;
    p.length = 0;
    p.edges_ids = calloc(g->nb_edges, sizeof(uint32_t)); // Initialize the
    // edges_ids array with 0s
    if (p.edges_ids == NULL) {
        perror("Failed to allocate memory for edges_ids in create_path");
        exit(1);
    }

    p.vertices_ids = calloc(g->nb_vertex, sizeof(uint32_t));
    if (p.vertices_ids == NULL) {
        free(p.edges_ids);
        perror("Failed to allocate memory for vertices_ids in create_path");
        exit(1);
    }

    int current = v;
    p.vertices_ids[v] = 1;
    p.vertices_ids[u] = 1;

    while (current != u) {
        const int predecessor = g->predecessors[u * g->nb_vertex + current];
        if (predecessor == -1) {
            free(p.edges_ids);
            free(p.vertices_ids);
            fprintf(
                stderr,
                "No predecessor for %d in the path from %d to %d in create_path, should not "
                "happen if the function is called correctly",
                current, u, v);
            exit(1);
        }

        const int edge_index = g->edge_indices[current * g->nb_vertex + predecessor];
        if (edge_index == -1) {
            free(p.edges_ids);
            free(p.vertices_ids);
            fprintf(
                stderr,
                "No edge between %d and its predecessor %d in create_path, should not happen "
                "if the function is called correctly",
                current, predecessor);
            exit(1);
        }

        p.edges_ids[edge_index] = 1; // Add the edge to the path
        p.vertices_ids[g->edges[edge_index].u] = 1;
        p.vertices_ids[g->edges[edge_index].v] = 1;
        p.length++;
        current = predecessor; // Move to the predecessor for the next iteration
    }
    return p;
}

/**
 * @brief Add a cycle to the list of Horton cycles in the graph.
 *
 * This function adds a given cycle to the horton_cycles array in the graph structure. If the
 * current number of Horton cycles reaches the capacity of the array, it doubles the capacity and
 * reallocates memory for the array before adding the new cycle.
 *
 * @param g Pointer to the graph structure where the cycle will be added.
 * @param c The cycle to be added to the list of Horton cycles.
 */
void add_horton_cycle(Graph *g, const Path c) {
    // Double the capacity of the horton_cycles array if the current capacity is reached and reallocate memory
    if (g->nb_horton_cycles == g->capacity_horton_cycles) {
        const int new_cap = g->capacity_horton_cycles * 2; // double the capacity
        Path *tmp = realloc(g->horton_cycles, new_cap * sizeof(Path));
        if (!tmp) {
            perror("Failed to reallocate memory for Graph.horton_cycles in add_horton_cycle");
            exit(1);
        }
        g->horton_cycles = tmp;
        g->capacity_horton_cycles = new_cap;
    }

    g->horton_cycles[g->nb_horton_cycles++] = c;
}

/**
 * @brief Find all cycles generated by Horton's algorithm and store them in the graph structure.
 *
 * This function implements the core of Horton's algorithm to find all cycles that can be formed by
 * combining the shortest paths from a root vertex to two vertices u and v with the edge (u,v).
 * It iterates over all vertices as potential roots and all edges to find valid combinations that
 * form cycles. The generated cycles are stored in the horton_cycles array of the graph structure.
 *
 * @param g Pointer to the graph structure for which to find Horton cycles.
 */
void find_horton_cycles(Graph *g) {
    if (g->predecessors == NULL || g->edges == NULL) {
        printf("find_horton_cycles called but the shortest paths have not been computed or "
            "there is no edge in the graph\n");
        return;
    }

    // For each vertex r
    for (int v = 0; v < g->nb_vertex; v++) {

        // For each edge (u,v)
        for (int e = 0; e < g->nb_edges; e++) {
            const int x = g->edges[e].u;
            const int y = g->edges[e].v;

            if (x == v || y == v) continue;

            const int pred_x = g->predecessors[v * g->nb_vertex + x]; // predecessor of u in the
            // shortest path from r to u
            const int pred_y = g->predecessors[v * g->nb_vertex + y]; // predecessor of v in the
            // shortest path from r to v

            // If there is no path from r to u or from r to v, we cannot create a cycle with the
            // edge (u,v)
            if (pred_x == -1 || pred_y == -1) {
                fprintf(stderr, "No path from %d to %d or from %d to %d, skipping edge (%d,%d)"
                                " for root %d\n", v, x, v, y, x, y, v);
                continue;
            }

            const Path path_v_to_x = create_path(g, v, x);
            const Path path_v_to_y = create_path(g, v, y);

            uint32_t *expected = calloc(g->nb_vertex, sizeof(uint32_t));
            if (!expected) {
                perror("Failed to allocate memory for expected in find_horton_cycles");
                exit(1);
            }
            expected[v] = 1;

            uint32_t *actual = calloc(g->nb_vertex, sizeof(uint32_t));
            if (!actual) {
                free(expected);
                perror("Failed to allocate memory for actual in find_horton_cycles");
                exit(1);
            }

            for (int i = 0; i < g->nb_vertex; i++) {
                actual[i] = path_v_to_x.vertices_ids[i] && path_v_to_y.vertices_ids[i];
            }

            if (memcmp(expected, actual, g->nb_vertex * sizeof(uint32_t)) == 0) {
                Path c;
                uint32_t *edges_of_c = calloc(g->nb_edges, sizeof(uint32_t));
                uint32_t *vertices_of_c = calloc(g->nb_vertex, sizeof(uint32_t));
                if (!edges_of_c) {
                    free(expected);
                    free(actual);
                    perror("Failed to allocate memory for edges_of_c in find_horton_cycles");
                    exit(1);
                }
                if (!vertices_of_c) {
                    free(expected);
                    free(actual);
                    free(edges_of_c);
                    perror("Failed to allocate memory for vertices_of_c in find_horton_cycles");
                    exit(1);
                }
                for (int i = 0; i < g->nb_edges; i++) {
                    edges_of_c[i] = path_v_to_x.edges_ids[i] ^ path_v_to_y.edges_ids[i];
                }
                for (int i = 0; i < g->nb_vertex; i++) {
                    vertices_of_c[i] = path_v_to_x.vertices_ids[i] ^ path_v_to_y.vertices_ids[i];
                }
                c.edges_ids = edges_of_c;
                c.edges_ids[e] = 1;
                c.vertices_ids = vertices_of_c;
                c.length = path_v_to_x.length + path_v_to_y.length + 1;
                add_horton_cycle(g, c);
            }
        }
    }
}

// Global variable to store the number of edges for the qsort comparison function
static int nb_edges_for_qsort = 0;

/**
 * @brief Comparison function for sorting cycles by weight.
 *
 * This function is used as a comparison function for qsort to sort the cycles generated by
 * Horton's algorithm by their length.
 *
 * @param a Pointer to the first cycle to compare (cast to void* for qsort).
 * @param b Pointer to the second cycle to compare (cast to void* for qsort).
 * @return Negative value if c1 < c2, positive value if c1 > c2, 0 if c1 == c2 for sorting purposes.
 */
int compare_cycles(const void *a, const void *b) {
    const Path *c1 = a;
    const Path *c2 = b;

    // First criterion: compare by length
    if (c1->length != c2->length) return c1->length - c2->length;

    // Second criterion: compare the edges_ids arrays for tie-breaking (lexicographical order) to
    // ensure a deterministic order of cycles with the same length
    return memcmp(c1->edges_ids, c2->edges_ids,
                  nb_edges_for_qsort * sizeof(uint32_t));
}

/** @brief Sort the cycles generated by Horton's algorithm in the graph structure.
 *
 * This function sorts the horton_cycles array in the graph structure using qsort and the
 * compare_cycles function to order the cycles by their length and then by their weight.
 *
 * @param g Pointer to the graph structure whose Horton cycles will be sorted.
 */
void order_horton_cycles(const Graph *g) {
    nb_edges_for_qsort = g->nb_edges;
    if (g->nb_horton_cycles > 1) {
        qsort(g->horton_cycles,
              g->nb_horton_cycles,
              sizeof(Path),
              compare_cycles);
    }
}

/**
 * @brief Check if two cycles are the same by comparing their edges_ids arrays.
 *
 * @param c1 First cycle to compare.
 * @param c2 Second cycle to compare.
 * @param nb_edges Number of edges in the graph
 * @return 1 if the two cycles are the same, 0 otherwise.
 */
int are_same_cycles(const Path *c1, const Path *c2, const int nb_edges) {
    if (c1->length != c2->length) return 0;
    return memcmp(c1->edges_ids, c2->edges_ids, nb_edges * sizeof(uint32_t)) == 0;
}

/**
 * @brief Wrapper for the are_same_cycles function to be used as a comparison function for qsort.
 *
 * @param c1 First cycle to compare (cast to void* for qsort).
 * @param c2 Second cycle to compare (cast to void* for qsort).
 * @return 1 if the two cycles are the same, 0 otherwise.
 */
int are_same_cycles_qsort(const void *c1, const void *c2) {
    return are_same_cycles(c1, c2, nb_edges_for_qsort);
}


/**
 * @brief Check if two minimal cycle bases are the same by comparing their cycles.
 *
 * This function checks if two minimal cycle bases are the same by comparing their sorted cycles.
 * If all the cycles are the same (they have the same edges_ids), then the two bases are considered
 * the same.
 *
 * @param a First minimal basis to compare (cast to void* for qsort).
 * @param b Second minimal basis to compare (cast to void* for qsort).
 * @param nb_edges Number of edges in the graph.
 * @param basis_dimension Dimension of the cycle space of the graph (i.e., the number of cycles in
 * any cycle basis).
 * @return 1 if the two minimal bases are the same, 0 otherwise.
 */
int are_same_minimal_bases(const void *a, const void *b, const int nb_edges,
                           const int basis_dimension) {
    const Minimal_basis *mb1 = a;
    const Minimal_basis *mb2 = b;

    // Create temporary sorted copies of both bases
    Path *sorted1 = malloc(basis_dimension * sizeof(Path));
    Path *sorted2 = malloc(basis_dimension * sizeof(Path));
    if (!sorted1 || !sorted2) {
        free(sorted1);
        free(sorted2);
        perror("Failed to allocate memory for sorted1 or sorted2 in are_same_minimal_bases");
        exit(1);
    }

    // Copy and allocate edges_ids for sorted1
    for (int i = 0; i < basis_dimension; i++) {
        sorted1[i].length = mb1->cycles[i].length;
        sorted1[i].edges_ids = malloc(nb_edges * sizeof(uint32_t));
        if (!sorted1[i].edges_ids) {
            for (int j = 0; j < i; j++) free(sorted1[j].edges_ids);
            free(sorted1);
            free(sorted2);
            perror("Failed to allocate memory for sorted1.edges_ids in are_same_minimal_bases");
            exit(1);
        }
        memcpy(sorted1[i].edges_ids, mb1->cycles[i].edges_ids, nb_edges * sizeof(uint32_t));
    }

    // Copy and allocate edges_ids for sorted2
    for (int i = 0; i < basis_dimension; i++) {
        sorted2[i].length = mb2->cycles[i].length;
        sorted2[i].edges_ids = malloc(nb_edges * sizeof(uint32_t));
        if (!sorted2[i].edges_ids) {
            for (int j = 0; j < i; j++) free(sorted2[j].edges_ids);
            for (int j = 0; j < basis_dimension; j++) free(sorted1[j].edges_ids);
            free(sorted1);
            free(sorted2);
            perror("Failed to allocate memory for sorted2.edges_ids in are_same_minimal_bases");
            exit(1);
        }
        memcpy(sorted2[i].edges_ids, mb2->cycles[i].edges_ids, nb_edges * sizeof(uint32_t));
    }

    nb_edges_for_qsort = nb_edges; // Set the global variable for qsort comparison
    // Sort the cycles in both bases by their edges_ids
    qsort(sorted1, basis_dimension, sizeof(Path), compare_cycles);
    qsort(sorted2, basis_dimension, sizeof(Path), compare_cycles);

    // Compare sorted cycles
    int result = 1;
    for (int i = 0; i < basis_dimension; i++) {
        if (!are_same_cycles(&sorted1[i], &sorted2[i], nb_edges)) {
            // Print the difference between the two cycles
            fprintf(stderr, "Difference found in cycle %d:\n", i);
            fprintf(stderr, "Cycle in mb1: length = %d, edges = [", sorted1[i].length);
            for (int e = 0; e < nb_edges; e++) {
                if (sorted1[i].edges_ids[e] == 1) {
                    fprintf(stderr, " %d", e);
                }
            }
            fprintf(stderr, " ]\n");
            fprintf(stderr, "Cycle in mb2: length = %d, edges = [", sorted2[i].length);
            for (int e = 0; e < nb_edges; e++) {
                if (sorted2[i].edges_ids[e] == 1) {
                    fprintf(stderr, " %d", e);
                }
            }
            fprintf(stderr, " ]\n");
            result = 0;
            break;
        }
    }

    for (int i = 0; i < basis_dimension; i++) {
        free(sorted1[i].edges_ids);
        free(sorted2[i].edges_ids);
    }
    free(sorted1);
    free(sorted2);

    return result;
}

/**
 * @brief Check if a candidate minimal basis already exists in the list of minimal bases in the graph.
 *
 * This function checks if a given candidate minimal basis is already present in the list of minimal
 * bases stored in the graph structure. It compares the candidate basis with each existing basis using
 * the are_same_minimal_bases function.
 *
 * @param g Pointer to the graph structure containing the list of minimal bases.
 * @param candidate The candidate minimal basis to check for existence in the list.
 * @return 1 if the candidate basis already exists in the list, 0 otherwise.
 */
int already_exists_basis(const Graph *g, const Minimal_basis candidate) {
    if (g->minimals_basis == NULL || g->nb_minimal_bases == 0) {
        return 0;
    }

    // For each minimal basis of the graph
    for (int b = 0; b < g->nb_minimal_bases; b++) {
        const Minimal_basis mb = g->minimals_basis[b];
        if (are_same_minimal_bases(&mb, &candidate, g->nb_edges, g->basis_dimension)) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Compute a minimal cycle basis of the graph using a greedy cycle cover algorithm.
 *
 * This function computes a minimal cycle basis of the graph by iterating over the list of
 * Horton cycles (sorted by weight) and adding each cycle to the basis if it is linearly independent
 * of the cycles already in the basis. The linear independence is checked using a Gaussian
 * elimination-like process on the binary vectors representing the cycles.
 *
 * @param g Pointer to the graph structure.
 */
void greedy_cycle_cover(Graph *g) {
    if (g->horton_cycles == NULL || g->nb_horton_cycles == 0) {
        g->basis_dimension = 0;
        printf("No cycles to compute the minimal cycle basis in greedy_cycle_cover\n");
        return;
    }

    // Basis[e] will point to the reduced vector of the cycle in the base that has a pivot at edge e,
    // or NULL if no cycle in the base has a pivot at e
    uint32_t **basis = calloc(g->nb_edges, sizeof(uint32_t *));
    if (basis == NULL) {
        perror("Failed to allocate memory for basis in greedy_cycle_cover");
        exit(1);
    }

    // Array to store the cycles of the minimal cycle basis
    Path *mcb = calloc(g->nb_horton_cycles, sizeof(Path));
    if (mcb == NULL) {
        perror("Failed to allocate memory for mcb in greedy_cycle_cover");
        exit(1);
    }

    // Counter for the number of cycles in the minimal cycle basis
    int mcb_count = 0;

    // For each cycle in the sorted list of Horton cycles
    for (int i = 0; i < g->nb_horton_cycles; i++) {
        const Path candidate = g->horton_cycles[i];

        // Create a copy of the candidate cycle's edge vector to perform the reduction process
        uint32_t *vec = malloc(g->nb_edges * sizeof(uint32_t));
        if (!vec) {
            perror("Failed to allocate memory for vec in greedy_cycle_cover");
            exit(1);
        }
        memcpy(vec, candidate.edges_ids, g->nb_edges * sizeof(uint32_t));

        int independent = 0;

        // Gaussian elimination to check if the cycle is linearly independent of the current basis
        for (int e = 0; e < g->nb_edges; e++) {
            if (vec[e] == 1) {
                // If there is no cycle in the basis with a pivot at edge e, then vec is independent
                // of the current basis
                if (basis[e] == NULL) {
                    basis[e] = vec;
                    independent = 1;
                    break;
                }
                // Else there is a cycle in the basis with a pivot at edge e, we add it to vec (XOR)
                // to see if we can reduce vec to 0, which would mean that the candidate cycle is
                // not independent of the current basis
                for (int k = 0; k < g->nb_edges; k++) {
                    vec[k] ^= basis[e][k];
                }
            }
        }

        // If the candidat cycle is independent, we add it to the minimal cycle basis
        if (independent) {
            mcb[mcb_count].length = candidate.length;
            mcb[mcb_count].edges_ids = malloc(g->nb_edges * sizeof(uint32_t));

            if (mcb[mcb_count].edges_ids == NULL) {
                perror("Failed to allocate memory for edges_ids in greedy_cycle_cover");
                exit(1);
            }
            memcpy(mcb[mcb_count].edges_ids, candidate.edges_ids, g->nb_edges * sizeof(uint32_t));
            mcb_count++;
        } else {
            free(vec);
        }
    }

    for (int e = 0; e < g->nb_edges; e++) {
        if (basis[e] != NULL) {
            free(basis[e]);
        }
    }
    free(basis);

    if (mcb_count > 0) {
        Path *exact_mcb = realloc(mcb, mcb_count * sizeof(Path));
        if (exact_mcb) {
            mcb = exact_mcb;
        } else {
            perror("Failed to reallocate memory for mcb in greedy_cycle_cover");
            exit(1);
        }
    }

    const Minimal_basis mb = {.cycles = mcb};

    if (g->nb_minimal_bases != 0) {
        if (mcb_count != g->basis_dimension) {
            perror("Error in the basis dimension, two different dimensions for the cycle "
                "space of the graph, should not happen");
            exit(1);
        }
        if (already_exists_basis(g, mb)) {
            for (int k = 0; k < mcb_count; k++) {
                free(mcb[k].edges_ids);
            }
            free(mcb);
            return;
        }
    }

    g->nb_minimal_bases += 1;
    Minimal_basis *tmp_bases = realloc(g->minimals_basis,
                                       g->nb_minimal_bases * sizeof(Minimal_basis));
    if (tmp_bases == NULL) {
        perror("Failed to reallocate memory for Graph.minimals_basis in greedy_cycle_cover");
        exit(1);
    }
    g->minimals_basis = tmp_bases;
    g->minimals_basis[g->nb_minimal_bases - 1] = mb;
    g->basis_dimension = mcb_count;
}

/**
 * @brief Save the graph structure to a file in a specific format.
 *
 * This function saves the graph structure, including its vertices, edges, and minimal cycle bases,
 * to a file in a specific format. The first line of the file contains the number of vertices V,
 * number of edges E, number of minimal bases M, and the dimension of the cycle space D. The V next
 * lines contain the vertices with their ID and coordinates. The E next lines contain the edges with
 * their id, endpoints, and weight. The M*D remaining lines contain the cycles of each minimal basis
 * with the indices of the edges that belong to each cycle.
 *
 * @param g Pointer to the graph structure to be saved.
 * @param filename The name of the file where the graph will be saved.
 */
void save_graph(const Graph *g, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Failed to open file for writing in save_graph");
        exit(1);
    }

    if (g == NULL || g->vertices == NULL || g->edges == NULL) {
        fclose(file);
        perror("Invalid graph structure in save_graph");
        exit(1);
    }

    if (g->minimals_basis) {
        fprintf(file, "%d %d %d %d\n", g->nb_vertex, g->nb_edges, g->nb_minimal_bases,
                g->basis_dimension);
    } else {
        fprintf(file, "%d %d %d %d\n", g->nb_vertex, g->nb_edges, 0, 0);
    }

    // Write the vertices information
    for (int i = 0; i < g->nb_vertex; i++) {
        fprintf(file, "%d %lf %lf\n", g->vertices[i].id, g->vertices[i].x, g->vertices[i].y);
    }

    // Write the edges information
    for (int i = 0; i < g->nb_edges; i++) {
        fprintf(file, "%d %d %d %d %d\n", g->edges[i].id, g->edges[i].u, g->edges[i].v,
                g->edges[i].weight, g->edges[i].label);
    }

    if (g->minimals_basis) {
        // For each minimal basis
        for (int b = 0; b < g->nb_minimal_bases; b++) {
            // For each cycle in the minimal basis
            for (int i = 0; i < g->basis_dimension; i++) {
                // Write the index of the edge in the graph's edges array
                for (int e = 0; e < g->nb_edges; e++) {
                    if (g->minimals_basis[b].cycles[i].edges_ids[e] == 1) {
                        fprintf(file, "%d ", e);
                    }
                }
                fprintf(file, "\n");
            }
        }
    }
    fclose(file);
}

/**
 * @brief Load a graph structure from a file in a specific format. This function only loads the
 * vertices and edges of the graph, and does not load the minimal cycle bases. However, the first
 * line of the file should still contain two numbers that won't be used.
 *
 * This function loads a graph structure from a file that is formatted in a specific way. The first
 * line of the file should contain the number of vertices V, number of edges E, number of minimal
 * bases M, and the dimension of the cycle space D. The next V lines should contain the vertices with
 * their ID and coordinates. The next E lines should contain the edges with their id, endpoints, and
 * weight. The remaining lines should contain the cycles of each minimal basis with the indices of
 * the edges that belong to each cycle.
 *
 * @param g Pointer to the graph structure where the loaded graph will be stored.
 * @param filename The name of the file from which to load the graph.
 */
void load_graph(Graph *g, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Failed to open file for reading in load_graph");
        exit(1);
    }

    int nb_vertex, nb_edges, nb_minimal_bases, basis_dimension;
    if (fscanf(file, "%d %d %d %d\n", &nb_vertex, &nb_edges, &nb_minimal_bases,
               &basis_dimension) != 4) {
        fclose(file);
        perror("Failed to read graph metadata in load_graph");
        exit(1);
    }

    for (int i = 0; i < nb_vertex; i++) {
        int id;
        double x, y;
        if (fscanf(file, "%d %lf %lf\n", &id, &x, &y) != 3) {
            fclose(file);
            perror("Failed to read vertex information in load_graph");
            exit(1);
        }
        create_vertex(g, x, y);
    }

    for (int i = 0; i < nb_edges; i++) {
        int id, u, v, weight;
        if (fscanf(file, "%d %d %d %d\n", &id, &u, &v, &weight) != 4) {
            fclose(file);
            perror("Failed to read edge information in load_graph");
            exit(1);
        }
        create_edge(g, u, v, weight);
    }

    fclose(file);
}

/**
 * @brief Prepare the graph structure for Horton's algorithm by initializing the edge_indices,
 * predecessors, and distances matrices.
 *
 * This function initializes the edge_indices matrix to store the index of each edge in the graph's
 * edges array, and initializes the predecessors and distances matrices to store the results of
 * shortest path computations. The edge_indices matrix is filled based on the edges present in the
 * graph, while the predecessors and distances matrices are initialized to -1 to indicate that
 * shortest paths have not been computed yet.
 *
 * @param graph Pointer to the graph structure to be prepared for Horton's algorithm.
 */
void prepare_graph_matrices(Graph *graph) {
    if (graph->edges == NULL || graph->nb_vertex == 0) {
        printf(
            "The graph has no edges or vertices, skipping the preparation of matrices in prepare_graph_matrices\n");
        return;
    }

    const int n = graph->nb_vertex;
    const int size = n * n;

    // Free existing matrices
    if (graph->edge_indices) free(graph->edge_indices);
    if (graph->predecessors) free(graph->predecessors);
    if (graph->distances) free(graph->distances);

    graph->edge_indices = malloc(size * sizeof(int));
    graph->predecessors = malloc(size * sizeof(int));
    graph->distances = malloc(size * sizeof(int));
    if (!graph->edge_indices || !graph->predecessors || !graph->distances) {
        free(graph->edge_indices);
        free(graph->predecessors);
        free(graph->distances);
        perror("Failed to allocate memory for edge_indices, predecessors, or distances in "
            "prepare_graph_matrices");
        exit(1);
    }

    // Initialize the edge_indices, predecessors, and distances matrices with -1
    for (int i = 0; i < size; i++) {
        graph->edge_indices[i] = -1;
        graph->predecessors[i] = -1;
        graph->distances[i] = -1;
    }

    // Fill the edge_indices matrix
    for (int i = 0; i < graph->nb_edges; i++) {
        const uint32_t u = graph->edges[i].u;
        const uint32_t v = graph->edges[i].v;
        graph->edge_indices[u * n + v] = i;
        graph->edge_indices[v * n + u] = i;
    }
}

/**
 * @brief Reset the results of Horton's algorithm in the graph structure by freeing the memory
 * allocated for Horton cycles and resetting the count of Horton cycles to zero.
 *
 * This function is used to clear the results of a previous execution of Horton's algorithm before
 * running it again, ensuring that the graph structure is clean and ready for new computations.
 *
 * @param g Pointer to the graph structure whose Horton cycle results will be reset.
 */
void reset_graph_results(Graph *g) {
    if (g->horton_cycles) {
        for (int i = 0; i < g->nb_horton_cycles; i++) {
            free(g->horton_cycles[i].edges_ids);
        }
        g->nb_horton_cycles = 0;
    }
}

/**
 * @brief Main function to execute Horton's algorithm on a graph structure.
 *
 * This function orchestrates the execution of Horton's algorithm by first resetting any previous
 * results, preparing the necessary matrices for shortest path computations, computing all
 * shortest paths using BFS, finding all cycles, ordering them, and finally computing a
 * minimal cycle basis using a greedy cycle cover algorithm.
 *
 * @param g Pointer to the graph structure on which to execute Horton's algorithm.
 */
void horton(Graph *g) {
    reset_graph_results(g);
    prepare_graph_matrices(g);
    compute_all_shortest_paths(g);
    find_horton_cycles(g);
    order_horton_cycles(g);
    greedy_cycle_cover(g);
}

/**
 * @brief Execute Horton's algorithm multiple times on the same graph structure with different
 * permutations of edge label.
 *
 * This function allows for executing Horton's algorithm multiple times on the same graph structure
 * by applying different permutations of edge label in order to try to find different minimal
 * cycle bases. The function continues to generate and apply permutations until it reaches the
 * specified maximum number of permutations or exhausts all possible permutations.
 *
 * @param g Pointer to the graph structure on which to execute Horton's algorithm multiple times.
 * @param inv Pointer to an array representing the inversion table for generating permutations of
 * edge label. It should be an array of size equal to the number of edges in the graph, initialized
 * to represent the identity permutation (inv[i] = 0 for all i in [0,...,nb_edges-1]).
 * @param max_permutations Maximum number of permutations to apply and execute Horton's algorithm.
 */
void multiple_horton(Graph *g, int *inv, const int max_permutations) {
    if (inv == NULL || g->edges == NULL) {
        printf("No inversion table given or no edges in the graph, skipping multiple_horton\n");
        return;
    }

    int *perm = calloc(g->nb_edges, sizeof(int));
    if (!perm) {
        perror("Failed to allocate perm in multiple_horton");
        exit(1);
    }

    generate_random_inversion_table(inv, g->nb_edges);
    inversion_to_permutation(inv, perm, g->nb_edges);

    // Apply labels from first permutation
    for (int i = 0; i < g->nb_edges; i++) {
        g->edges[i].label = perm[i];
    }
    horton(g);

    int cnt = 1; // count the number of permutations processed
    while (cnt < max_permutations) {

        memset(perm, 0, g->nb_edges * sizeof(int));
        generate_random_inversion_table(inv, g->nb_edges);
        inversion_to_permutation(inv, perm, g->nb_edges);

        for (int i = 0; i < g->nb_edges; i++) {
            g->edges[i].label = perm[i]; // update the edge label according to the new permutation
        }

        printf("%d / %d\n", cnt++, max_permutations);
        horton(g);
    }

    free(perm);
}

/**
 * @brief Free the memory allocated for the graph structure and its components.
 *
 * This function is responsible for freeing all the memory allocated for the graph structure.
 *
 * @param g Pointer to the graph structure to be freed.
 */
void free_graph(Graph *g) {
    if (!g) return;

    free(g->vertices);
    free(g->edges);

    if (g->neighbors) {
        for (int i = 0; i < g->nb_vertex; i++) {
            free(g->neighbors[i].neighbors);
        }
        free(g->neighbors);
    }

    free(g->edge_indices);
    free(g->predecessors);
    free(g->distances);

    if (g->horton_cycles) {
        for (int i = 0; i < g->nb_horton_cycles; i++) {
            free(g->horton_cycles[i].edges_ids);
        }
        free(g->horton_cycles);
    }

    if (g->minimals_basis) {
        for (int b = 0; b < g->nb_minimal_bases; b++) {
            for (int i = 0; i < g->basis_dimension; i++) {
                free(g->minimals_basis[b].cycles[i].edges_ids);
            }
            free(g->minimals_basis[b].cycles);
        }
        free(g->minimals_basis);
    }
    free(g);
}
