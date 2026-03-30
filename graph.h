#ifndef TER_GRAPHES_PLANAIRES_GRAPH_H
#define TER_GRAPHES_PLANAIRES_GRAPH_H
#include <stdint.h>

/**
 * @brief Vertex of the graph.
 *
 * Represents a node in the graph. The coordinates (x, y) are optional and are only used for
 * visualization purposes. The label is currently not used in the algorithm.
 */
typedef struct {
    double x; /**< X coordinate in the plane. */
    double y; /**< Y coordinate in the plane. */
    int id; /**< Unique identifier (0 to nb_vertex-1), represents the index of the vertex in
    graph->vertices. */
    int label; /**< Label of the vertex (0 to nb_vertex-1). */
    int deleted; /**< Flag to indicate if the vertex has been deleted (1 if deleted, 0 otherwise). */
} Vertex;

/**
 * @brief Edge of the graph.
 *
 * Represents a connection between two vertices in the graph. The weight is currently not used in
 * the algorithm but can be useful for future extensions (e.g., weighted cycle bases). The label is
 * used to assign a unique identifier to each edge for the purpose of generating different minimal
 * cycle bases by permuting these labels.
 */
typedef struct {
    int u; /**< ID of the first vertex */
    int v; /**< ID of the second vertex */
    int label; /**< Label of the edge (0 to nb_edges-1), used for generating different minimal cycle
    bases by permuting these labels. */
    int id; /**< Unique identifier (0 to nb_edges-1), represents the index of the edge in
    graph->edges */
    int deleted;
    int is_outer; /**< Flag to indicate if the edge is an outer edge (1 if it is an outer edge, 0 otherwise). */
    int face_id[2];
} Edge;

/**
 * @brief List of neighbors for a vertex.
 *
 * Stores the list of neighbors of a vertex in the graph. It contains a pointer to a dynamic array
 * of neighbors id, the count of neighbors currently stored, and the capacity of the allocated
 * array.
 */
typedef struct {
    int *neighbors; /**< Dynamic array of neighbors id, where neighbors[i] is the id of the (i+1)th
    neighbor vertex in graph->vertices and the ith neighbor in polar angle */
    double *angles; /**< Dynamic array of angles to neighbors, where angles[i] is the angle between
    the edge connecting the vertex to its ith neighbor */
    int count; /**< Number of neighbors currently stored in the list */
    int capacity; /**< Capacity of the allocated neighbors array */
} Neighbor_list;

/**
 * @brief Path in the graph.
 *
 * Represents a Path in the graph. The edges_ids array is a binary vector of size nb_edges, where
 * edges_ids[i] = 1 if the edge with index i in graph->edges is part of the path, and 0 otherwise.
 * The length is the number of vertices in the path.
 */
typedef struct {
    uint32_t *edges_ids;
    /**< Binary vector of size nb_edges, where edges_ids[i] = 1 if the edge with index i in
     *graph->edges is part of the path, and 0 otherwise */
    uint32_t *vertices_ids; /**< Binary vector of size nb_vertex, where vertices_ids[i] = 1 if the
    vertex with id i is in the path. */
    uint32_t *edges_labels; /**< Binary vector of size nb_edges, where edges_labels[i] = 1 if the
    edge with label i is part of the path, and 0 otherwise. */
    int length; /**< Number of vertices in the path */
} Path;

/**
 * @brief Minimal cycle basis of the graph.
 *
 * Represents the cycles of a minimal cycle basis of the graph.
 */
typedef struct {
    Path *cycles; /**< Array of cycles in the minimal cycle basis, where cycles[i] is the (i+1)-th
    cycle of the basis */
    int is_faces; /**< 1 if this basis correspond to all faces (except the outer one), 0 otherwise */
    int is_faces_outer; /**< 1 if this basis correspond to all faces (including the outer one) except one, 0 otherwise */
} Minimal_basis;

/**
 * @brief Graph structure.
 *
 * Represents the graph with its vertices, edges, and the data structures needed for Horton's
 * algorithm and the storage of the differents minimal cycle bases.
 */
typedef struct {
    int nb_vertex; /**< Number of vertices in the graph */
    int capacity_vertices; /**< Capacity of the allocated vertices array */
    Vertex *vertices; /**< Dynamic array of vertices, where vertices[i] is the vertex with id i in
    the graph */

    int nb_edges; /**< Number of edges in the graph */
    int capacity_edges; /**< Capacity of the allocated edges array */
    Edge *edges; /**< Dynamic array of edges, where edges[i] is the edge with id i in the graph */

    int *edge_indices; /**< Dynamic array to store the ids of edges between pairs of vertices.
    edge_indices[(u * nb_vertex) + v] gives the index of the edge connecting vertex u and vertex v
    in graph->edges, or -1 if there is no edge between u and v. */

    Neighbor_list *neighbors; /**< Dynamic array of neighbor lists, where neighbors[i] is the list
    of neighbors of the vertex with id i in the graph */

    int *predecessors; /**< Flattened 2D array to store the predecessors for path reconstruction.
    predecessors[(i * nb_vertex) + j] stores the predecessor of vertex j in the shortest path from
    vertex i to vertex j */

    int *distances; /**< Flattened 2D array to store the distances for path reconstruction.
    distances[(i * nb_vertex) + j] stores the shortest distance from vertex i to vertex j */

    Path *faces; /**< Dynamic array of faces of the graph, where faces[i] is the (i+1)-th face of
    the graph */
    int face_basis; /**< Index of the face basis in minimals_basis, or -1 if the face basis has not been computed. */
    int face_basis_outer; /**< Index of the face basis in minimals_basis including the outer face,
    or -1 if the face basis has not been computed. */

    int nb_horton_cycles; /**< Number of cycles generated by Horton's algorithm,
    stored in horton_cycles */
    int capacity_horton_cycles; /**< Capacity of the allocated horton_cycles array */
    Path *horton_cycles; /**< Dynamic array of cycles generated by Horton's algorithm,
    where horton_cycles[i] is the (i+1)-th cycle generated */
    int basis_dimension; /**< Dimension of the cycle space of the graph (i.e., the number of cycles
    in any cycle basis). */
    Minimal_basis *minimals_basis; /**< Dynamic array of differents minimal cycle bases found, where
    minimals_basis[i] is the (i+1)-th minimal cycle basis founded */
    int nb_minimal_bases; /**< Number of differents minimal cycle bases found */
    int nb_faces; /**< Number of faces in the graph */
    int outer_face; /**< index of the outer face in Graph->faces */
} Graph;

Graph *create_graph();

void create_vertex(Graph *g, double x, double y);

void create_edge(Graph *g, int v1_id, int v2_id);

void save_graph(const Graph *g, const char *filename);

void delete_graph(Graph *g);

void multiple_horton(Graph *g, int *inv, int max_permutations);

void load_graph(Graph *g, const char *filename);

void horton(Graph *g);

void find_faces(Graph *g);

void split_edges(Graph *g, const int *edge_ids, int number_edge_to_split, int number_vertex_to_add);

void split_edge(Graph *g, int edge_id, int number_vertex_to_add);

void delete_vertex(Graph *g, int v_id);

void delete_edge(Graph *g, int e_id);

void move_vertex(Graph *g, int v, double new_x, double new_y);

void sort_neighbors_by_atan(const Graph *g, int v_id);

#endif //TER_GRAPHES_PLANAIRES_GRAPH_H