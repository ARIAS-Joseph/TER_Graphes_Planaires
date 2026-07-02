#ifndef DFS_GRAPH_H
#define DFS_GRAPH_H

#include <stdint.h>
#include "graph.h"

typedef enum {
    TREE_EDGE,
    BACK_EDGE
} DfsEdgeType;

typedef struct {
    int low1;
    int vertex;
    int *edges;
    int size;
} SingularSet;

typedef struct {
    int id;

    int from;
    int to;

    DfsEdgeType type;

    int parent_edge;

    int *successors;
    int nb_successors;

    int subtree_start;
    int subtree_end;

    int low1;
    int low2;

    int phi;
    int singular;
    int is_reference;

    int label;
} DfsEdge;

typedef struct {
    int id;

    int parent;
    int parent_edge;

    int dfs_num;

    int first_child;
    int nb_children;

    int *children;

    int *out_edges; /* Phi(v) : indices dans dfs->edges, triés par phi croissant */
    int nb_out_edges;

    int ancestor_start;
    int ancestor_end;

    int ref_edge;

} DfsVertex;

typedef struct {
    int *edges;
    int  size;
} SameGroup;

typedef struct {
    SameGroup *groups;
    int        nb_groups;
} SAME_t;

typedef struct {
    int *partner;   /* DIFF.partner[i] = indice du groupe partenaire, ou -1 */
    int  nb_groups;
} DIFF_t;

typedef struct {
    int vertices_count;
    int edge_count;
    int edge_capacity;

    int root;

    DfsVertex *vertices;
    DfsEdge *edges;

    int *E0;
    int E0_size;

    SingularSet *singular_sets;
    int nb_singular_sets;

    int *same_parent;   /* union-find : same_parent[e] = e pour les racines (E0 seulement) */
    int *same_partner;  /* same_partner[root] = root du groupe DIFF-partenaire, ou -1 */

    SAME_t SAME;          /* resultat final, lisible */
    DIFF_t DIFF;

    SAME_t SAME_prime;
    DIFF_t DIFF_prime;

    int *edge_to_group;   /* edge_to_group[e] = indice du groupe SAME de e,
                              ou -1 si e n'est pas dans E0 */

    /* FIX (bug #4) : indice SAME -> indice SAME' (ou -1 si le groupe a disparu
       de SAME', c.-a-d. qu'il ne contenait que des arcs singuliers non-
       representants). Necessaire pour que apply_partition_bits lise le bon
       espace d'indices. */
    int *same_to_prime;

    int *edge_indices;

} DfsGraph;

typedef struct Block {
    int *attachments;
    int size;
    int capacity;

    int representative;
} Block;

typedef struct {
    Block left;
    Block right;
} BlockPair;

typedef struct {
    BlockPair *pairs;

    int size;
    int capacity;
} AttList;

typedef struct {
    int **rotation_system;
    int *rotation_lengths;  /* rotation_lengths[v] = nb d'arcs dans rotation_system[v] */
    int n;                  /* nombre de sommets */
} PlanarEmbedding;

typedef struct {
    PlanarEmbedding *embeddings;
    int count;
} EmbeddingSet;

typedef struct {
    int *pair_g1;   /* racines des paires DIFF non-triviales (gauche), indices SAME' */
    int *pair_g2;   /* racines des paires DIFF non-triviales (droite), indices SAME' */
    int npairs;     /* d */
    int *free_g;    /* racines des groupes libres, indices SAME' */
    int nfree;      /* s */
} ReducedSystem;


typedef struct {
    int *left;   int nleft;   /* aretes du cote de x (cote inchange)   */
    int *right;  int nright;  /* aretes du cote oppose a x             */
} SingularVariant;

typedef struct {
    int *edges;  /* M'(v) : aretes sortantes de v, dans l'ordre final */
    int n;
} PartialMap;

typedef struct {
    int *order;  /* M(v) : ordre horaire complet (incl. arete entrante) */
    int n;
} VertexRotation;

typedef struct {
    DfsGraph *dfs;
    int *edge_set_id;      /* edge_set_id[e] = indice du singular_set de e, ou -1 */
    EmbeddingSet *result;
} EmbedCtx;

/* side[e] = 0 (LL) ou 1 (RR), valide uniquement pour e dans E0 */
typedef struct { int *side; } DynPartition;

typedef struct {
    DfsGraph *dfs;
    DynPartition *part;
    SingularVariant **chosen;  /* chosen[s] = variante choisie pour singular_sets[s] */
    int *edge_set_id;
    EmbeddingSet *result;
} VariantCtx;

int is_ancestor(const DfsGraph *dfs, int u, int v);
int is_descendant_edge(const DfsGraph *dfs, int e1, int e2);

DfsGraph *build_dfs_graph(Graph *graph);

void compute_low_values(const DfsGraph *dfs);

void build_phi_lists(DfsGraph * dfs);

void build_singular_sets(DfsGraph *dfs);

void compute_same_diff(DfsGraph *dfs);

void build_SAME_DIFF_prime(DfsGraph *dfs);

EmbeddingSet enumerate_embeddings(DfsGraph *dfs);

EmbeddingSet generate_all_embeddings(Graph *g);

void free_embedding_set(EmbeddingSet *es);

void free_dfs_graph(DfsGraph *dfs);

int is_biconnected(const Graph *g);

/* dfs_graph.h - ajouter les prototypes */
long long factorial(int n);
long long theoretical_embedding_count(const DfsGraph *dfs);
int      verify_embedding_count(const DfsGraph *dfs, int actual_count);

#endif