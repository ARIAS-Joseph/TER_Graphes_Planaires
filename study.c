/**
 * @file study.c
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <stdint.h>
#include <tgmath.h>

#ifdef _WIN32
  #include <direct.h>
  #define make_dir(path) _mkdir(path)
#else
  #include <sys/stat.h>
  #define make_dir(path) mkdir((path), 0755)
#endif

#include "graph.h"
#include "planar_graph_creator.h"
#include "dfs_graph.h"
#include "dif_faces.h"

/**
 * @brief Compute the total weight (sum of cycle lengths) of the first
 *        stored MCB. All MCBs of a graph share the same total weight.
 */
static int compute_mcb_weight(const Graph *g) {
    if (g->nb_minimal_bases == 0) return -1;
    int w = 0;
    for (int c = 0; c < g->minimals_basis[0].dimension; c++)
        w += g->minimals_basis[0].cycles[c].length;
    return w;
}

int popcount64(uint64_t x) {
    int c = 0; while (x) { x &= x - 1; c++; } return c;
}

/**
 * @brief Return the number of time internal faces can be a MCB.
 *
 * We try every face as the candidate outer face (excluded from the internal
 * set). For each candidate, we sum the edge counts of the remaining D faces.
 * If this sum equals @p mcb_weight, the internal faces are a minimal basis.
 *
 * @param faces      Face list for one Cai embedding (from trace_faces()).
 * @param D          g->basis_dimension  (= E − V + 1).
 * @param mcb_weight Sum of cycle lengths of any known MCB.
 */
static int embedding_faces_are_mcb(const FaceList *faces, const int D, const int mcb_weight) {

    if (faces->count != D + 1) {
        return -1;
    }

    int nb_valid_outer = 0;

    for (int skip = 0; skip < faces->count; skip++) {
        int total = 0;
        for (int f = 0; f < faces->count; f++) {
            if (f == skip) continue;
            total += popcount64(faces->edges_masks[f]);
        }
        if (total == mcb_weight) nb_valid_outer++;
    }

    return nb_valid_outer;
}

/**
 * @brief Test whether every edge appears in at most 2 cycles in @p mb.
 *
 * @param g   Graph (for edge count / deletion flags).
 * @param mb  Minimal basis to inspect.
 * @return    1 if no edge appears in 3 or more cycles, 0 otherwise.
 */
static int basis_edge_at_most_2_cycles(const Graph *g, const Minimal_basis *mb) {
    for (int e = 0; e < g->nb_edges; e++) {
        if (g->edges[e].deleted) continue;
        int count = 0;
        for (int c = 0; c < mb->dimension; c++)
            count += mb->cycles[c].edges_ids[e];
        if (count > 2) return 0;
    }
    return 1;
}

/**
 * @brief Return 1 iff every minimal basis of @p g satisfies the at-most-2
 *        cycles per edge condition.
 */
static int all_bases_edge_at_most_2(const Graph *g) {
    for (int b = 0; b < g->nb_minimal_bases; b++)
        if (!basis_edge_at_most_2_cycles(g, &g->minimals_basis[b]))
            return 0;
    return 1;
}

/* =========================================================================
 * MCB vs embedding faces comparison
 * ========================================================================= */

static uint64_t cycle_to_edges_mask(const Graph *g, const Path *cycle) {
    uint64_t mask = 0;
    for (int e = 0; e < g->nb_edges; e++) {
        if (g->edges[e].deleted) continue;
        if (cycle->edges_ids[e]) mask |= (1ULL << e);
    }
    return mask;
}

/**
 * @brief Test whether @p mb matches the internal faces of @p faces.
 *
 * Returns 1 iff the D MCB edge-masks equal the D non-outer face edge-masks
 * for some choice of outer face.
 */
static int mcb_matches_embedding(const Graph *g, const Minimal_basis *mb, const FaceList faces) {

    int D = mb->dimension;
    if (faces.count != D + 1) {
        printf("mcb_matches_embedding: face count %d != basis dimension + 1 (%d)\n",
               faces.count, D + 1);
        return -1;
    }

    uint64_t *mcb = malloc(D * sizeof(uint64_t));
    if (!mcb) { perror("mcb_matches_embedding"); exit(EXIT_FAILURE); }
    for (int c = 0; c < D; c++)
        mcb[c] = cycle_to_edges_mask(g, &mb->cycles[c]);

    int result = 0;
    for (int skip = 0; skip < faces.count && !result; skip++) {
        int all_found = 1;
        for (int c = 0; c < D && all_found; c++) {
            int found = 0;
            for (int f = 0; f < faces.count && !found; f++) {
                if (f == skip) continue;
                if (faces.edges_masks[f] == mcb[c]) found = 1;
            }
            if (!found) all_found = 0;
        }
        if (all_found) result = 1;
    }

    free(mcb);
    return result;
}

/* =========================================================================
 * File-system helpers
 * ========================================================================= */

static void ensure_dir(const char *path) {
    if (make_dir(path) == -1 && errno != EEXIST) {
        perror(path); exit(EXIT_FAILURE);
    }
}

static void save_to_folder(const char *folder, const Graph *g, const char *filename) {
    ensure_dir(folder);
    char path[1024];
    printf("Saving graph to %s/graph_%s.txt\n", folder, filename);
    snprintf(path, sizeof(path), "%s/graph_%s.txt", folder, filename);
    save_graph(g, path);
}

/**
 * @brief Compte le nombre de composantes connexes du graphe H des chaînes entre
 *        u et v, en utilisant la réduction polynomiale suivante.
 *
 * Soit H le graphe dont les sommets sont les chaînes simples de u à v dans G
 * (y compris l'arête directe (u,v)), et dont les arêtes relient deux chaînes
 * qui partagent une arête de G ou un sommet intérieur (≠ u, v).
 *
 * Réduction (démontrée dans has_four_chain_components) :
 *
 *   nb_composantes(H)  =  1  +  |{ C composante de G−{u,v}  |
 *                                    C a un voisin de u  ET  un voisin de v }|
 *
 * Complexité : O(V + E).
 *
 * @param g  Graphe G.
 * @param u  Premier sommet de l'arête (u,v).
 * @param v  Second  sommet de l'arête (u,v).
 * @return   Nombre de composantes connexes de H.
 */
static int count_chain_components(const Graph *g, const int u, const int v) {
    const int n = g->nb_vertices;

    /* ── Étape 1 : BFS pour étiqueter les composantes de G − {u, v} ──────
     *
     * On retire u et v du graphe en les ignorant pendant le parcours.
     * comp[w] = indice de la composante de w dans G−{u,v}, ou -1 si w n'est
     * pas encore visité (ou si w = u ou w = v).                            */
    int *comp  = malloc(n * sizeof(int));
    int *queue = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) comp[i] = -1;

    int nb_comp = 0;

    for (int s = 0; s < n; s++) {
        /* Ignorer u, v et les sommets supprimés. */
        if (s == u || s == v || g->vertices[s].deleted) continue;
        if (comp[s] != -1) continue;   /* déjà visité */

        /* BFS depuis s dans G − {u, v}. */
        int head = 0, tail = 0;
        queue[tail++] = s;
        comp[s] = nb_comp;

        while (head < tail) {
            const int cur = queue[head++];
            for (int i = 0; i < g->neighbors[cur].count; i++) {
                const int w = g->neighbors[cur].neighbors[i];
                /* Ne pas traverser u ni v ni les sommets supprimés. */
                if (w == u || w == v || g->vertices[w].deleted) continue;
                if (comp[w] == -1) {
                    comp[w] = nb_comp;
                    queue[tail++] = w;
                }
            }
        }
        nb_comp++;   /* nouvelle composante découverte */
    }

    free(queue);

    /* Cas dégénéré : aucun sommet hors de {u,v} → seule l'arête directe. */
    if (nb_comp == 0) {
        free(comp);
        return 1;
    }

    /* ── Étape 2 : déterminer quelles composantes touchent u ET v ─────────
     *
     * Une composante C contribue à H si et seulement si elle possède au
     * moins un voisin de u ET au moins un voisin de v dans G.
     * (Les voisins de u ou v qui sont u ou v eux-mêmes sont exclus.)      */
    int *touches_u = calloc(nb_comp, sizeof(int));
    int *touches_v = calloc(nb_comp, sizeof(int));

    for (int i = 0; i < g->neighbors[u].count; i++) {
        const int w = g->neighbors[u].neighbors[i];
        if (w == v || g->vertices[w].deleted) continue;
        if (comp[w] != -1) touches_u[comp[w]] = 1;
    }

    for (int i = 0; i < g->neighbors[v].count; i++) {
        const int w = g->neighbors[v].neighbors[i];
        if (w == u || g->vertices[w].deleted) continue;
        if (comp[w] != -1) touches_v[comp[w]] = 1;
    }

    /* ── Étape 3 : compter ────────────────────────────────────────────────
     *
     * On part de 1 (l'arête directe (u,v) est toujours une composante
     * isolée de H : elle n'a aucun sommet intérieur et son arête (u,v)
     * n'est utilisée par aucune autre chaîne).
     * Chaque composante de G−{u,v} touchant à la fois u et v ajoute 1.   */
    int h_components = 1;
    for (int c = 0; c < nb_comp; c++)
        if (touches_u[c] && touches_v[c]) h_components++;

    free(comp);
    free(touches_u);
    free(touches_v);

    return h_components;
}

/**
 * @brief Teste si le graphe G admet une arête (u,v) telle que le graphe des
 *        chaînes H ait au moins 4 composantes connexes.
 *
 * Condition sur H :
 *   Les sommets de H sont les chaînes simples de u à v dans G (arête directe
 *   incluse). Deux chaînes y sont adjacentes si elles partagent une arête de G
 *   ou un sommet intérieur (≠ u, v).
 *
 * Hypothèse (à vérifier expérimentalement) :
 *   Si une telle arête existe dans G, alors pour tout plongement planaire de G
 *   les faces internes ne forment jamais une base de cycles minimale.
 *
 * Réduction polynomiale utilisée (voir count_chain_components) :
 *   nb_composantes(H) = 1 + |{ composantes de G−{u,v} touchant u et v }|
 *
 * Complexité totale : O(E · (V + E)).
 *
 * @param g  Graphe G (2-connexe planaire).
 * @return   Indice de la première arête (u,v) témoin, ou -1 si aucune.
 */
int has_four_chain_components(const Graph *g)
{
    for (int e = 0; e < g->nb_edges; e++) {
        if (g->edges[e].deleted) continue;

        const int u = g->edges[e].u;
        const int v = g->edges[e].v;

        const int h_comp = count_chain_components(g, u, v);

        if (h_comp >= 4) {
            /* Arête témoin trouvée : on peut afficher les détails ici. */
            return e;
        }
    }

    return -1;   /* aucune arête de G ne remplit la condition */
}

static int count_components_after_path(const Graph *g, const int *path, const int path_len)
{
    const int n = g->nb_vertices;

    int *removed = calloc(n, sizeof(int));

    /* supprimer uniquement les sommets internes */
    for (int i = 1; i < path_len - 1; i++)
        removed[path[i]] = 1;

    int *comp = malloc(n * sizeof(int));
    int *queue = malloc(n * sizeof(int));

    for (int i = 0; i < n; i++)
        comp[i] = -1;

    int nb_comp = 0;

    for (int s = 0; s < n; s++) {

        if (removed[s] || g->vertices[s].deleted)
            continue;

        if (comp[s] != -1)
            continue;

        int head = 0;
        int tail = 0;

        queue[tail++] = s;
        comp[s] = nb_comp;

        while (head < tail) {

            int cur = queue[head++];

            for (int i = 0; i < g->neighbors[cur].count; i++) {

                int w = g->neighbors[cur].neighbors[i];

                if (removed[w] || g->vertices[w].deleted)
                    continue;

                if (comp[w] == -1) {
                    comp[w] = nb_comp;
                    queue[tail++] = w;
                }
            }
        }

        nb_comp++;
    }

    free(queue);
    free(comp);
    free(removed);

    return nb_comp;
}

static int dfs_paths(const Graph *g, const int current, const int target, int *visited, int *path,
    const int path_len) {

    if (current == target) {

        if (count_components_after_path(g, path, path_len) == 3) {
            // printf("Found path: %d -> %d (len %d). Path:\n", current, target, path_len);
            // for (int i = 0; i < path_len; i++) {
            //     printf("%d ", path[i]);
            // }
            // printf("\n");
            return 1;
        }

        return 0;
    }

    for (int i = 0; i < g->neighbors[current].count; i++) {

        int nxt = g->neighbors[current].neighbors[i];

        if (g->vertices[nxt].deleted) continue;

        if (visited[nxt]) continue;

        visited[nxt] = 1;
        path[path_len] = nxt;

        if (dfs_paths(g, nxt, target, visited, path, path_len + 1)) {
            return 1;
        }

        visited[nxt] = 0;
    }

    return 0;
}

int exists_three_component_separator(const Graph *g) {
    const int n = g->nb_vertices;

    int *visited = calloc(n, sizeof(int));
    int *path = malloc(n * sizeof(int));

    for (int u = 0; u < n; u++) {

        if (g->vertices[u].deleted) continue;

        for (int v = u + 1; v < n; v++) {

            if (g->vertices[v].deleted) continue;

            for (int i = 0; i < n; i++) visited[i] = 0;

            visited[u] = 1;
            path[0] = u;

            if (dfs_paths(g, u, v, visited, path, 1)) {
                free(path);
                free(visited);
                return 1;
            }
        }
    }

    free(path);
    free(visited);

    return 0;
}

/**
 * @brief Teste si le graphe G−{u,v} reste connexe (BFS en ignorant u et v).
 *
 * Sert de brique pour is_3_connected() : G est 3-connexe ssi pour toute
 * paire {u,v}, G−{u,v} est connexe.
 *
 * @param g  Graphe planaire 2-connexe.
 * @param u  Premier sommet supprimé.
 * @param v  Second  sommet supprimé.
 * @return   1 si G−{u,v} est connexe, 0 sinon.
 */
static int is_connected_without(const Graph *g, int u, int v) {

    const int n = g->nb_vertices;
    int *visited = calloc(n, sizeof(int));
    int *queue   = malloc(n * sizeof(int));

    /* Premier sommet actif (ni u, ni v, ni supprimé). */
    int start = -1;
    for (int i = 0; i < n; i++) {
        if (i != u && i != v && !g->vertices[i].deleted) { start = i; break; }
    }

    int connected = 1;

    if (start == -1) goto cleanup;  /* 0 ou 1 sommet restant : trivial */

    /* BFS depuis start dans G − {u, v}. */
    int head = 0, tail = 0;
    queue[tail++] = start;
    visited[start] = 1;

    while (head < tail) {
        const int cur = queue[head++];
        for (int i = 0; i < g->neighbors[cur].count; i++) {
            const int w = g->neighbors[cur].neighbors[i];
            if (w == u || w == v || g->vertices[w].deleted || visited[w]) continue;
            visited[w] = 1;
            queue[tail++] = w;
        }
    }

    /* Tous les sommets actifs doivent avoir été atteints. */
    for (int i = 0; i < n; i++) {
        if (i != u && i != v && !g->vertices[i].deleted && !visited[i]) {
            connected = 0;
            break;
        }
    }

cleanup:
    free(visited);
    free(queue);
    return connected;
}

/**
 * @brief Teste la 3-connexité de G.
 *
 * G est 3-connexe ssi il est 2-connexe (vérifié par is_biconnected) ET
 * pour toute paire de sommets {u,v}, G−{u,v} est connexe.
 *
 * Complexité : O(V² × (V+E)).  Suffisant pour V ≤ 20.
 *
 * @param g  Graphe planaire 2-connexe.
 * @return   1 si G est 3-connexe, 0 sinon.
 */
static int is_3_connected(const Graph *g) {

    if (!is_biconnected(g)) return 0;

    const int n = g->nb_vertices;
    for (int u = 0; u < n; u++) {
        if (g->vertices[u].deleted) continue;
        for (int v = u + 1; v < n; v++) {
            if (g->vertices[v].deleted) continue;
            if (!is_connected_without(g, u, v)) return 0;
        }
    }
    return 1;
}

/**
 * @brief Compte le nombre de bits à 1 d'un masque 64 bits.
 */
static int popcount_local(uint64_t x) {
    int c = 0; while (x) { x &= x - 1; c++; } return c;
}

static int is_necklace(const Graph *g,
                       uint64_t outer_edges)
{
    const int n = g->nb_vertices;

    int *deg = calloc(n, sizeof(int));
    int *alive = calloc(n, sizeof(int));
    int *queue = malloc(n * sizeof(int));

    /* degrés dans l'arbre T */
    for (int e = 0; e < g->nb_edges; e++) {

        if (g->edges[e].deleted)
            continue;

        if (outer_edges & (1ULL << e))
            continue;               /* arête du cycle extérieur */

        int u = g->edges[e].u;
        int v = g->edges[e].v;

        deg[u]++;
        deg[v]++;
    }

    for (int v = 0; v < n; v++)
        alive[v] = !g->vertices[v].deleted;

    /* effeuillage */

    int head = 0, tail = 0;

    for (int v = 0; v < n; v++)
        if (alive[v] && deg[v] == 1)
            queue[tail++] = v;

    while (head < tail) {

        int v = queue[head++];

        if (!alive[v] || deg[v] != 1)
            continue;

        alive[v] = 0;

        for (int e = 0; e < g->nb_edges; e++) {

            if (g->edges[e].deleted)
                continue;

            if (outer_edges & (1ULL << e))
                continue;

            int u = -1;

            if (g->edges[e].u == v)
                u = g->edges[e].v;
            else if (g->edges[e].v == v)
                u = g->edges[e].u;
            else
                continue;

            if (!alive[u])
                continue;

            deg[u]--;

            if (deg[u] == 1)
                queue[tail++] = u;
        }
    }

    /* les sommets restants doivent former un chemin */

    int remain = 0;
    int deg1 = 0;
    int deg2 = 0;

    for (int v = 0; v < n; v++) {

        if (!alive[v])
            continue;

        remain++;

        int d = 0;

        for (int e = 0; e < g->nb_edges; e++) {

            if (g->edges[e].deleted)
                continue;

            if (outer_edges & (1ULL << e))
                continue;

            int u = -1;

            if (g->edges[e].u == v)
                u = g->edges[e].v;
            else if (g->edges[e].v == v)
                u = g->edges[e].u;
            else
                continue;

            if (alive[u])
                d++;
        }

        if (d == 1)
            deg1++;
        else if (d == 2)
            deg2++;
        else {
            free(deg);
            free(alive);
            free(queue);
            return 0;
        }
    }

    free(deg);
    free(alive);
    free(queue);

    if (remain == 2)
        return 1;      /* K4 */

    return (deg1 == 2 && deg2 == remain - 2);
}

/**
 * @brief Teste si le graphe G est un graphe de Halin.
 *
 * Utilise la caractérisation de Wikipedia (méthode 2) :
 * G est un graphe de Halin si et seulement s'il est
 *   (a) planaire,
 *   (b) 3-connexe,
 *   (c) il existe une face dont le nombre de sommets est égal au rang
 *       cyclomatique  D = m − n + 1.
 *
 * Justification de (c) : dans un graphe de Halin avec k feuilles, l'arbre
 * couvrant a n−1 arêtes, le cycle extérieur en ajoute k, donc m = n−1+k et
 * D = k.  La face extérieure (cycle des feuilles) a exactement k sommets. ✓
 *
 * La planéité est garantie par le pipeline (embs->count > 0).
 * Pour un graphe 3-connexe, Whitney garantit l'unicité du plongement, donc
 * il suffit de tester le premier ; on parcourt tous par robustesse.
 *
 * Complexité : O(V² × (V+E))  dominée par is_3_connected.
 *
 * @param g     Graphe planaire 2-connexe avec g->nb_edges et g->nb_vertices
 *              à jour (pas d'arêtes/sommets marqués deleted).
 * @param dfs   Graphe DFS de g (déjà construit).
 * @param embs  Ensemble de tous les plongements planaires de g.
 * @param faces
 * @return      1 si G est un graphe de Halin, 0 sinon.
 */
int is_halin_graph(const Graph *g, const DfsGraph *dfs, const EmbeddingSet *embs, const FaceList *faces) {

    /* (a) Planéité */
    if (embs->count == 0) return 0;

    /* (b) 3-connexité */
    if (!is_3_connected(g)) return 0;

    /* (c) Existence d'une face de taille D = m − n + 1 */
    const int D = g->nb_edges - g->nb_vertices + 1;

    for (int ei = 0; ei < embs->count; ei++) {

        for (int f = 0; f < faces[ei].count; f++) {
            if (popcount_local(faces[ei].vertices_masks[f]) == D) {

                if (is_necklace(g, faces[ei].edges_masks[f])) {
                    return 0;
                }

                return 1;
            }
        }
    }

    return 0;
}

int is_wheel_graph(const Graph *g) {
/**je t'aime */
    int v_center = -1;
    for (int v=0; v < g->nb_vertices; v++) {
        if (g->neighbors[v].count == g->nb_vertices - 1) {
            v_center = v;
            break;
        }
    }

    if (v_center == -1) {return 0;}

    int all_v_deg2 = 1;
    for (int v=0; v < g->nb_vertices; v++) {
        if (v == v_center) {continue;}
        int d = g->neighbors[v].count;
        d -= g->edge_indices[v * g->nb_vertices + v_center] == - 1 ? 0 : 1;
        if (d != 2) {
            all_v_deg2 = 0;
            break;
        }
    }

    Graph* temp = create_graph();
    for (int v=0; v < g->nb_vertices; v++) {
        if (v == v_center) {continue;}
        create_vertex(temp, 0, 0);
    }

    for (int e=0; e < g->nb_edges; e++) {
        if (g->edges[e].u == v_center || g->edges[e].v == v_center) {continue;}
        create_edge(temp, g->edges[e].u, g->edges[e].v);
    }

    const int two_connex = is_biconnected(temp);
    delete_graph(temp);

    if (all_v_deg2 == 1 && two_connex) { return 1; }
    return 0;
}

/* =========================================================================
 * Core study loop
 * ========================================================================= */

void study_graph(const int nb_vertex, const int nb_tests, int *id, FILE *f) {

    const int min_edges = nb_vertex-1;
    const int max_edges = 3 * nb_vertex - 6 > 64 ? 64 : 3 * nb_vertex - 6;
    if (min_edges > max_edges) return;

    for (int trial = 0; trial < nb_tests; trial++) {

        const int target_edges = min_edges + rand() % (max_edges - min_edges);

        Graph *g = create_planar_graph(nb_vertex, target_edges);
        reduce_graph(g);
        if (!g) {
            trial--;
            continue;
        }

        if (g->nb_edges < min_edges || !is_biconnected(g)) { delete_graph(g); trial--; continue; }

        int *inv_e = calloc(g->nb_edges,    sizeof(int));
        int *inv_v = calloc(g->nb_vertices, sizeof(int));
        if (!inv_e || !inv_v) {
            free(inv_e); free(inv_v); delete_graph(g); trial--; continue;
        }

        multiple_horton(g, inv_e, inv_v, 10000);
        free(inv_e); free(inv_v);

        const int expected_dim = g->nb_edges - g->nb_vertices + 1;
        if (g->minimals_basis[0].dimension != expected_dim) {
            delete_graph(g); trial--; continue;
        }

        printf("n=%d m=%d id=%d", g->nb_vertices, g->nb_edges, *id);

        const int edge_max_2_cycle_all_bcm = all_bases_edge_at_most_2(g);

        DfsGraph *dfs = build_dfs_graph(g);
        compute_low_values(dfs);
        build_phi_lists(dfs);
        build_singular_sets(dfs);
        compute_same_diff(dfs);
        build_SAME_DIFF_prime(dfs);

        const EmbeddingSet embs = enumerate_embeddings(dfs);

        printf("\ntrial=%d/%d embeddings=%d\n", trial + 1, nb_tests, embs.count);

        const int mcb_weight = compute_mcb_weight(g);

        FaceList *cached_faces = malloc(embs.count * sizeof(FaceList));
        for (int ei = 0; ei < embs.count; ei++)
            cached_faces[ei] = trace_faces(dfs, &embs.embeddings[ei]);

        FaceList all_faces;
        face_list_init(&all_faces);
        int at_least_one_map_is_MCB = 0;
        int is_outerplanar = 0;
        const int is_wheel = is_wheel_graph(g);
        int is_halin = is_wheel ? 1 : is_halin_graph(g, dfs, &embs, cached_faces);

        const uint64_t all_verts = (1ULL << g->nb_vertices) - 1;

        int error = 0;
        int nb_map_is_bcm = 0;

        for (int ei = 0; ei < embs.count; ei++) {

            FaceList faces = cached_faces[ei];
            const int expected_faces = dfs->edge_count - dfs->vertices_count + 2;

            if (faces.count != expected_faces) {
                fml_free(&faces);
                error = 1;
                break;
            }

            const int face_are_mcb = embedding_faces_are_mcb(&faces, g->basis_dimension, mcb_weight);
            if (face_are_mcb == -1) { error = 1; break; }

            nb_map_is_bcm += face_are_mcb;
            if (face_are_mcb > 0) at_least_one_map_is_MCB = 1;

            if (!is_outerplanar) {
                for (int i = 0; i < faces.count; i++) {
                    if (faces.vertices_masks[i] == all_verts) {
                        is_outerplanar = 1;
                        at_least_one_map_is_MCB = 1;
                        break;}
                }
            }

            for (int i = 0; i < faces.count; i++)
                face_list_push(&all_faces, faces.vertices_masks[i], faces.edges_masks[i]);
        }

        if (error) {
            free_embedding_set(&embs);
            free_dfs_graph(dfs);
            delete_graph(g);
            trial--;
            continue;
        }

        fml_free(&all_faces);

        int nb_bases_with_map = 0;
        for (int mb = 0; mb < g->nb_minimal_bases; mb++) {
            for (int ei = 0; ei < embs.count; ei++) {
                FaceList faces = cached_faces[ei];
                if (mcb_matches_embedding(g, &g->minimals_basis[mb], faces)) {
                    nb_bases_with_map++;
                    break;
                }
            }
        }
        const int all_bcm_are_maps = nb_bases_with_map == g->nb_minimal_bases;
        int has_max_edges = is_outerplanar ? g->nb_edges == 2 * g->nb_vertices - 3 : g->nb_edges == 3*g->nb_vertices - 6;
        float ratio = 100.0f * nb_map_is_bcm / (float)(embs.count * (g->basis_dimension + 1));


        int has_sparse_not_map = 0;

        for (int mb = 0; mb < g->nb_minimal_bases; mb++) {

            if (!basis_edge_at_most_2_cycles(g, &g->minimals_basis[mb])) continue;

            int found_embedding = 0;

            for (int ei = 0; ei < embs.count; ei++) {
                const FaceList faces = cached_faces[ei];
                if (mcb_matches_embedding(g, &g->minimals_basis[mb], faces)) {
                    found_embedding = 1;
                    break;
                }
            }

            if (!found_embedding) {
                has_sparse_not_map = 1;
                break;
            }
        }

        char filename[1024];
        snprintf(filename, sizeof(filename), "%d_v(%d)_e(%d)_outer(%d)_halin(%d)_wheel(%d)"
                                             "_wheelSize(%d)_maxEdges(%d)_sparseMCBNotMap(%d)_ratio(%d)_nbMCB(%d)",
                                             *id, dfs->vertices_count,
                                             dfs->edge_count, is_outerplanar, is_halin, is_wheel, is_wheel ? g->nb_vertices - 1 : -1, has_max_edges, has_sparse_not_map,
                                             (int)(100.0f * nb_map_is_bcm / (float)(embs.count * (g->basis_dimension + 1))),
                                             g->nb_minimal_bases);

        if (!at_least_one_map_is_MCB) {
            save_to_folder("no_map_is_BCM", g, filename);
            const int witness_edge = has_four_chain_components(g);
            if (witness_edge == -1) save_to_folder("no_map_is_BCM_no_4chains", g, filename);
            if (edge_max_2_cycle_all_bcm) save_to_folder("no_map_BCM_but_sparse_basis", g, filename);
        }

        if (edge_max_2_cycle_all_bcm) {
            save_to_folder("All_BCM_are_sparse", g, filename);
            if (!all_bcm_are_maps) save_to_folder("All_BCM_are_sparse_but_not_all_BCM_are_map", g, filename);
            else save_to_folder("All_BCM_are_sparse_and_all_BCM_are_map", g, filename);
        }

        if (all_bcm_are_maps) {
            save_to_folder("All_BCM_are_map", g, filename);
            if (is_halin_graph(g, dfs, &embs, cached_faces))
                save_to_folder("All_BCM_are_map_not_halin", g, filename);
        }

        if (nb_map_is_bcm == embs.count * (g->basis_dimension + 1)) {
            save_to_folder("All_maps_are_BCM", g, filename);
            if (!all_bcm_are_maps) save_to_folder("All_maps_are_BCM_but_not_all_BCM_are_map", g, filename);
            else save_to_folder("All_maps_are_BCM_and_all_BCM_are_map", g, filename);
        }

        if (nb_map_is_bcm == 1) save_to_folder("1_map_is_BCM", g, filename);

        if (nb_map_is_bcm > 1) {
            save_to_folder("several_maps_are_BCM", g, filename);
            if (!all_bcm_are_maps) save_to_folder("several_maps_are_BCM_but_not_all_BCM_are_map", g, filename);
            else save_to_folder("several_maps_are_BCM_and_all_BCM_are_map", g, filename);
        }

        fprintf(f,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%d\n",
        *id,
        dfs->vertices_count,
        dfs->edge_count,
        is_outerplanar,
        is_halin,
        is_wheel,
        is_wheel ? g->nb_vertices - 1 : -1,
        has_max_edges,
        has_sparse_not_map,
        nb_map_is_bcm == embs.count * (g->basis_dimension + 1),
        all_bcm_are_maps,
        nb_map_is_bcm == 0,
        nb_map_is_bcm == 1,
        ratio,
        g->nb_minimal_bases);

        fflush(f);

        (*id)++;

        free_embedding_set(&embs);
        free_dfs_graph(dfs);
        delete_graph(g);
        for (int ei = 0; ei < embs.count; ei++)
            fml_free(&cached_faces[ei]);
        free(cached_faces);
        g = NULL;
        }
}

int main(void) {
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);
    srand((unsigned int)time(NULL));

    FILE *f = fopen("graphs.csv", "w");

    fprintf(f, "id,|V|,|E|,is_outerplanar,is_halin,is_wheel,size_wheel,has_max_edges,has_sparse_MCB_not_map,All_maps_are_MCB,All_BCM_are_map,None_map_MCB,Only_one_map_MCB, nb_map_is_MCB / nb_map, nb_MCB\n");

    int id = 0;
    for (int n = 20; n <= 25; n++) {
        printf( "n = %d ...\n", n);
        study_graph(n, 1000, &id, f);
    }

    fclose(f);

    printf( "Done. %d graphs processed.\n", id);
    return 0;
}