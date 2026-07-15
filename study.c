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

/* =========================================================================
 * Maximum number of embeddings processed per graph.
 * Graphs with more embeddings are skipped: their enumeration would stall
 * the program for an impractical amount of time.
 * ========================================================================= */
#define MAX_EMBEDDINGS 5000

/* =========================================================================
 * Embedding count from Cai's formula  2^(d+s) * prod h(x)
 *
 * Computed BEFORE enumerate_embeddings() so we can skip expensive graphs.
 *
 * Needs only the exported fields of DfsGraph; no internal static function
 * is required. We inline a minimal union-find path-follower (no compression
 * needed here — correctness is all we need, not speed).
 * ========================================================================= */

/** Non-compressing union-find root (same_parent is a public DfsGraph field). */
static int uf_root(const int *parent, int x) {
    while (parent[x] != x) x = parent[x];
    return x;
}

/* =========================================================================
 * Check whether the internal faces of a Cai embedding form a minimal
 * cycle basis, WITHOUT relying on Horton having stored the right MCB.
 *
 * Key property: the D internal faces of a 2-connected planar graph always
 * span the cycle space (they are D independent cycles). They form a
 * MINIMAL cycle basis iff their total edge count equals the total weight
 * of any known MCB (all MCBs share the same total weight).
 *
 * This replaces mcb_matches_embedding(), which incorrectly required the
 * face-basis MCB to already be present in g->minimals_basis. With only
 * 1000 Horton permutations, MCBs matching non-geometric embeddings are
 * routinely missed.
 * ========================================================================= */

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

    printf("faces ptr=%p count=%d\n", (void*)faces, faces->count);

    printf("MCB input face count = %d\n", faces->count);

    for(int i=0;i<faces->count;i++)
    {
        printf("face %d length=%d\n",
               i,
               popcount64(faces->edges_masks[i]));
    }

    if (faces->count != D + 1) {
        printf("embedding_faces_are_mcb: face count %d != basis dimension + 1 (%d)\n",
               faces->count, D + 1);
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
 * @brief Compute the number of distinct planar embeddings of the graph
 *        encoded in @p dfs using Cai's counting formula.
 *
 * Formula (Theorem 1, Cai 1992):
 * @verbatim
 *   N = 2^(d+s)  *  product_{x in RS} h(x)
 *
 *   d   = number of non-trivial pairs in DIFF'
 *   s   = number of free (unpartnered) groups in SAME'
 *   RS  = set of representative singular edges
 *   g(x)= size of the singular set containing x
 *   h(x)= g(x)!          if x is "bound"  (singular set ⊆ one SAME group)
 *        = (g(x)+1)! / 2  if x is "free"
 * @endverbatim
 *
 * Must be called after build_SAME_DIFF_prime().
 *
 * @param dfs DFS graph with SAME'/DIFF' already built.
 * @return    Number of distinct planar embeddings, or INT64_MAX if the
 *            intermediate product overflows a 64-bit signed integer.
 */
static int64_t compute_nb_embeddings(const DfsGraph *dfs) {

    /* --- Step 1: count d (non-trivial DIFF' pairs) and s (free groups) --- */
    int ng = dfs->SAME_prime.nb_groups;
    int d = 0, s = 0;
    int *done = calloc(ng, sizeof(int));
    if (!done) { perror("compute_nb_embeddings: calloc"); exit(EXIT_FAILURE); }

    for (int i = 0; i < ng; i++) {
        if (done[i]) continue;
        int p = dfs->DIFF_prime.partner[i];
        if (p == -1) {
            s++;
            done[i] = 1;
        } else {
            d++;
            done[i] = done[p] = 1;
        }
    }
    free(done);

    /* --- Step 2: compute prod h(x) over all singular sets (= RS) --------- */
    int64_t h_prod = 1;

    for (int si = 0; si < dfs->nb_singular_sets; si++) {
        const SingularSet *set = &dfs->singular_sets[si];
        int g = set->size;   /* g(x) */

        /* A singular set is "bound" iff all its edges share the same
         * SAME-group root in the union-find.                          */
        int bound = 1;
        int r0 = uf_root(dfs->same_parent, set->edges[0]);
        for (int j = 1; j < g && bound; j++)
            if (uf_root(dfs->same_parent, set->edges[j]) != r0)
                bound = 0;

        /* h(x) = g! (bound) or (g+1)!/2 (free) */
        int64_t h = 1;
        int top = bound ? g : (g + 1);
        for (int j = 2; j <= top; j++) {
            /* Overflow guard: if h would exceed INT64_MAX, cap it */
            if (h > INT64_MAX / j) { h_prod = INT64_MAX; goto done_prod; }
            h *= j;
        }
        if (!bound) h /= 2;

        if (h_prod > INT64_MAX / h) { h_prod = INT64_MAX; goto done_prod; }
        h_prod *= h;
    }

done_prod:;
    /* --- Step 3: 2^(d+s) * h_prod --------------------------------------- */
    int shift = d + s;
    if (shift >= 62 || (INT64_MAX >> shift) < h_prod) return INT64_MAX;
    return ((int64_t)1 << shift) * h_prod;
}

/* =========================================================================
 * MCB edge-multiplicity filter
 *
 * Returns 1 iff in EVERY minimal basis found, every edge appears in AT MOST
 * 2 cycles of that basis. This is a property of the MCB alone (independent
 * of the embeddings).
 * ========================================================================= */

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

/* Regarde si une arete apparait tantot 2 fois tantot3 fois */
static int edge_two_and_three_appearence(const Graph *g, const Minimal_basis *mbs) {
    int *occurences = calloc(g->nb_edges, sizeof(int));
    for (int mb =0; mb < g->nb_minimal_bases; mb++) {
        for (int e = 0; e < g->nb_edges; e++) {
            if (g->edges[e].deleted) continue;
            int count = 0;
            for (int c = 0; c < g->minimals_basis[mb].dimension; c++)
                count += mbs[mb].cycles[c].edges_ids[e];
            if (occurences[e] == 0) occurences[e] = count;
            if ((occurences[e] <= 2 && count > 2) || (occurences[e] > 2 && count <= 2)) {
                free(occurences);
                return 1;
            }
        }
    }
    free(occurences);
    return 0;
}

/* Regarde si tantot toutes les aretes sont max 2 fois et parfois non */
static int diff_edge_occur(const Graph *g, const Minimal_basis *mbs) {
    int all_leq_2 = 1;
    int bcm_all_leq_2 = -1;
    for (int mb =0; mb < g->nb_minimal_bases; mb++) {
        for (int e = 0; e < g->nb_edges; e++) {
            if (g->edges[e].deleted) continue;
            int count = 0;
            for (int c = 0; c < g->minimals_basis[mb].dimension; c++)
                count += mbs[mb].cycles[c].edges_ids[e];
            if (count > 2) all_leq_2 = 0;
        }
        if (bcm_all_leq_2 == -1) {
            bcm_all_leq_2 = all_leq_2;
        } else if (bcm_all_leq_2 != all_leq_2) {
            return 1;
        }

    }
    return 0;
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
 * Never-co-facial search
 * ========================================================================= */

static int count_minimal_never_cofacial(int n, const FaceList *all_faces) {
    if (n > 64) return 0;
    SearchCtx ctx = { all_faces, NULL, 0, 0 };
    for (int k = 2; k <= n; k++)
        search_rec(n, 0, k, 0, &ctx);
    int count = ctx.nb_minimal;
    free(ctx.minimal);
    return count;
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
            printf("Found path: %d -> %d (len %d). Path:\n", current, target, path_len);
            for (int i = 0; i < path_len; i++) {
                printf("%d ", path[i]);
            }
            printf("\n");
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
 * @brief Teste si chaque base de cycles minimale trouvée correspond aux faces
 *        internes d'au moins un plongement planaire.
 *
 * Pour chaque base B dans g->minimals_basis[], on cherche un plongement de
 * `embs` tel que l'ensemble des D faces internes de ce plongement soit
 * exactement l'ensemble des D cycles de B (comparaison par masque d'arêtes).
 *
 * On utilise deux critères complémentaires :
 *  1. Critère de poids  (rapide, O(F)) :
 *     Les faces internes forment une BCM ⟺ leur poids total = mcb_weight.
 *     Suffit à réfuter rapidement les plongements non-BCM.
 *
 *  2. Comparaison directe des masques (O(D²)) :
 *     Si le poids correspond, on vérifie que les D masques de faces sont
 *     exactement les D masques des cycles de B (en ignorant la face externe,
 *     essayée à tour de rôle).
 *
 * @param g          Graphe avec ses BCM calculées par multiple_horton().
 * @param dfs        Graphe DFS associé à g (déjà construit).
 * @param embs       Ensemble de tous les plongements planaires de g.
 * @param mcb_weight Somme des longueurs des cycles d'une BCM
 *                   (invariant : identique pour toutes les BCM).
 * @return           1 si TOUTES les bases correspondent à un plongement,
 *                   0 si au moins une base ne correspond à aucun plongement.
 */
static int all_bases_are_maps(const Graph *g, const DfsGraph *dfs, const EmbeddingSet *embs,
    const int mcb_weight) {

    const int D = g->basis_dimension;

    /* Précalcul des masques d'arêtes de chaque plongement, pour chaque choix
     * de face externe possible. On stocke (nb_embeddings × (D+1)) masques. */

    /* Pour éviter de recalculer trace_faces plusieurs fois, on le fait une
     * seule fois par plongement et on stocke la FaceList. */
    FaceList *all_face_lists = malloc(embs->count * sizeof(FaceList));
    if (!all_face_lists) { perror("all_bases_are_maps"); exit(EXIT_FAILURE); }

    for (int ei = 0; ei < embs->count; ei++)

        all_face_lists[ei] = trace_faces(dfs, &embs->embeddings[ei]);

    int all_match = 1;

    for (int b = 0; b < g->nb_minimal_bases && all_match; b++) {

        const Minimal_basis *mb = &g->minimals_basis[b];

        /* Calcul des masques d'arêtes des D cycles de la base b. */
        uint64_t *mcb_masks = malloc(D * sizeof(uint64_t));
        if (!mcb_masks) { perror("all_bases_are_maps: mcb_masks"); exit(EXIT_FAILURE); }

        for (int c = 0; c < D; c++) {
            uint64_t mask = 0;
            for (int e = 0; e < g->nb_edges; e++) {
                if (!g->edges[e].deleted && mb->cycles[c].edges_ids[e])
                    mask |= (1ULL << e);
            }
            mcb_masks[c] = mask;
        }

        /* Cherche un plongement dont les faces internes = cycles de b. */
        int base_has_matching_embedding = 0;

        for (int ei = 0; ei < embs->count && !base_has_matching_embedding; ei++) {

            const FaceList *faces = &all_face_lists[ei];
            if (faces->count != D + 1) continue;

            /* ── Critère 1 : poids total des faces internes ───────────── */
            for (int skip = 0; skip < faces->count && !base_has_matching_embedding; skip++) {

                int weight_ok = 0;
                {
                    int total = 0;
                    for (int f = 0; f < faces->count; f++) {
                        if (f == skip) continue;
                        uint64_t m = faces->edges_masks[f];
                        int cnt = 0;
                        while (m) { m &= m - 1; cnt++; }
                        total += cnt;
                    }
                    weight_ok = (total == mcb_weight);
                }
                if (!weight_ok) continue;

                /* ── Critère 2 : bijection exacte entre masques ─────────
                 *
                 * On vérifie que { edges_masks[f] | f ≠ skip }
                 *              = { mcb_masks[c]   | c = 0..D-1 }.
                 * Chaque masque de face doit apparaître exactement une fois
                 * dans les masques de cycles (pas de doublon autorisé).    */
                int matched_cycle[64]; /* matched_cycle[c] = 1 si mcb_masks[c] déjà couplé */
                for (int c = 0; c < D; c++) matched_cycle[c] = 0;

                int bijection_ok = 1;
                for (int f = 0; f < faces->count && bijection_ok; f++) {
                    if (f == skip) continue;
                    int found = 0;
                    for (int c = 0; c < D && !found; c++) {
                        if (!matched_cycle[c] &&
                            faces->edges_masks[f] == mcb_masks[c]) {
                            matched_cycle[c] = 1;
                            found = 1;
                        }
                    }
                    if (!found) bijection_ok = 0;
                }

                if (bijection_ok) base_has_matching_embedding = 1;
            }
        }

        if (!base_has_matching_embedding) {
            /* Cette base ne correspond à aucun plongement. */
            all_match = 0;
        }

        free(mcb_masks);
    }

    /* Libération des FaceList précalculées. */
    for (int ei = 0; ei < embs->count; ei++)
        fml_free(&all_face_lists[ei]);
    free(all_face_lists);

    return all_match;
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

        if (g->nb_edges < min_edges) { delete_graph(g); trial--; continue; }

        if (!is_biconnected(g)) { delete_graph(g); trial--; continue; }

        /* ── Step 1 : minimal cycle bases ──────────────────────────── */
        int *inv_e = calloc(g->nb_edges,    sizeof(int));
        int *inv_v = calloc(g->nb_vertices, sizeof(int));
        if (!inv_e || !inv_v) {
            free(inv_e); free(inv_v); delete_graph(g); trial--; continue;
        }

        multiple_horton(g, inv_e, inv_v, 5000);
        free(inv_e); free(inv_v);

        if (g->nb_minimal_bases == 0 || g->minimals_basis[0].dimension == 0 || g->face_basis !=-1 || g->nb_face_basis_outer > 0) {
            delete_graph(g); trial--; continue;
        }

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

        /* ── Step 3 : enumerate embeddings and trace faces ───────────── */
        const EmbeddingSet embs = enumerate_embeddings(dfs);
        const int good_embedding_nb = verify_embedding_count(dfs, embs.count);
        if (embs.count == 0 || !good_embedding_nb) {
            free_embedding_set(&embs);   /* FIX : manquant → fuite à chaque skip */
            free_dfs_graph(dfs);
            delete_graph(g);
            trial--;
            continue;
        }

        /* Juste après enumerate_embeddings :
 * Compter combien de plongements ont Euler = 2. */
        // int euler_ok_count = 0;
        // for (int ei = 0; ei < embs.count; ei++) {
        //     FaceList faces_tmp = trace_faces(dfs, &embs.embeddings[ei]);
        //     /* trace_faces_silent = trace_faces sans fprintf, juste le compte */
        //     if (dfs->vertices_count - dfs->edge_count + faces_tmp.count == 2)
        //         euler_ok_count++;
        //     fml_free(&faces_tmp);
        // }
        // if (euler_ok_count == 0) {
        //     fprintf(stderr,
        //         "*** TOUS LES PLONGEMENTS NON-PLANAIRES "
        //         "(n=%d m=%d trial=%d) — labels corrompus ? ***\n",
        //         g->nb_vertices, g->nb_edges, trial);
        //
        //     /* Test d'isolation : reconstruire le DFS depuis zéro
        //      * sur le même graphe et vérifier les labels. */
        //     DfsGraph *dfs2 = build_dfs_graph(g);
        //     compute_low_values(dfs2);
        //     build_phi_lists(dfs2);
        //     build_singular_sets(dfs2);
        //     compute_same_diff(dfs2);
        //     build_SAME_DIFF_prime(dfs2);
        //     EmbeddingSet embs2 = enumerate_embeddings(dfs2);
        //
        //     int euler_ok2 = 0;
        //     for (int ei = 0; ei < embs2.count; ei++) {
        //         FaceList f2 = trace_faces(dfs2, &embs2.embeddings[ei]);
        //         if (dfs2->vertices_count - dfs2->edge_count + f2.count == 2)
        //             euler_ok2++;
        //         fml_free(&f2);
        //     }
        //     fprintf(stderr,
        //         "    Reconstruction iso : %d/%d plongements planaires\n",
        //         euler_ok2, embs2.count);
        //
        //     free_embedding_set(&embs2);
        //     free_dfs_graph(dfs2);
        // }

        printf("\ntrial=%d/%d embeddings=%d\n", trial + 1, nb_tests, embs.count);

        const int mcb_weight = compute_mcb_weight(g);

        FaceList all_faces;
        face_list_init(&all_faces);
        int mcb_matches = 0;
        int is_outerplanar = 0;

        const uint64_t all_verts = (g->nb_vertices < 64) ? ((1ULL << g->nb_vertices) - 1) : ~(uint64_t)0;

        int error = 0;
        int nb_map_is_bcm = 0;

        FaceList *cached_faces = malloc(embs.count * sizeof(FaceList));
        for (int ei = 0; ei < embs.count; ei++)
            cached_faces[ei] = trace_faces(dfs, &embs.embeddings[ei]);
            printf("trace_faces returned count=%d\n", cached_faces->count);

        for (int ei = 0; ei < embs.count; ei++) {

            FaceList faces = cached_faces[ei];

            printf("embedding %d : faces=%d D=%d\n",
       ei,
       faces.count,
       g->basis_dimension);

            printf("before MCB count=%d\n", faces.count);

            printf("before MCB ptr=%p count=%d\n", (void*)&faces, faces.count);

            int expected_faces = dfs->edge_count - dfs->vertices_count + 2;

            if (faces.count != expected_faces) {
                fml_free(&faces);
                continue;
            }

            int face_are_mcb = embedding_faces_are_mcb(&faces, g->basis_dimension, mcb_weight);
            if (face_are_mcb == -1) { error = 1; break; }
            nb_map_is_bcm += face_are_mcb;
            if (face_are_mcb > 0) mcb_matches = 1;

            if (!is_outerplanar) {
                for (int f = 0; f < faces.count; f++) {
                    if (faces.vertices_masks[f] == all_verts) {
                        is_outerplanar = 1;
                        mcb_matches = 1;
                        break;}
                }
            }

            for (int f = 0; f < faces.count; f++)
                face_list_push(&all_faces, faces.vertices_masks[f], faces.edges_masks[f]);
        }

        if (error) {
            free_embedding_set(&embs);
            free_dfs_graph(dfs);
            delete_graph(g);
            trial--;
            continue;
        }

        // const int nb_no_cofacial = count_minimal_never_cofacial(g->nb_vertices, &all_faces);
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
        const int all_bcm_are_maps = (nb_bases_with_map == g->nb_minimal_bases);


        char filename[1024];
        snprintf(filename, sizeof(filename), "%d_v(%d)_e(%d)_outer(%d)_mapMCB(%d)_allMapsMCB(%d)"
                                             "_MCBAlwaysMap(%d)_nbMCB(%d)_nbMCBMap(%d)_nbBaseswMap(%d)_nbMap(%d)_ratio(%d)_triangule(%d)",
                                             *id, dfs->vertices_count,
                                             dfs->edge_count, is_outerplanar, mcb_matches,
                                             nb_map_is_bcm == embs.count*g->basis_dimension, all_bcm_are_maps, g->nb_minimal_bases,
                                             nb_map_is_bcm, nb_bases_with_map, embs.count*g->basis_dimension+1, (int)(100.0f * nb_map_is_bcm / (float)(embs.count * (g->basis_dimension + 1))),
                                             g->nb_edges == 3*g->nb_vertices-6);


        for (int mb = 0; mb < g->nb_minimal_bases; mb++) {

            if (!basis_edge_at_most_2_cycles(g, &g->minimals_basis[mb]))
                continue;

            int found_embedding = 0;

            if (!mcb_matches) {
                save_to_folder("no_mcb_matches_but_oneMCB_2_edges_max", g, filename);
            }

            for (int ei = 0; ei < embs.count; ei++) {
                const FaceList faces = cached_faces[ei];

                if (mcb_matches_embedding(g, &g->minimals_basis[mb], faces)) {
                    found_embedding = 1;
                    break;
                }
            }

            if (!found_embedding) {
                save_to_folder("BCM_2_edges_not_map", g, filename);
                break;
            }
        }

        if (!mcb_matches) {
            save_to_folder("BCM_no_face", g, filename);
            const int witness_edge = has_four_chain_components(g);
            if (witness_edge == -1) save_to_folder("BCM_contre_exemple", g, filename);
            if (edge_max_2_cycle_all_bcm) save_to_folder("BCM_no_face_interesting", g, filename);
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

        if (nb_map_is_bcm > 1) save_to_folder("several_maps_are_BCM", g, filename);

        // if (nb_no_cofacial == 0 && !is_outerplanar) save_to_folder("pour_chloe", g, filename);

        fprintf(f,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
        *id,
        dfs->vertices_count,
        dfs->edge_count,
        is_outerplanar,
        mcb_matches,
        edge_max_2_cycle_all_bcm,
        (nb_map_is_bcm == embs.count * (g->basis_dimension + 1)),
        all_bcm_are_maps,
        g->nb_minimal_bases,
        nb_map_is_bcm,
        embs.count * (g->basis_dimension + 1),
        (int)(100.0f * nb_map_is_bcm / (float)(embs.count * (g->basis_dimension + 1))),
        g->nb_edges == 3*g->nb_vertices-6);

        fflush(f);

        (*id)++;

        /* ── Cleanup ──────────────────────────────────────────────── */
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

    FILE *f = fopen("graphs.csv", "a");

    // fprintf(f, "id,|V|,|E|,is_outerplanar,faces_can_be_MCB, all_bases_edges_at_most_2, "
    //            "all_maps_are_MCB, MCB_is_always_map, nb_BCM, nb_map_is_BCM, nb_map, nb_map_is_BCM/nb_map, triangule\n");

    int id = 0;
    for (int n = 15; n <= 30; n++) {
        printf( "n = %d ...\n", n);
        study_graph(n, 50, &id, f);
    }

    fclose(f);

    printf( "Done. %d graphs processed.\n", id);
    return 0;
}