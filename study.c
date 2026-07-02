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
 * @brief Return 1 iff the internal faces of this embedding form a minimal
 *        cycle basis.
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
        printf("embedding_faces_are_mcb: face count %d != basis dimension + 1 (%d)\n",
               faces->count, D + 1);
        return -1;
    }   /* Euler sanity check */

    int max = 0;
    int total = 0;
    int skipped = -1;
    for (int f = 0; f < faces->count; f++) {
        int lenght_face = popcount64(faces->edges_masks[f]);
        total += lenght_face;
        if (lenght_face > max) {
            max = lenght_face;
            skipped = f;
        }
    }

    total -= max;

    if (total == mcb_weight) return 1;
    // printf("%d %d\n", total, mcb_weight);
    // for (int f = 0; f < faces->count; f++) {
    //     if (f == skipped) continue;
    //     print_face_mask(faces->vertices_masks[f], 64);
    // }
    return 0;
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
static int mcb_matches_embedding(const Graph *g, const Minimal_basis *mb, const FaceList *faces) {

    int D = mb->dimension;
    if (faces->count != D + 1) return 0;

    uint64_t *mcb = malloc(D * sizeof(uint64_t));
    if (!mcb) { perror("mcb_matches_embedding"); exit(EXIT_FAILURE); }
    for (int c = 0; c < D; c++)
        mcb[c] = cycle_to_edges_mask(g, &mb->cycles[c]);

    int result = 0;
    for (int skip = 0; skip < faces->count && !result; skip++) {
        int all_found = 1;
        for (int c = 0; c < D && all_found; c++) {
            int found = 0;
            for (int f = 0; f < faces->count && !found; f++) {
                if (f == skip) continue;
                if (faces->edges_masks[f] == mcb[c]) found = 1;
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

static void save_to_folder(const char *folder, const Graph *g, int id) {
    ensure_dir(folder);
    char path[1024];
    snprintf(path, sizeof(path),
             "%s/graph_%dBCM_%dv_%de_id%d.txt",
             folder, g->nb_minimal_bases, g->nb_vertices, g->nb_edges, id);
    save_graph(g, path);
}

/* =========================================================================
 * Core study loop
 * ========================================================================= */

void study_graph(const int nb_vertex, const int nb_tests, int *id, FILE *f) {

    const int min_edges = nb_vertex+1;
    const int max_edges = 3 * nb_vertex - 6 > 64 ? 64 : 3 * nb_vertex - 6;
    if (min_edges > max_edges) return;

    for (int trial = 0; trial < nb_tests; trial++) {

        /* ── Generate ──────────────────────────────────────────────── */
        const int target_edges =
            min_edges + rand() % (max_edges - min_edges + 1);

        Graph *g = create_planar_graph(nb_vertex, target_edges);
    if (!g) {
        trial--;
        continue;
    }

    /* ── Vérification d'intégrité post-génération ──────────────────────
    * Teste directement l'hypothèse delete_vertex/compact_graph :
    * la somme des degrés doit valoir exactement 2*nb_edges, et chaque
    * arête de g->edges[] doit apparaître exactement deux fois dans les
    * listes de voisins (une fois par extrémité), sans doublon ni absence.
    */

    int sum_degrees = 0;
    for (int v = 0; v < g->nb_vertices; v++)
        sum_degrees += g->neighbors[v].count;

    if (sum_degrees != 2 * g->nb_edges) {
        fprintf(stderr,
            "*** INTEGRITY FAIL (id=%d): sum_degrees=%d != 2*nb_edges=%d ***\n",
            *id, sum_degrees, 2 * g->nb_edges);
        abort();
    }

    /* Vérifie que chaque arête (u,v) de g->edges[] est bien présente
     * dans neighbors[u] ET neighbors[v], et réciproquement que chaque
     * entrée de neighbors[] correspond à une vraie arête existante. */
    for (int e = 0; e < g->nb_edges; e++) {
        if (g->edges[e].deleted) {trial--; continue;}
        int u = g->edges[e].u, v = g->edges[e].v;

        int found_u_in_v = 0, found_v_in_u = 0;
        for (int i = 0; i < g->neighbors[u].count; i++)
            if (g->neighbors[u].neighbors[i] == v) found_v_in_u = 1;
        for (int i = 0; i < g->neighbors[v].count; i++)
            if (g->neighbors[v].neighbors[i] == u) found_u_in_v = 1;

        if (!found_v_in_u || !found_u_in_v) {
            fprintf(stderr,
                "*** INTEGRITY FAIL (id=%d): edge %d (%d,%d) missing from "
                "neighbor lists (found_v_in_u=%d found_u_in_v=%d) ***\n",
                *id, e, u, v, found_v_in_u, found_u_in_v);
            abort();
        }
    }

    /* Vérifie l'inverse : chaque entrée de neighbors[v] correspond à
     * une arête réellement existante dans g->edges[] (pas de "fantôme"). */
    for (int v = 0; v < g->nb_vertices; v++) {
        for (int i = 0; i < g->neighbors[v].count; i++) {
            int w = g->neighbors[v].neighbors[i];
            int found = 0;
            for (int e = 0; e < g->nb_edges; e++) {
                if (g->edges[e].deleted) {trial--; continue;}
                if ((g->edges[e].u == v && g->edges[e].v == w) ||
                    (g->edges[e].u == w && g->edges[e].v == v)) {
                    found = 1; break;
                }
            }
            if (!found) {
                fprintf(stderr,
                    "*** INTEGRITY FAIL (id=%d): neighbors[%d] contains "
                    "phantom edge to %d with no matching g->edges[] entry ***\n",
                    *id, v, w);
                abort();
            }
        }
    }
    if (g->nb_edges < min_edges) { delete_graph(g); trial--; continue; }

    /* 2-connectivity is required by Cai's algorithm.
     * Check BEFORE Horton to avoid wasting time on non-2-connected graphs. */
    if (!is_biconnected(g)) { delete_graph(g); trial--; continue; }

    /* ── Step 1 : minimal cycle bases ──────────────────────────── */
    int *inv_e = calloc(g->nb_edges,    sizeof(int));
    int *inv_v = calloc(g->nb_vertices, sizeof(int));
    if (!inv_e || !inv_v) {
        free(inv_e); free(inv_v); delete_graph(g); trial--; continue;
    }
    multiple_horton(g, inv_e, inv_v, 1000);
    free(inv_e); free(inv_v);

    if (g->nb_minimal_bases == 0 || g->minimals_basis[0].dimension == 0 || g->face_basis !=-1 || g->nb_face_basis_outer > 0) {
        delete_graph(g); trial--; continue;
    }

    const int expected_dim = g->nb_edges - g->nb_vertices + 1;
    if (g->minimals_basis[0].dimension != expected_dim) {
        delete_graph(g); trial--; continue;
    }

    printf("n=%d m=%d id=%d", g->nb_vertices, g->nb_edges, *id);

    // print neighbors
    // printf("=== Neighbors original graph ===");
    // for (int v=0; v < g->nb_vertices; v++) {
    //     printf("\n%d:", v);
    //     for (int n=0; n < g->neighbors[v].count; n++) {
    //         printf(" %d", g->neighbors[v].neighbors[n]);
    //     }
    //     printf("\n");
    // }

    /* ── Step 1b : MCB edge-multiplicity filter ─────────────────── */
    const int edges_ok = all_bases_edge_at_most_2(g);

    /* ── Step 2 : build DFS structures ──────────────────────────── */
    DfsGraph *dfs = build_dfs_graph(g);

    // printf( "=== DFS parent_edge check ===\n");
    // for (int v = 0; v < dfs->vertices_count; v++) {
    //     int pe = dfs->vertices[v].parent_edge;
    //     if (pe == -1) {
    //         printf( "  v=%d dfs_num=%d  root\n",
    //                 v, dfs->vertices[v].dfs_num);
    //     } else {
    //         int from = dfs->edges[pe].from;
    //         int to = dfs->edges[pe].to;
    //         int ok = (from == v || to == v);
    //         printf("  v=%d dfs_num=%d  parent_edge=%d (%d -> %d) %s\n",
    //                 v, dfs->vertices[v].dfs_num,
    //                 pe, from, to, ok ? "OK" : "*** WRONG ***");
    //     }
    // }

    compute_low_values(dfs);
    build_phi_lists(dfs);
    // printf( "=== out_edges (Phi(v)) ===\n");
    // for (int v = 0; v < dfs->vertices_count; v++) {
    //     printf( "  Phi(%d): ", v);
    //     for (int i = 0; i < dfs->vertices[v].nb_out_edges; i++) {
    //         int e = dfs->vertices[v].out_edges[i];
    //         printf( "arc%d(%d -> %d) ", e,
    //                 dfs->edges[e].from, dfs->edges[e].to);
    //     }
    //     printf( "\n");
    // }

    build_singular_sets(dfs);
    compute_same_diff(dfs);
    build_SAME_DIFF_prime(dfs);

    // Vérifie que chaque voisin est bien un ID valide
    for (int i = 0; i < g->nb_vertices; i++) {
        for (int j = 0; j < g->neighbors[i].count; j++) {
            int neighbor_id = g->neighbors[i].neighbors[j];
            if (neighbor_id < 0 || neighbor_id >= g->nb_vertices) {
                printf("ERREUR: ID invalide %d trouvé dans les voisins de %d\n", neighbor_id, i);
                exit(1);
            }
        }
    }

    /* ── Step 3 : enumerate embeddings and trace faces ───────────── */
    const EmbeddingSet embs = enumerate_embeddings(dfs);
    int good_embedding_nb = verify_embedding_count(dfs, embs.count);
    if (embs.count == 0 || !good_embedding_nb) {
        free_embedding_set((EmbeddingSet *)&embs);  /* FIX : fuite mémoire */
        free_dfs_graph(dfs);
        delete_graph(g);
        trial--;
        continue;
    }

    /* DEBUG: validate every arc in every rotation */
    for (int ei = 0; ei < embs.count; ei++) {
        for (int v = 0; v < dfs->vertices_count; v++) {
            for (int ai = 0; ai < embs.embeddings[ei].rotation_lengths[v]; ai++) {
                int arc = embs.embeddings[ei].rotation_system[v][ai];
                int from = dfs->edges[arc].from;
                int to = dfs->edges[arc].to;
                if (from != v && to != v)
                    printf(
                        "*** INVALID emb[%d] rot[%d][%d]=arc%d (%d→%d) NOT incident to %d!\n",
                        ei, v, ai, arc, from, to, v);
            }
        }
    }

    /* À ajouter juste après EmbeddingSet embs = enumerate_embeddings(dfs); dans study.c */

    /* Validation stricte : chaque arc de chaque rotation doit être incident
     * au sommet qui le porte. Si ce n'est jamais le cas, le bug n'est PAS
     * dans trace_faces mais dans enumerate_embeddings (rotation système
     * corrompu en amont, donc heap corruption toujours présente). */
    int rotation_corrupted = 0;
    for (int ei = 0; ei < embs.count; ei++) {
        PlanarEmbedding *emb = &embs.embeddings[ei];
        for (int v = 0; v < emb->n; v++) {
            for (int ai = 0; ai < emb->rotation_lengths[v]; ai++) {
                int arc = emb->rotation_system[v][ai];
                if (arc < 0 || arc >= dfs->edge_count) {
                    fprintf(stderr,
                        "*** emb[%d] v=%d ai=%d: arc=%d HORS LIMITES (edge_count=%d) ***\n",
                        ei, v, ai, arc, dfs->edge_count);
                    rotation_corrupted = 1;
                    trial--;
                    continue;
                }
                int from = dfs->edges[arc].from;
                int to   = dfs->edges[arc].to;
                if (from != v && to != v) {
                    fprintf(stderr,
                        "*** emb[%d] v=%d ai=%d: arc%d (%d->%d) NON INCIDENT à %d ***\n",
                        ei, v, ai, arc, from, to, v);
                    rotation_corrupted = 1;
                }
            }
        }
    }
    if (rotation_corrupted) {
        printf("*** ROTATION SYSTEM CORROMPU pour ce graphe (id=%d) ***\n", *id);
    }

    printf("\ntrial=%d/%d embeddings=%d\n",
            trial + 1, nb_tests, embs.count);

    /* Pre-compute MCB weight once: sum of cycle lengths in any MCB.
     * Used by embedding_faces_are_mcb() to decide minimality. */
    const int mcb_weight = compute_mcb_weight(g);

    FaceList all_faces;
    face_list_init(&all_faces);
    int mcb_matches = 0;
    int is_outerplanar = 0;

    const uint64_t all_verts = (g->nb_vertices < 64)
                               ? ((1ULL << g->nb_vertices) - 1)
                               : ~(uint64_t)0;

    int error = 0;
    for (int ei = 0; ei < embs.count; ei++) {
        // printf("=== Embedding %d ===\n", ei + 1);
        FaceList faces = trace_faces(dfs, &embs.embeddings[ei]);

        if (!mcb_matches) {
            int face_are_mcb = embedding_faces_are_mcb(&faces, g->minimals_basis[0].dimension, mcb_weight);
            if (face_are_mcb)
                mcb_matches = 1;
            if (face_are_mcb == -1) {
                fml_free(&faces);
                error = 1;
                break;
            }
        }

        if (!is_outerplanar) {
            // printf("=== Checking for outerplanarity ===\n");
            for (int f = 0; f < faces.count; f++) {
                // print_face_mask(faces.vertices_masks[f], g->nb_vertices);
                if (faces.vertices_masks[f] == all_verts) {
                    is_outerplanar = 1;
                    mcb_matches = 1;
                    // printf("outerplanar");
                    break;}
            }
        }

        for (int f = 0; f < faces.count; f++)
            face_list_push(&all_faces,
                           faces.vertices_masks[f],
                           faces.edges_masks[f]);

        fml_free(&faces);
    }

    /* Si une erreur est survenue, on abandonne complètement ce graphe */
    if (error) {
        free_embedding_set(&embs);
        free_dfs_graph(dfs);
        delete_graph(g);
        trial--;
        continue;
    }

    /* ── Step 4 : never-co-facial search ───────────────────────── */
    const int nb_no_cofacial = count_minimal_never_cofacial(g->nb_vertices,
                                                       &all_faces);
    fml_free(&all_faces);

    /* ── Step 5 : save ─────────────────────────────────────────── */

    int edge_change_nb_appearence = 0;
    int diff_edge = 0;
    if (!mcb_matches) {
        save_to_folder("BCM_no_face", g, *id);
        edge_change_nb_appearence = edge_two_and_three_appearence(g, g->minimals_basis);
        diff_edge = diff_edge_occur(g, g->minimals_basis);
        if (edges_ok) save_to_folder("BCM_no_face_interesting", g, *id);
        if (edge_change_nb_appearence) save_to_folder("No_BCM_face_1_edge_can_appear_2_or_3", g, *id);
        if (diff_edge) save_to_folder("No_BCM_face_all_edges_appear_2_somtimes_3", g, *id);
    }

    if (nb_no_cofacial == 0 && !is_outerplanar)
        save_to_folder("pour_chloe", g, *id);

    fprintf(f,"%d,%d,%d,%d,%d,%d,%d,%d\n",
    *id,
    nb_vertex,
    dfs->edge_count,
    is_outerplanar,
    mcb_matches,
    edges_ok,
    edge_change_nb_appearence,
    diff_edge);

    (*id)++;

    /* ── Cleanup ──────────────────────────────────────────────── */
    free_embedding_set(&embs);
    free_dfs_graph(dfs);
    delete_graph(g);
    g = NULL;
    }
}

int main(void) {
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);
    srand((unsigned int)time(NULL));

    FILE *f = fopen("graphs.csv", "w");

    fprintf(f,
    "id,|V|,|E|,is_outerplanar,faces_can_be_MCB, max_appearence_edge_always_2,"
    "1_edge_can_appear_2_or_3,all_edges_appear_2_somtimes_3\n");

    int id = 0;
    for (int n = 6; n <= 14; n++) {
        printf( "n = %d ...\n", n);
        study_graph(n, 1000, &id, f);
    }

    fclose(f);

    printf( "Done. %d graphs processed.\n", id);
    return 0;
}