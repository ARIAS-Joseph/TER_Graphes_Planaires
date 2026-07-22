/**
 * @file dif_faces_vertices.c
 * @brief Find minimal vertex sets that never appear together on the same face across all planar
 * embeddings of a 2-connected planar graph of 64 edges maximum.
 *
 * This file implements two main features:
 *  1. Face tracing from a CW rotation system
 *  2. Search for minimal "never co-facial" vertex sets: sets S such that no single face of any
 *  planar embedding contains all vertices of S.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "dfs_graph.h"
#include "dif_faces.h"

/**
 * @brief Initialize an empty FaceList with an initial capacity of 16 faces.
 *
 * @param l Pointer to the FaceList to initialize. Must not be NULL.
 */
void face_list_init(FaceList *l) {
    l->count = 0;
    l->capacity = 16;
    l->vertices_masks = malloc(l->capacity * sizeof(uint64_t));
    l->edges_masks = malloc(l->capacity * sizeof(uint64_t));
    if (!l->vertices_masks || !l->edges_masks) {
        perror("fml_init");
        exit(EXIT_FAILURE);
    }
}


/**
 * @brief Append a bitmask to a FaceList, growing the array if necessary.
 *
 * @param l Pointer to the target FaceList.
 * @param v_mask Bitmask of vertices in the face.
 * @param e_mask Bitmask of edges in the face.
 */
void face_list_push(FaceList *l, const uint64_t v_mask, const uint64_t e_mask) {

    if (l->count >= l->capacity) {
        l->capacity *= 2;
        uint64_t *temp_v_mask = realloc(l->vertices_masks, l->capacity * sizeof(uint64_t));
        uint64_t *temp_e_mask = realloc(l->edges_masks, l->capacity * sizeof(uint64_t));
        if (!temp_v_mask || !temp_e_mask) {
            perror("face_list_push");
            exit(EXIT_FAILURE);
        }

        l->vertices_masks = temp_v_mask;
        l->edges_masks = temp_e_mask;
    }

    l->vertices_masks[l->count] = v_mask;
    l->edges_masks[l->count] = e_mask;
    l->count++;
}

/**
 * @brief Free the internal buffer of a FaceList and reset its fields.
 *
 * After this call, l->masks is NULL and l->count and l->capacity are both 0.
 *
 * @param l Pointer to the FaceList to free.
 */
void fml_free(FaceList *l) {
    free(l->vertices_masks);
    l->vertices_masks = NULL;
    free(l->edges_masks);
    l->edges_masks = NULL;
    l->count = l->capacity = 0;
}

/**
 * @brief Return the endpoint of DFS arc e that is not v.
 *
 * Each DFS arc has two endpoints stored as "from" and "to". Given one endpoint v, this function
 * returns the other one.
 *
 * @param dfs Pointer to the DFS graph.
 * @param e Index of the DFS arc.
 * @param v Known endpoint.
 * @return The opposite endpoint of v according to arc e.
 */
int other_endpoint(const DfsGraph *dfs, int e, int v)
{
    if (dfs->edges[e].from == v)
        return dfs->edges[e].to;

    if (dfs->edges[e].to == v)
        return dfs->edges[e].from;

    // printf("ERROR: vertex %d is not incident to edge %d (%d,%d)\n",
    //        v, e,
    //        dfs->edges[e].from,
    //        dfs->edges[e].to);
    abort();
}

/**
 * @brief Return the position of arc e in the rotation system of vertex v.
 *
 * The rotation system emb->rotation_system[v] lists all arcs incident to v in clockwise (CW) order.
 *
 * @param emb Pointer to the planar embedding.
 * @param v Vertex whose rotation list is used.
 * @param e Index of the arc to locate.
 * @return Index of e in rotation_system[v], or -1 if not found (should not happen on a valid
 * embedding).
 */
int arc_position(const PlanarEmbedding *emb, const int v, const int e) {
    for (int j = 0; j < emb->rotation_lengths[v]; j++) {
        if (emb->rotation_system[v][j] == e) {
            return j;
        }
    }
    return -1;
}

/**
 * @brief Trace all faces of a planar embedding using the vertices rotation system.
 *
 * The embedding is encoded as a clockwise rotation system: for each vertex v,
 * emb->rotation_system[v] lists the indices of all incident arcs in clockwise order around v.
 *
 * A dart (directed half-edge) is an arc seen from one of its endpoints: dart(v, ai) means "standing
 * at vertex v, taking the arc at position ai in my CW rotation list".
 *
 * Faces are the cycles of the face permutation phi defined by:
 *   phi(v, ai)  →  (w, (bi − 1 + deg_w) mod deg_w)
 *   where  e = rotation_system[v][ai] (arc taken from v)
 *          w = other endpoint of e (arrival vertex)
 *          bi = position of e in rotation_system[w]
 *          deg_w = rotation_lengths[w]
 *
 * Since every dart belongs to exactly one face, each dart is marked as visited once processed to
 * avoid retracing the same face.
 *
 * @param dfs Pointer to the DFS graph.
 * @param emb Pointer to the planar embedding whose faces are to be traced.
 * @return A FaceList containing two bitmasks per face found (on for the vertices and one for the
 * edges).
 */
FaceList trace_faces(const DfsGraph *dfs, const PlanarEmbedding *emb) {

    const int n = emb->n;
    const int m = dfs->edge_count;

    // for (int v = 0; v < n; v++) {
    //     printf("rotation[%d] :", v);
    //     for (int i = 0; i < emb->rotation_lengths[v]; i++)
    //         printf(" %d", emb->rotation_system[v][i]);
    //     printf("\n");
    // }

    /* ── Validation 1 : somme des degrés = 2m ──────────────────────────── */
    int total_rotation = 0;
    for (int v = 0; v < n; v++) total_rotation += emb->rotation_lengths[v];
    // if (total_rotation != 2 * m)
    //     printf(
    //         "[trace_faces] SOMME ROTATIONS %d != 2*m=%d\n",
    //         total_rotation, 2 * m);

    /* ── Validation 2 : chaque arc apparaît exactement à from ET à to ──── */
    /* (l'ancienne version remplissait arc_occurrences mais n'imprimait rien) */
    int *arc_at_from = calloc(m, sizeof(int));
    int *arc_at_to = calloc(m, sizeof(int));

    for (int v = 0; v < n; v++) {
        for (int ai = 0; ai < emb->rotation_lengths[v]; ai++) {
            int e = emb->rotation_system[v][ai];
            if (e < 0 || e >= m) {
                // printf(
                //     "[trace_faces] ARC HORS BORNES : rotation[%d][%d]=%d "
                //     "(edge_count=%d)\n", v, ai, e, m);
                continue;
            }
            if (dfs->edges[e].from == v) arc_at_from[e]++;
            else if (dfs->edges[e].to == v) arc_at_to[e]++;
            // else
                // printf(
                //     "[trace_faces] ARC MAL PLACE : arc %d (from=%d to=%d) "
                //     "dans rotation[%d]\n",
                //     e, dfs->edges[e].from, dfs->edges[e].to, v);
        }
    }

    // for (int e = 0; e < m; e++) {
    //     if (arc_at_from[e] != 1)
    //         printf(
    //             "[trace_faces] Arc %d (from=%d to=%d) : "
    //             "present %d fois côté FROM (attendu 1)\n",
    //             e, dfs->edges[e].from, dfs->edges[e].to, arc_at_from[e]);
    //     if (arc_at_to[e] != 1)
    //         printf(
    //             "[trace_faces] Arc %d (from=%d to=%d) : "
    //             "present %d fois côté TO (attendu 1)\n",
    //             e, dfs->edges[e].from, dfs->edges[e].to, arc_at_to[e]);
    // }

    free(arc_at_from);
    free(arc_at_to);

    /* ── Tracé des faces ────────────────────────────────────────────────── */
    int **visited = malloc(n * sizeof(int *));
    for (int v = 0; v < n; v++)
        visited[v] = calloc(emb->rotation_lengths[v], sizeof(int));

    FaceList faces;
    face_list_init(&faces);

    /* Garde-fou contre boucle infinie : chaque dart est visité au plus 1 fois */
    const int max_iters = 2 * m + 4;

    for (int start_v = 0; start_v < n; start_v++) {
        for (int start_ai = 0; start_ai < emb->rotation_lengths[start_v]; start_ai++) {

            if (visited[start_v][start_ai]) continue;

            uint64_t face_mask_vertices = 0;
            uint64_t face_mask_edges = 0;
            int v = start_v;
            int ai = start_ai;
            int iters = 0;
            int ok = 1;

            int face_len = 0;

            do {
                face_len++;
                if (++iters > max_iters) {
                    ok = 0;
                    break;
                }

                visited[v][ai] = 1;
                face_mask_vertices |= (1ULL << v);

                const int e = emb->rotation_system[v][ai];

                // printf("v=%d ai=%d e=%d from=%d to=%d\n",
                //        v,
                //        ai,
                //        e,
                //        dfs->edges[e].from,
                //        dfs->edges[e].to);

                const int w = other_endpoint(dfs,e,v);

                // printf("   w=%d\n", w);
                const int bi = arc_position(emb, w, e);

                if (bi < 0) {
                    // printf("edge %d not found in rotation of vertex %d\n", e, v);
                    abort();
                }

                if (bi == -1) {
                    ok = 0;
                    break;
                }

                const int deg_w = emb->rotation_lengths[w];
                const int next_ai = (bi + 1) % deg_w;

                int deg_w2 = 0;

                for (int e = 0; e < dfs->edge_count; e++) {
                    if (dfs->edges[e].from == w || dfs->edges[e].to == w) {
                        deg_w2++;
                    }
                }

                // printf("w=%d e=%d bi=%d deg rotation=%d def defs=%d\n", w, e, bi, deg_w, deg_w2);

                assert(bi >= 0);
                assert(bi < deg_w);
                assert(emb->rotation_system[w][bi] == e);
                assert(deg_w2 == deg_w);

                int idx = dfs->edge_indices[v * dfs->vertices_count + w];

                if (idx < 0 || idx >= 64) {
                    // printf("BAD EDGE INDEX: v=%d w=%d idx=%d\n",
                    //        v, w, idx);
                    abort();
                }

                face_mask_edges |= (1ULL << idx);

                // printf("dart=%d  %d->%d  edgeidx=%d\n", e, v, w, dfs->edge_indices[v*n+w]);

                v = w;
                ai = next_ai;

            } while (v != start_v || ai != start_ai);

            if (ok) {
                // printf("FACE %d length=%d\n", faces.count+1, popcount64(face_mask_edges));
                face_list_push(&faces, face_mask_vertices, face_mask_edges);
                // printf("FACE %d darts=%d edges=%d\n", faces.count + 1, face_len, popcount64(face_mask_edges));
            }
            /* Si !ok : on ne pousse pas la face tronquée, et les darts
             * non-visités seront repris par la boucle externe normalement. */
        }
    }

    int visited_count = 0;

    for (int v = 0; v < dfs->vertices_count; ++v)
        for (int s = 0; s < emb->rotation_lengths[v]; ++s)
            if (visited[v][s])
                visited_count++;

    // printf("visited darts = %d / %d\n", visited_count, 2 * dfs->edge_count);

    for (int v = 0; v < n; v++) free(visited[v]);
    free(visited);

    int total = 0;
    for (int v = 0; v < n; v++)
        total += emb->rotation_lengths[v];

    // printf("total darts=%d expected=%d\n", total, 2 * dfs->edge_count);

    int sum_face_lengths = 0;

    for (int f = 0; f < faces.count; f++) {
        uint64_t m = faces.edges_masks[f];

        int len = 0;
        while (m) {
            m &= m - 1;
            len++;
        }

        sum_face_lengths += len;
    }

    // printf("faces=%d sum=%d expected=%d\n", faces.count, sum_face_lengths,2 * dfs->edge_count);

    int total2 = 0;
    // for (int i = 0; i < faces.count; i++) {
    //     printf("stored face %d length=%d\n", i,
    //            popcount64(faces.edges_masks[i]),
    //     total2 += popcount64(faces.edges_masks[i]));
    // }
    // printf("total=%d\n", total2);

    int expected_faces = dfs->edge_count - dfs->vertices_count + 2;

    // if (faces.count != expected_faces) {
    //     printf("INVALID EMBEDDING : faces=%d expected=%d\n",
    //            faces.count,
    //            expected_faces);
    // }

    return faces;
}

/**
 * @brief Print the vertices of a face bitmask (enclosed in braces).
 *
 * @param mask Bitmask of the face.
 * @param nb_vertices Total number of vertices in the graph.
 */
void print_face_mask(const uint64_t mask, const int nb_vertices) {
    printf("{ ");
    for (int v = 0; v < nb_vertices; v++)
        if (mask & (1ULL << v)) printf("%d ", v);
    printf("}\n");
}

/**
 * @brief Test if the candidate vertex set is never entirely contained in any recorded face.
 *
 * @param candidate Bitmask of the vertex set to test.
 * @param faces All face bitmasks accumulated over every embedding.
 * @return 1 if the candidate set is valid (never co-facial), 0 otherwise.
 */
int is_valid(const uint64_t candidate, const FaceList *faces) {
    for (int i = 0; i < faces->count; i++)
        if ((faces->vertices_masks[i] & candidate) == candidate) {
            return 0;
        }
    return 1;
}

/**
 * @brief Test whether candidate contains a known minimal set as a subset.
 *
 * @param candidate Bitmask of the set to test.
 * @param minimal Array of already found minimal vertices bitmasks.
 * @param nb_minimal Number of entries in minimal.
 * @return 1 if candidate contains at least one known minimal subset, 0 otherwise.
 */
int has_minimal_subset(const uint64_t candidate, const uint64_t *minimal, const int nb_minimal) {
    for (int i = 0; i < nb_minimal; i++)
        if ((candidate & minimal[i]) == minimal[i]) return 1;
    return 0;
}

/**
 * @brief Record a new minimal never-co-facial set.
 *
 * @param ctx Search context holding the minimal set list.
 * @param mask Bitmask of the new minimal set to record.
 */
void record_minimal(SearchCtx *ctx, const uint64_t mask) {

    if (ctx->nb_minimal >= ctx->cap_minimal) {
        ctx->cap_minimal = ctx->cap_minimal ? ctx->cap_minimal * 2 : 8;
        uint64_t *temp = realloc(ctx->minimal, ctx->cap_minimal * sizeof(uint64_t));
        if (!temp) {
            perror("record_minimal");
            exit(EXIT_FAILURE);
        }
        ctx->minimal = temp;
    }

    ctx->minimal[ctx->nb_minimal++] = mask;
}

/**
 * @brief Recursively enumerate all subsets of size "remaining" from vertices in [start, n], testing
 * each subset for minimality.
 *
 * Subsets are built by choosing vertices in increasing order, so each combination is visited
 * exactly once.
 *
 * As soon as the partial set "next" already contains a known minimal subset, the entire subtree
 * rooted at "next" is skipped.
 *
 * At a leaf (remaining == 0), "current" is a complete candidate of the target size. It is recorded
 * as a new minimal set only if:
 * - no subset of "current" is already a known minimal;
 * - it is never entirely contained in any face.
 *
 * @param n Total number of vertices.
 * @param start Lowest vertex index still available for selection.
 * @param remaining Number of vertices still to be added to reach target size.
 * @param current Bitmask of vertices chosen so far.
 * @param ctx Search context (all faces of all embeddings + no-co-facial sets).
 */
void search_rec(const int n, const int start, const int remaining, const uint64_t current,
    SearchCtx *ctx) {

    if (remaining == 0) {
        if (!has_minimal_subset(current, ctx->minimal, ctx->nb_minimal)
            && is_valid(current, ctx->faces)) {

            record_minimal(ctx, current);
        }

        return;
    }

    /* Upper bound on v: at least "remaining" vertices must be left after v, so v <= n - remaining.
     */
    for (int v = start; v <= n - remaining; v++) {
        const uint64_t next = current | (1ULL << v);
        if (has_minimal_subset(next, ctx->minimal, ctx->nb_minimal)) {
            continue; /* early stop: whole subtree is non-minimal */
        }
        search_rec(n, v + 1, remaining - 1, next, ctx);
    }
}

/**
 * @brief Find and print all minimal vertex sets that are never co-facial across any of the planar
 * embeddings recorded in all_faces.
 *
 * A vertex set S is valid if there is no face F of any embedding where S is entirely contained in
 * F.
 *
 * A valid set S is minimal if no subset of S of size > 1 is already valid.
 *
 * The search iterates over subset sizes k = 2, 3, ... and uses search_rec() with early stop to
 * avoid revisiting supersets of already found minimals sets.
 *
 * @param nb_vertices Number of vertices in the graph (must be < 65).
 * @param nb_edges Number of edges in the graph (must be < 65).
 * @param all_faces All FaceList from all planar embeddings.
 */
void find_minimal_never_cofacial(const int nb_vertices, const int nb_edges, const FaceList *all_faces) {
    if (nb_vertices < 2) {
        // printf("\nFewer than 2 vertices: search not applicable.\n"); return;
    }
    if (nb_vertices > 64 || nb_edges > 64) {
        // printf("\nMore than 64 vertices and/or edges: uint64_t encoding is insufficient.\n");
        return;
    }

    SearchCtx ctx = { all_faces, NULL, 0, 0 };

    for (int k = 2; k <= nb_vertices; k++) {
        search_rec(nb_vertices, 0, k, 0, &ctx);
    }

    // printf("\n");
    // if (ctx.nb_minimal == 0) {
    //     printf("No vertex set of size >= 2 satisfies the condition");
    // } else {
    //     printf("%d minimal set(s) found :\n", ctx.nb_minimal);
    //     for (int i = 0; i < ctx.nb_minimal; i++) {
    //         print_face_mask(ctx.minimal[i], nb_vertices);
    //     }
    // }

    free(ctx.minimal);
}


// int main(void) {
//
//     /* Example graph */
//     Graph *g = create_graph();
//     for (int i = 0; i < 4; ++i) create_vertex(g, 0.0, 0.0);
//
//     create_edge(g, 0, 1);
//     create_edge(g, 0, 3);
//     create_edge(g, 1, 2);
//     create_edge(g, 2, 3);
//
//
//     // create_edge(g, 0, 1);
//     // create_edge(g, 0, 2);
//     // create_edge(g, 0, 4);
//     // create_edge(g, 0, 5);
//     // create_edge(g, 0,3);
//     // create_edge(g, 1, 2);
//     // create_edge(g, 2,6);
//     // create_edge(g, 2, 4);
//     // create_edge(g, 3, 4);
//     // create_edge(g, 3, 6);
//     // create_edge(g, 4,5);
//
//     /* Pipeline to compute all embeddings of G */
//     DfsGraph *dfs = build_dfs_graph(g);
//     compute_low_values(dfs);
//     build_phi_lists(dfs);
//     build_singular_sets(dfs);
//     compute_same_diff(dfs);
//     build_SAME_DIFF_prime(dfs);
//     const EmbeddingSet result = enumerate_embeddings(dfs);
//
//     FaceList all_faces;
//     face_list_init(&all_faces);
//
//     /* Creation of the expected outer face if G is outerplanar */
//     const uint64_t all_verts = ((1ULL << g->nb_vertices) - 1);
//     int is_outerplanar = 0;
//
//     for (int i = 0; i < result.count; i++) { /* for every embedding of G */
//         printf("=== Embedding %d ===\n", i + 1);
//
//         FaceList faces = trace_faces(dfs, &result.embeddings[i]);
//
//         for (int f = 0; f < faces.count; f++) {
//             printf("  Face %d : ", f);
//             print_face_mask(faces.vertices_masks[f], dfs->vertices_count);
//             /* Copy the face vertices mask into a global list */
//             face_list_push(&all_faces, faces.vertices_masks[f], faces.edges_masks[f]);
//         }
//         printf("\n");
//
//         if (!is_outerplanar) {
//             for (int f = 0; f < faces.count; f++) {
//                 if (faces.vertices_masks[f] == all_verts) {
//                     is_outerplanar = 1;
//                     break;
//                 }
//             }
//         }
//
//         fml_free(&faces);
//     }
//
//     printf("%s\n", is_outerplanar ? "outerplanar" : "not outerplanar");
//
//     find_minimal_never_cofacial(dfs->vertices_count, dfs->edge_count, &all_faces);
//
//     fml_free(&all_faces);
//     free_embedding_set((EmbeddingSet *)&result);
//     free_dfs_graph(dfs);
//     delete_graph(g);
//
//     return 0;
// }