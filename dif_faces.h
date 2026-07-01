#ifndef DIF_FACES_H
#define DIF_FACES_H

#include <stdint.h>
#include "dfs_graph.h"

/**
 * @brief Dynamic array of bitmasks representing the faces of a planar embedding.
 *
 * Each entry at index @c i encodes one face:
 *  - @c vertices_masks[i] : bitmask of vertices on that face (bit v = 1 iff vertex v is on the face).
 *  - @c edges_masks[i]    : bitmask of edges on that face (bit e = 1 iff edge e is on the face).
 *
 * Both arrays are kept strictly in sync: they always have the same @c count and
 * the same @c capacity, and are grown together by @c fml_push().
 *
 * The bitmask encoding makes subset testing O(1):
 *   S ⊆ face  iff  (face_mask & S_mask) == S_mask.
 */
typedef struct {
    uint64_t *vertices_masks; /**< Vertex-set bitmask, one per face.  */
    uint64_t *edges_masks; /**< Edge-set bitmask, one per face.    */
    int count; /**< Number of faces currently stored.  */
    int capacity; /**< Allocated capacity (doubled on overflow). */
} FaceList;

/**
 * @brief Context passed through the recursive never-co-facial search.
 */
typedef struct {
    const FaceList *faces; /**< All faces accumulated from every embedding. */
    uint64_t *minimal; /**< Bitmasks of minimal never-co-facial sets.   */
    int nb_minimal; /**< Number of minimal sets found so far.        */
    int cap_minimal; /**< Allocated capacity of @c minimal[].         */
} SearchCtx;

/* ── FaceList lifecycle ─────────────────────────────────────────────────── */

void face_list_init(FaceList *l);

/**
 * @brief Atomically append one face (vertex mask + edge mask) to a FaceList.
 *
 * Both internal arrays (@c vertices_masks and @c edges_masks) are grown
 * together when the capacity is exceeded, and @c count is incremented
 * exactly once. This is the only way to add a face; never call the old
 * split push_vertices / push_edges pair.
 *
 * @param l Target FaceList.
 * @param v_mask Vertex-set bitmask for the face.
 * @param e_mask Edge-set bitmask for the face (pass 0 if not needed).
 */
void face_list_push(FaceList *l, uint64_t v_mask, uint64_t e_mask);

void fml_free(FaceList *l);

/* ── Rotation-system helpers ────────────────────────────────────────────── */

int other_endpoint(const DfsGraph *dfs, int e, int v);
int arc_position(const PlanarEmbedding *emb, int v, int e);

/* ── Face tracing ───────────────────────────────────────────────────────── */

FaceList trace_faces(const DfsGraph *dfs, const PlanarEmbedding *emb);
void print_face_mask(uint64_t mask, int nb_vertices);

/* ── Never-co-facial search ─────────────────────────────────────────────── */

int  is_valid(uint64_t candidate, const FaceList *faces);
int  has_minimal_subset(uint64_t candidate, const uint64_t *minimal, int nb_minimal);
void record_minimal(SearchCtx *ctx, uint64_t mask);
void search_rec(int n, int start, int remaining, uint64_t current, SearchCtx *ctx);
int  popcount64(uint64_t x);
int  cmp_u64(const void *a, const void *b);
void find_minimal_never_cofacial(int nb_vertices, int nb_edges, const FaceList *all_faces);

#endif /* DIF_FACES_H */