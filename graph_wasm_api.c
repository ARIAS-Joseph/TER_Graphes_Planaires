#include <emscripten.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "graph.h"
#include "planar_graph_creator.h"


static Graph *g = NULL;

/**
 * Free all minimal bases stored in the graph without touching its topology.
 * Mirrors the reset_bases() helper in main_cli.c.
 */
static void reset_bases(void) {
    if (!g) return;

    /* Libérer les bases minimales */
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
        g->minimals_basis   = NULL;
        g->nb_minimal_bases = 0;
        g->basis_dimension  = 0;
    }

    /*
     * Libérer aussi les données de faces et remettre nb_faces à 0.
     *
     * horton() ne rappelle find_faces() que si nb_faces == 0.
     * Si on ne remet pas nb_faces à 0 ici, les faces restent celles de
     * l'ancien graphe après une mutation (add/delete/split), ce qui fausse
     * la détection "cette base est-elle la base des faces ?".
     *
     * On libère les tableaux internes des Path (edges_ids, etc.) mais PAS
     * le tableau g->faces lui-même — il est alloué une seule fois dans
     * create_graph / create_edge et réutilisé.
     */
    for (int i = 0; i < g->nb_faces; i++) {
        free(g->faces[i].edges_ids);
        free(g->faces[i].edges_labels);
        free(g->faces[i].vertices_ids);
        g->faces[i].edges_ids    = NULL;
        g->faces[i].edges_labels = NULL;
        g->faces[i].vertices_ids = NULL;
    }
    g->nb_faces         = 0;
    g->face_basis       = -1;
    free(g->face_basis_outer);
    g->face_basis_outer    = NULL;
    g->nb_face_basis_outer = 0;
}

/**
 * Run multiple_horton with a freshly-zeroed inversion table.
 * n <= 0  →  no-op.
 */
static void run_horton(int n) {
    if (!g || n <= 0 || g->nb_edges < 1) return;
    int *inv = calloc(g->nb_edges, sizeof(int));
    if (!inv) return;
    multiple_horton(g, inv, n);
    free(inv);
}

/**
 * Must be called once from JS after the WASM module is ready.
 * Creates the initial (empty) graph and seeds the PRNG.
 */
EMSCRIPTEN_KEEPALIVE
void wasm_init() {
    if (g) delete_graph(g);
    g = create_graph();
}

/**
 * Generate a random planar graph, then run Horton n_horton times.
 */
EMSCRIPTEN_KEEPALIVE
void wasm_generate_planar(int nb_v, int nb_e, int n_horton) {
    srand((unsigned)time(NULL));
    if (g) delete_graph(g);
    g = create_planar_graph(nb_v, nb_e);
    run_horton(n_horton);
}

/**
 * Generate a random outer-planar graph, then run Horton n_horton times.
 */
EMSCRIPTEN_KEEPALIVE
void wasm_generate_outer_planar(int nb_v, int nb_e, int n_horton) {
    srand((unsigned)time(NULL));
    if (g) delete_graph(g);
    g = create_outer_planar_graph(nb_v, nb_e);
    run_horton(n_horton);
}

/**
 * @brief Write the graph text to the MEMFS temp file.
 *
 * IMPORTANT: do NOT call this via ccall() with a 'string' argument for
 * large files.  ccall() uses allocateUTF8OnStack() which is limited to
 * the WASM stack size (typically 64 KB).  A file with 1000 bases can be
 * several megabytes — this overflows the stack and causes "index out of
 * bounds" errors.
 *
 * The correct approach is to write the file from JS using Module.FS.writeFile()
 * and then call wasm_load_from_file() with no arguments.
 * wasm_load_from_text() is kept only for small graphs / backward compatibility.
 *
 * @param text  NUL-terminated string in the save_graph file format.
 */
EMSCRIPTEN_KEEPALIVE
void wasm_load_from_text(const char *text) {
    FILE *f = fopen("/tmp/wasm_graph.txt", "w");
    if (!f) return;
    fputs(text, f);
    fclose(f);
    if (g) delete_graph(g);
    g = create_graph();
    load_graph(g, "/tmp/wasm_graph.txt");
}

/**
 * @brief Load a graph from the MEMFS temp file /tmp/wasm_graph.txt.
 *
 * The file must have been written beforehand using Module.FS.writeFile() from
 * JavaScript.  This avoids the WASM stack overflow that occurs when ccall()
 * tries to copy a large string via allocateUTF8OnStack().
 *
 * JS usage:
 *   Module.FS.writeFile('/tmp/wasm_graph.txt', text);
 *   Module._wasm_load_from_file();
 */
EMSCRIPTEN_KEEPALIVE
void wasm_load_from_file(void) {
    if (g) delete_graph(g);
    g = create_graph();
    load_graph(g, "/tmp/wasm_graph.txt");
}

/**
 * Run Horton's algorithm n times on the current graph.
 */
EMSCRIPTEN_KEEPALIVE
void wasm_run_horton(int n) {
    run_horton(n);
}

/**
 * Add a vertex at (x, y) and re-run Horton n_horton times.
 */
EMSCRIPTEN_KEEPALIVE
void wasm_add_vertex(double x, double y, int n_horton) {
    if (!g) return;
    reset_bases();
    create_vertex(g, x, y);
    run_horton(n_horton);
}

/**
 * Delete vertex vid (and all its incident edges), compact, then re-run Horton.
 */
EMSCRIPTEN_KEEPALIVE
void wasm_delete_vertex(int vid, int n_horton) {
    if (!g || vid < 0 || vid >= g->nb_vertex) return;
    /*
     * reset_bases() est indispensable avant toute mutation :
     *   1. Il remet nb_faces à 0, ce qui force find_faces() dans horton()
     *      à recalculer les faces du nouveau graphe.  Sans ça, mb_is_faces()
     *      compare contre les faces de l'ancien graphe et ne reconnaît jamais
     *      la base des faces.
     *   2. Il invalide les anciennes bases dont les edges_ids sont dimensionnés
     *      pour l'ancien nb_edges.  Après compactage, nb_edges diminue et
     *      already_exists_basis() lirait hors bornes → faux doublons.
     */
    reset_bases();
    delete_vertex(g, vid);   /* handles incident edges + compact_graph */
    run_horton(n_horton);
}

/**
 * Add edge (u, v) and re-run Horton n_horton times.
 */
EMSCRIPTEN_KEEPALIVE
void wasm_add_edge(int u, int v, int n_horton) {
    if (!g) return;
    /*
     * Ajouter une arête augmente nb_edges. Les edges_ids des anciennes bases
     * ont taille old_nb_edges. already_exists_basis() compare en utilisant
     * g->nb_edges (nouveau) → lecture hors bornes sur les anciennes bases.
     * On les invalide avant.
     */
    reset_bases();
    create_edge(g, u, v);
    run_horton(n_horton);
}

/**
 * Delete edge eid and re-run Horton n_horton times.
 */
EMSCRIPTEN_KEEPALIVE
void wasm_delete_edge(int eid, int n_horton) {
    if (!g || eid < 0 || eid >= g->nb_edges) return;
    /* Même raison que wasm_delete_vertex : reset_bases() force find_faces()
     * et invalide les anciennes bases dont les tailles sont devenues fausses
     * après compactage. */
    reset_bases();
    delete_edge(g, eid);   /* includes compact_graph */
    run_horton(n_horton);
}

/**
 * Move vertex vid to (x, y) and refresh the neighbour angle sort.
 * Does NOT re-run Horton (position changes don't affect cycle bases).
 */
EMSCRIPTEN_KEEPALIVE
void wasm_move_vertex(int vid, double x, double y) {
    if (!g || vid < 0 || vid >= g->nb_vertex) return;

    move_vertex(g, vid, x, y);   /* updates angles, neighbour sort, find_faces */

    /* Reset face-basis tracking */
    g->face_basis = -1;
    free(g->face_basis_outer);
    g->face_basis_outer    = NULL;
    g->nb_face_basis_outer = 0;

    /* Recompute is_faces / is_faces_outer for every stored basis against the
     * new face decomposition — the old flags are stale after a move. */
    for (int b = 0; b < g->nb_minimal_bases; b++) {
        const int result = mb_is_faces(g, &g->minimals_basis[b]);
        g->minimals_basis[b].is_faces       = (result == 1) ? 1 : 0;
        g->minimals_basis[b].is_faces_outer = (result == 2) ? 1 : 0;

        if (result == 1) {
            g->face_basis = b;
        } else if (result == 2) {
            int *tmp = realloc(g->face_basis_outer,
                               (g->nb_face_basis_outer + 1) * sizeof(int));
            if (tmp) {
                g->face_basis_outer = tmp;
                g->face_basis_outer[g->nb_face_basis_outer] = b;
                g->nb_face_basis_outer++;
            }
        }
    }
}

/**
 * Split a set of edges (each divided into k+1 segments) then re-run Horton.
 *
 * @param eids      Pointer to an int array allocated in WASM linear memory
 *                  by the JS caller (via Module._malloc + Module.HEAP32).
 * @param n_eids    Length of the eids array.
 * @param k         Number of vertices to insert per edge.
 * @param n_horton  Horton iterations to run afterwards.
 */
EMSCRIPTEN_KEEPALIVE
void wasm_split_edges(int *eids, int n_eids, int k, int n_horton) {
    if (!g) return;
    /*
     * split_edge() supprime une arête et en crée k+1 nouvelles → nb_edges
     * augmente. Les anciennes bases ont des edges_ids[old_nb_edges] et
     * référencent des IDs d'arêtes supprimées. already_exists_basis() lirait
     * hors bornes, et readGraph() JS trouverait des IDs introuvables dans
     * d.edges → crash reconstructCyclePath.
     */
    reset_bases();
    for (int i = 0; i < n_eids; i++)
        split_edge(g, eids[i], k);
    run_horton(n_horton);
}

/** Total number of vertex slots (including deleted vertices). */
EMSCRIPTEN_KEEPALIVE int wasm_nb_vertices(void) { return g ? g->nb_vertex : 0;  }

/** Total number of edge slots (including deleted edges). */
EMSCRIPTEN_KEEPALIVE int wasm_nb_edges(void) { return g ? g->nb_edges : 0;  }

/** Number of distinct minimal bases found so far. */
EMSCRIPTEN_KEEPALIVE int wasm_nb_bases(void) { return g ? g->nb_minimal_bases : 0;  }

/** Dimension of the cycle space (= number of cycles in each basis). */
EMSCRIPTEN_KEEPALIVE int wasm_basis_dimension(void) { return g ? g->basis_dimension : 0;  }

/**
 * Index of the face basis (0-based) or -1 if none found yet.
 */
EMSCRIPTEN_KEEPALIVE
int wasm_face_basis(void) {
    if (!g || g->face_basis < 0) return -1;
    return g->face_basis;
}

/**
 * @brief Number of bases that equal the outer face + D-1 interior faces.
 *
 * Returns 0 if none have been found yet.
 */
EMSCRIPTEN_KEEPALIVE
int wasm_nb_face_basis_outer(void) {
    if (!g) return 0;
    return g->nb_face_basis_outer;
}

/**
 * @brief Return the i-th (0-based) entry in the face_basis_outer array.
 *
 * Each entry is a 1-based basis index.  Returns -1 if @p i is out of range.
 */
EMSCRIPTEN_KEEPALIVE
int wasm_face_basis_outer_at(int i) {
    if (!g || i < 0 || i >= g->nb_face_basis_outer) return -1;
    return g->face_basis_outer[i];
}

/** 1 if the vertex has been deleted, 0 otherwise. */
EMSCRIPTEN_KEEPALIVE int wasm_vertex_deleted(int id) { return g->vertices[id].deleted; }

/** X coordinate of vertex id. */
EMSCRIPTEN_KEEPALIVE double wasm_vertex_x(int id) { return g->vertices[id].x; }

/** Y coordinate of vertex id. */
EMSCRIPTEN_KEEPALIVE double wasm_vertex_y(int id) { return g->vertices[id].y; }

/** 1 if the edge has been deleted, 0 otherwise. */
EMSCRIPTEN_KEEPALIVE int wasm_edge_deleted(int id) { return g->edges[id].deleted; }

/** ID of the first endpoint of edge id. */
EMSCRIPTEN_KEEPALIVE int wasm_edge_u(int id) { return g->edges[id].u; }

/** ID of the second endpoint of edge id. */
EMSCRIPTEN_KEEPALIVE int wasm_edge_v(int id) { return g->edges[id].v; }

/** Label of edge id (used for cycle basis identification). */
EMSCRIPTEN_KEEPALIVE int wasm_edge_label(int id) { return g->edges[id].label; }

/**
 * Return a raw pointer to the edges_ids binary vector of cycle c in basis b.
 *
 * The vector has g->nb_edges uint32_t entries; entry i is 1 if the edge
 * whose **id** equals i belongs to the cycle, 0 otherwise.
 *
 * JS reads it as:
 *   const ptr = Module._wasm_cycle_edges_ids_ptr(b, c);
 *   for (let i = 0; i < nbEdgesTotal; i++) {
 *       if (Module.HEAPU32[ptr/4 + i]) edgeIds.push(i);
 *   }
 *
 * The pointer is valid as long as the graph is not mutated.
 */
/**
 * @brief Read one entry of a cycle's edges_ids vector.
 *
 * Returns 1 if the edge with index @p edge_idx belongs to cycle @p c of basis
 * @p b, 0 otherwise (or if any index is out of range).
 *
 * This accessor is the safe alternative to reading Module.HEAPU32 directly
 * from JavaScript.  Direct HEAPU32 reads become stale when WASM memory grows
 * (ALLOW_MEMORY_GROWTH=1 replaces the underlying ArrayBuffer), causing
 * out-of-bounds errors.  Routing through this function always reads from the
 * current WASM linear memory.
 *
 * @param b        Basis index (0-based).
 * @param c        Cycle index within the basis (0-based).
 * @param edge_idx Edge slot index (0-based, up to nb_edges - 1).
 * @return 1 if the edge is part of the cycle, 0 otherwise.
 */
EMSCRIPTEN_KEEPALIVE
int wasm_cycle_edge_at(int b, int c, int edge_idx) {
    if (!g || b < 0 || b >= g->nb_minimal_bases ||
        c < 0 || c >= g->basis_dimension ||
        edge_idx < 0 || edge_idx >= g->nb_edges) return 0;
    return (int)g->minimals_basis[b].cycles[c].edges_ids[edge_idx];
}

/**
 * @brief Write one int32 value into a caller-allocated WASM buffer.
 *
 * This accessor lets JavaScript pass integer arrays to C functions without
 * touching Module.HEAP32 directly.  Direct HEAP32 writes are unsafe when WASM
 * memory grows between _malloc() and the write, because the ArrayBuffer is
 * replaced and the JS-side HEAP32 view becomes stale.
 *
 * Usage from JS:
 * @code
 * const ptr = M._malloc(n * 4);
 * for (let i = 0; i < n; i++) M._wasm_write_int(ptr, i, values[i]);
 * M._wasm_some_function(ptr, n);
 * M._free(ptr);
 * @endcode
 *
 * @param base_ptr  Byte address of the beginning of the buffer (from _malloc).
 * @param index     0-based element index.
 * @param value     Integer value to write.
 */
EMSCRIPTEN_KEEPALIVE
void wasm_write_int(int *base_ptr, int index, int value) {
    base_ptr[index] = value;
}