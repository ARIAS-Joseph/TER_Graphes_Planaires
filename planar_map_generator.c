#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"
#include "dfs_graph.h"

static void dfs_biconnected(const Graph *g, int u, int parent, int *time, int *disc, int *low,
    int *visited, int *has_articulation)
{
    visited[u] = 1;
    disc[u] = low[u] = ++(*time);

    int children = 0;

    Neighbor_list neigh = g->neighbors[u];

    for (int i = 0; i < neigh.count; i++) {
        int v = neigh.neighbors[i];
        if (g->vertices[v].deleted) continue;

        if (!visited[v]) {
            children++;

            dfs_biconnected(g, v, u, time,
                            disc, low, visited, has_articulation);

            if (low[v] < low[u])
                low[u] = low[v];

            if (parent != -1 && low[v] >= disc[u])
                *has_articulation = 1;
        }
        else if (v != parent) {
            if (disc[v] < low[u])
                low[u] = disc[v];
        }
    }

    if (parent == -1 && children > 1)
        *has_articulation = 1;
}


int is_biconnected(const Graph *g)
{
    if (g->nb_vertices < 2)
        return 0;

    int n = g->nb_vertices;

    int *disc = calloc(n, sizeof(int));
    int *low = calloc(n, sizeof(int));
    int *visited = calloc(n, sizeof(int));

    int time = 0;
    int has_articulation = 0;

    dfs_biconnected(g, 0, -1,
                    &time,
                    disc,
                    low,
                    visited,
                    &has_articulation);

    /* graphe non connexe ? */
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            free(disc);
            free(low);
            free(visited);
            return 0;
        }
    }

    free(disc);
    free(low);
    free(visited);

    return !has_articulation;
}

int is_ancestor(const DfsGraph *dfs, const int u, const int v) {
    return dfs->vertices[u].ancestor_start <= dfs->vertices[v].ancestor_start &&
        dfs->vertices[v].ancestor_start <= dfs->vertices[u].ancestor_end;
}

int is_descendant_edge(const DfsGraph *dfs, const int e1, const int e2) {
    if (e1 == e2) {
        return 1;
    }
    const DfsEdge *a = &dfs->edges[e1];
    if (a->type == BACK_EDGE) {
        return 0;
    }
    const int src = dfs->edges[e2].from;
    return dfs->vertices[src].ancestor_start >= a->subtree_start
    && dfs->vertices[src].ancestor_start <= a->subtree_end;
}

static void add_dfs_edge(DfsGraph *dfs, const int from, const int to, const int id, const DfsEdgeType type) {

    DfsEdge *e = &dfs->edges[id];

    e->id = id;
    e->from = from;
    e->to = to;
    e->type = type;

    e->parent_edge = -1;

    e->successors = NULL;
    e->nb_successors = 0;

    e->subtree_start = -1;
    e->subtree_end = -1;

    e->label = -1;

    e->singular = 0;
    e->is_reference = 0;
    e->phi = -1;
    e->low1 = -1;
    e->low2 = -1;

    dfs->edge_count++;
}

static int dfs_visit_call_count = 0;

static void dfs_visit(Graph *graph, DfsGraph *dfs, const int u, int *visited, int *time) {

    if (visited == NULL || time == NULL) {
        printf( "dfs_visit: visited or time is NULL\n");
        exit(EXIT_FAILURE);
    }

    visited[u] = 1;

    /* À ajouter au tout début de dfs_visit, dans planar_map_generator.c,
 * juste après `visited[u] = 1;` */
    if (dfs_visit_call_count++ == 0) {
        /* Premier appel = visite de la racine : dump complet ici */
        printf("=== Neighbors AS SEEN BY dfs_visit (live) ===\n");
        for (int vv = 0; vv < graph->nb_vertices; vv++) {
            printf("%d:", vv);
            for (int i = 0; i < graph->neighbors[vv].count; i++)
                printf(" %d", graph->neighbors[vv].neighbors[i]);
            printf("\n");
        }
    }

    DfsVertex *vu = &dfs->vertices[u];

    vu->dfs_num = (*time)++;
    vu->ancestor_start = vu->dfs_num;

    const Neighbor_list *neighbors = &graph->neighbors[u];

    int *out_edges = malloc(neighbors->count * sizeof(int));
    int nb_out_edges = 0;

    for (int i = 0; i < neighbors->count; i++) {

        const int v = neighbors->neighbors[i];
        if (graph->vertices[v].deleted) continue;

        if (!visited[v]) {

            vu->children[vu->nb_children++] = v;
            dfs->vertices[v].parent = u;

            const int edge_id = graph->edge_indices[u * graph->nb_vertices + v];
            add_dfs_edge(dfs, u, v, edge_id, TREE_EDGE);
            dfs->vertices[v].parent_edge = edge_id;
            /* DEBUG */
            printf( "  dfs_visit: set parent_edge[%d] = arc %d (%d -> %d)\n",
                    v, edge_id,
                    dfs->edges[edge_id].from,
                    dfs->edges[edge_id].to);
            out_edges[nb_out_edges++] = edge_id;

            dfs_visit(graph, dfs, v, visited, time);

        } else if (v != vu->parent && dfs->vertices[v].dfs_num < vu->dfs_num) {

            const int edge_id = graph->edge_indices[u * graph->nb_vertices + v];
            add_dfs_edge(dfs, u, v, edge_id, BACK_EDGE);
            out_edges[nb_out_edges++] = edge_id;
        }
    }

    vu->ancestor_end = *time - 1;

    if (vu->parent_edge != -1) {
        DfsEdge *pe = &dfs->edges[vu->parent_edge];
        pe->successors = out_edges;
        pe->nb_successors = nb_out_edges;
        pe->subtree_start = vu->ancestor_start;
        pe->subtree_end = vu->ancestor_end;
    } else {
        free(out_edges);
    }
}

DfsGraph *build_dfs_graph(Graph *graph) {

    dfs_visit_call_count = 0;

    prepare_graph_matrices(graph);

    DfsGraph *dfs = calloc(1, sizeof(DfsGraph));
    if (!dfs) exit(EXIT_FAILURE);

    /* CANARI DEBUG : si edge_count != 0 ici, le calloc n'a pas donné une
     * zone réellement vierge -- preuve directe de heap corruption en amont. */
    if (dfs->edge_count != 0 || dfs->vertices_count != 0) {
        fprintf(stderr,
            "*** HEAP CORRUPTION DETECTED: fresh calloc'd DfsGraph has "
            "edge_count=%d vertices_count=%d (expected 0,0) ***\n",
            dfs->edge_count, dfs->vertices_count);
        abort();   /* force un crash immédiat et localisé ici, plutôt que de
                      laisser la corruption se propager silencieusement plus loin */
    }

    const int n = graph->nb_vertices;

    dfs->vertices_count = n;
    dfs->edge_capacity = 2 * graph->nb_edges;
    dfs->edge_indices = graph->edge_indices;

    dfs->root = 0;

    dfs->vertices = calloc(n, sizeof(DfsVertex));

    if (!dfs->vertices) exit(EXIT_FAILURE);

    for (int i = 0; i < n; i++) {
        dfs->vertices[i].id = graph->vertices[i].id;
        dfs->vertices[i].parent = -1;
        dfs->vertices[i].parent_edge = -1;
        dfs->vertices[i].dfs_num = -1;
        dfs->vertices[i].ref_edge = -1;

        dfs->vertices[i].children = malloc(n * sizeof(int));
        dfs->vertices[i].nb_children = 0;
    }

    dfs->edges = malloc(dfs->edge_capacity * sizeof(DfsEdge));
    dfs->edge_count = 0;

    int *visited = calloc(n, sizeof(int));
    int time = 1;

    dfs_visit(graph, dfs, dfs->root, visited, &time);

    free(visited);

    return dfs;
}

void free_dfs_graph(DfsGraph *dfs) {
    if (!dfs) {
        return;
    }

    for (int v = 0; v < dfs->vertices_count; ++v) {

        free(dfs->vertices[v].children);
        free(dfs->vertices[v].out_edges);
    }

    for (int e = 0; e < dfs->edge_count; ++e) {
        free(dfs->edges[e].successors);
    }

    for (int i = 0; i < dfs->nb_singular_sets; ++i) {
        free(dfs->singular_sets[i].edges);
    }

    free(dfs->singular_sets);

    free(dfs->E0);

    free(dfs->same_parent);
    free(dfs->same_partner);
    free(dfs->edge_to_group);

    if (dfs->SAME.groups) {
        for (int i = 0; i < dfs->SAME.nb_groups; ++i) {
            free(dfs->SAME.groups[i].edges);
        }
        free(dfs->SAME.groups);
    }

    free(dfs->DIFF.partner);

    if (dfs->SAME_prime.groups) {
        for (int i = 0; i < dfs->SAME_prime.nb_groups; ++i) {
            free(dfs->SAME_prime.groups[i].edges);
        }
        free(dfs->SAME_prime.groups);
    }

    free(dfs->DIFF_prime.partner);

    /* FIX (bug #4) : nouveau champ, oublie de liberation sinon. */
    free(dfs->same_to_prime);

    free(dfs->vertices);
    free(dfs->edges);

    free(dfs);
}

static void merge_low(int *best1, int *best2, const int val) {
    if (val == *best1 || val == *best2) {
        return;
    }
    if (val < *best1) {
        *best2 = *best1;
        *best1 = val;
    } else if (val < *best2) {
        *best2 = val;
    }
}

void y_e(const DfsGraph *dfs, const int e) {

    DfsEdge *edge = &dfs->edges[e];
    const int max = dfs->vertices_count + 1;

    if (edge->type == BACK_EDGE) {
        edge->low1 = dfs->vertices[edge->to].dfs_num;
        edge->low2 = max;
        return;
    }

    int best1 = max;
    int best2 = max;

    for (int i = 0; i < edge->nb_successors; i++) {
        const DfsEdge *succ = &dfs->edges[edge->successors[i]];
        merge_low(&best1, &best2, succ->low1);
        merge_low(&best1, &best2, succ->low2);
    }

    edge->low1 = best1;
    edge->low2 = best2;
}

static void print_low_values(const DfsGraph *dfs)
{
    for(int e=0;e<dfs->edge_count;e++)
    {
        printf("v=%2d low1=%2d low2=%2d\n",
            e,
            dfs->edges[e].low1,
            dfs->edges[e].low2);
    }
}

void compute_low_values(const DfsGraph *dfs)  {

    for (int i = 0; i < dfs->edge_count; i++) {
        if (dfs->edges[i].type == BACK_EDGE) {
            y_e(dfs, i);
        }
    }

    const int n = dfs->vertices_count;
    int *order = malloc(n * sizeof(int));
    for (int v = 0; v < n; v++) {
        order[dfs->vertices[v].dfs_num-1] = v;
    }

    for (int rank = n - 1; rank >= 1; rank--) {
        const int v = order[rank];
        y_e(dfs, dfs->vertices[v].parent_edge);
    }

    free(order);

    printf("\n================ LOW ================\n");
    print_low_values(dfs);
}

static int phi_of(const DfsGraph *dfs, const int e) {
    DfsEdge *edge = &dfs->edges[e];
    const int v = dfs->vertices[edge->from].dfs_num;  /* "v" au sens DFS */

    /* FIX (bug #5) : on fixe explicitement les DEUX branches, plutot que de
       ne mettre `singular` a 1 que dans un seul cas et de laisser l'autre
       cas a sa valeur non-initialisee (malloc, pas calloc, sur dfs->edges). */
    edge->singular = edge->low2 >= v ? 1 : 0;

    return edge->singular ? 2 * edge->low1 : 2 * edge->low1 + 1;
}

static void print_phi_lists(const DfsGraph *dfs)
{
    printf("\nPhi(v)\n");

    for(int v=0;v<dfs->vertices_count;v++){

        printf("%2d : ",v);

        for(int i=0;i<dfs->vertices[v].nb_out_edges;i++)
            printf("%d ",dfs->vertices[v].out_edges[i]);

        printf("\n");
    }
}

void build_phi_lists(DfsGraph *dfs) {

    const int n = dfs->vertices_count;
    const int m = dfs->edge_count;

    int max_phi = 0;
    for (int i = 0; i < m; i++) {
        dfs->edges[i].phi = phi_of(dfs, i);
        if (dfs->edges[i].phi > max_phi) max_phi = dfs->edges[i].phi;
    }

    /* 2. tri par seau (chainage), O(m + max_phi) */
    int *bucket_head = malloc((max_phi + 1) * sizeof(int));
    int *bucket_next = malloc(m * sizeof(int));
    for (int p = 0; p <= max_phi; p++) bucket_head[p] = -1;

    for (int i = 0; i < m; i++) {
        const int p = dfs->edges[i].phi;
        bucket_next[i] = bucket_head[p];
        bucket_head[p] = i;
        /* l'ordre entre deux arcs de meme phi n'a pas besoin d'etre
           fixe : par le Lemme 3 de l'article, deux arcs singuliers de
           meme low1 peuvent toujours etre echanges dans M'(v). */
    }

    /* 3. on compte le degre sortant de chaque sommet pour allouer
       chaque Phi(v) une seule fois */
    int *count = calloc(n, sizeof(int));
    for (int i = 0; i < m; i++) count[dfs->edges[i].from]++;
    for (int v = 0; v < n; v++) {
        dfs->vertices[v].out_edges = malloc(count[v] * sizeof(int));
        dfs->vertices[v].nb_out_edges = 0;
    }
    free(count);

    /* 4. on parcourt les seaux par phi croissant : chaque Phi(v) se
       remplit alors automatiquement dans le bon ordre */
    for (int p = 0; p <= max_phi; p++) {
        for (int i = bucket_head[p]; i != -1; i = bucket_next[i]) {
            DfsVertex *vv = &dfs->vertices[dfs->edges[i].from];
            vv->out_edges[vv->nb_out_edges++] = i;
        }
    }
    free(bucket_head);
    free(bucket_next);

    /* 5. e_{v,ref} = premier element de Phi(v) ; E0 = tout le reste */
    for (int i = 0; i < m; i++) dfs->edges[i].is_reference = 0;

    for (int v = 0; v < n; v++) {
        DfsVertex *vv = &dfs->vertices[v];

        /* FIX (bug #6) : garde absente jusqu'ici. Sur un graphe vraiment
           2-connexe, ce cas ne devrait jamais se produire (cf. argument
           dans la revue), mais sans cette garde, out_edges[0] sur un
           tableau alloue a 0 element est un acces hors-limites silencieux. */
        if (vv->nb_out_edges == 0) {
            printf(
                    "build_phi_lists: sommet (dfs_num=%d) sans arc sortant "
                    "-- le graphe n'est peut-etre pas 2-connexe.\n",
                    vv->dfs_num);
            continue;
        }

        vv->ref_edge = vv->out_edges[0];
        dfs->edges[vv->out_edges[0]].is_reference = 1;
    }

    int e0_count = 0;
    for (int i = 0; i < m; i++) if (!dfs->edges[i].is_reference) e0_count++;

    dfs->E0 = malloc(e0_count * sizeof(int));
    dfs->E0_size = 0;
    for (int i = 0; i < m; i++) {
        if (!dfs->edges[i].is_reference) {
            dfs->E0[dfs->E0_size++] = i;
        }
    }

    printf("\n================ PHI ================\n");
    print_phi_lists(dfs);

}

static void print_singular_sets(const DfsGraph *dfs)
{
    for(int s=0;s<dfs->nb_singular_sets;s++)
    {
        printf("S%d : ",s);

        SingularSet *set=&dfs->singular_sets[s];

        for(int i=0;i<set->size;i++)
        {
            int e=set->edges[i];

            printf("%d(%d-%d) ",
                e,
                dfs->edges[e].from,
                dfs->edges[e].to);
        }

        printf("\n");
    }
}

void build_singular_sets(DfsGraph *dfs) {

    dfs->nb_singular_sets = 0;

    int capacity = 4;

    dfs->singular_sets = malloc(capacity * sizeof(SingularSet));

    if (!dfs->singular_sets) {
        perror("malloc singular_sets");
        exit(EXIT_FAILURE);
    }

    for (int v = 0; v < dfs->vertices_count; ++v) {
        const DfsVertex *vertex = &dfs->vertices[v];
        int current_set = -1;

        if (vertex->nb_out_edges == 0) {
            continue;
        }

        int current_low1 = -1;

        for (int i = 0; i < vertex->nb_out_edges; ++i) {

            const int edge_id = vertex->out_edges[i];

            if (!(dfs->edges[edge_id].singular && !dfs->edges[edge_id].is_reference)) {
                continue;
            }

            const DfsEdge *e = &dfs->edges[edge_id];

            if (e->low1 != current_low1) {

                current_low1 = e->low1;

                if (dfs->nb_singular_sets >= capacity) {

                    capacity *= 2;

                    SingularSet *tmp = realloc(
                        dfs->singular_sets,
                        capacity * sizeof(SingularSet)
                    );

                    if (!tmp) {
                        perror("realloc singular_sets");
                        exit(EXIT_FAILURE);
                    }

                    dfs->singular_sets = tmp;
                }

                current_set = dfs->nb_singular_sets++;

                SingularSet *set = &dfs->singular_sets[current_set];

                set->vertex = v;
                set->low1 = current_low1;

                set->edges = malloc(
                    vertex->nb_out_edges * sizeof(int)
                );

                if (!set->edges) {
                    perror("malloc set->edges");
                    exit(EXIT_FAILURE);
                }

                set->size = 0;
            }

            SingularSet *set = &dfs->singular_sets[current_set];

            set->edges[set->size++] = edge_id;
        }
    }

    printf("\n================ SINGULAR SETS ================\n");
    print_singular_sets(dfs);

}

static int uf_find(int *parent, int x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

/* Union SAME pure (sans toucher au partner) */
static void uf_union(int *parent, const int x, const int y) {
    const int rx = uf_find(parent, x);
    const int ry = uf_find(parent, y);
    if (rx != ry) parent[rx] = ry;
}

/* partner(g) -- renvoie -1 si g est "libre" (pas de contrainte DIFF) */
static int get_partner(int *parent, const int *partner, int g) {
    if (g == -1) return -1;
    const int r = uf_find(parent, g);
    if (partner[r] == -1) return -1;
    return uf_find(parent, partner[r]);
}

/*
 * Remplace les deux paires DIFF (a1,b1) et (a2,b2)
 * par la paire unique (a1 U a2, b1 U b2).
 * -1 represente l'ensemble vide.
 */
static void merge_pairs(int *parent, int *partner, const int a1, const int b1, const int a2,
    const int b2) {
    int newA, newB;

    if (a1 == -1) newA = (a2 == -1) ? -1 : uf_find(parent, a2);
    else if (a2 == -1) newA = uf_find(parent, a1);
    else { uf_union(parent, a1, a2); newA = uf_find(parent, a1); }

    if (b1 == -1) newB = (b2 == -1) ? -1 : uf_find(parent, b2);
    else if (b2 == -1) newB = uf_find(parent, b1);
    else { uf_union(parent, b1, b2); newB = uf_find(parent, b1); }

    if (newA != -1) partner[newA] = newB;
    if (newB != -1) partner[newB] = newA;
}

/*
 * Cas "Bi entrelace X uniquement" (Etape 4, fusion en croix) :
 * remplace (Bi,U) et (X,Y) par (U U X, Bi U Y).
 */
static void merge_pairs_cross(int *parent, int *partner, const int Bi, const int U, const int X,
    const int Y) {
    int newFirst, newSecond;

    if (U == -1) newFirst = (X == -1) ? -1 : uf_find(parent, X);
    else if (X == -1) newFirst = uf_find(parent, U);
    else { uf_union(parent, U, X); newFirst = uf_find(parent, U); }

    if (Bi == -1) newSecond = (Y == -1) ? -1 : uf_find(parent, Y);
    else if (Y == -1) newSecond = uf_find(parent, Bi);
    else { uf_union(parent, Bi, Y); newSecond = uf_find(parent, Bi); }

    if (newFirst != -1) partner[newFirst] = newSecond;
    if (newSecond != -1) partner[newSecond] = newFirst;
}

static void block_init(Block *b) {
    b->size = 0;
    b->capacity = 4;
    b->attachments = malloc(b->capacity * sizeof(int));
    b->representative = -1;
}

static void block_push(Block *b, const int val) {
    if (b->size >= b->capacity) {
        b->capacity *= 2;
        int* tpm = realloc(b->attachments, b->capacity * sizeof(int));
        if (!tpm) {
            perror("block_push: realloc");
            exit(EXIT_FAILURE);
        }
        b->attachments = tpm;
    }
    b->attachments[b->size++] = val;
}

static void block_remove_value(Block *b, const int val) {
    int w = 0;
    for (int r = 0; r < b->size; r++) {
        if (b->attachments[r] != val) b->attachments[w++] = b->attachments[r];
    }
    b->size = w;
}

static void block_concat(Block *dst, const Block *src) {
    for (int i = 0; i < src->size; i++) block_push(dst, src->attachments[i]);
}

static void block_free(Block *b) {
    free(b->attachments);
    b->attachments = NULL;
    b->size = b->capacity = 0;
}

/* first()/last() : conventions du papier (n+1 / 0 si vide) */
static int block_last(const Block *b, const int n) {
    (void)n;
    if (b->size == 0) return 0;
    return b->attachments[b->size - 1];
}

/* X est-il un bloc "normal" vis-a-vis de l'arc dont low1=low1_e, sommet=w_num ? */
static int block_is_normal(const Block *b, const int low1_e, const int w_num) {
    for (int i = 0; i < b->size; i++) {
        int y = b->attachments[i];
        if (low1_e < y && y < w_num) return 1;
    }
    return 0;
}

static void attlist_init(AttList *l) {
    l->size = 0;
    l->capacity = 4;
    l->pairs = malloc(l->capacity * sizeof(BlockPair));
}

static void attlist_push_back(AttList *l, const BlockPair p) {
    if (l->size >= l->capacity) {
        l->capacity *= 2;
        BlockPair* tpm = realloc(l->pairs, l->capacity * sizeof(BlockPair));
        if (!tpm) {
            perror("attlist_push_back: realloc");
            exit(EXIT_FAILURE);
        }
        l->pairs = tpm;
    }
    l->pairs[l->size++] = p;
}

static void attlist_free(AttList *l) {
    for (int i = 0; i < l->size; i++) {
        block_free(&l->pairs[i].left);
        block_free(&l->pairs[i].right);
    }
    free(l->pairs);
    l->pairs = NULL;
    l->size = l->capacity = 0;
}

/* Supprime les occurrences de `val` dans tous les blocs de l, puis
   retire les paires devenues entierement vides. */
static void attlist_remove_value(AttList *l, const int val) {
    int w = 0;
    for (int r = 0; r < l->size; r++) {
        block_remove_value(&l->pairs[r].left, val);
        block_remove_value(&l->pairs[r].right, val);
        if (l->pairs[r].left.size == 0 && l->pairs[r].right.size == 0) {
            block_free(&l->pairs[r].left);
            block_free(&l->pairs[r].right);
        } else {
            l->pairs[w++] = l->pairs[r];
        }
    }
    l->size = w;
}

static void build_back_edge_att(const DfsGraph *dfs, AttList *att, const int e) {
    const int y = dfs->vertices[dfs->edges[e].to].dfs_num;
    attlist_init(&att[e]);

    BlockPair p;
    block_init(&p.left);
    block_push(&p.left, y);
    block_init(&p.right);

    attlist_push_back(&att[e], p);
}

/* applique les regles SAME/DIFF de l'Etape 2 pour le successeur ei,
   met a jour att[ei] (consomme), et renvoie le bloc fusionne Bi */
static Block step2_build_Bi(const DfsGraph *dfs, AttList *att, const int ei, const int w_num) {
    const int low1_ei = dfs->edges[ei].low1;

    int Bi_root = uf_find(dfs->same_parent, ei); /* ei in E0 */

    AttList *li = &att[ei];
    for (int p = 0; p < li->size; p++) {
        const Block *Lx = &li->pairs[p].left;
        const Block *Ly = &li->pairs[p].right;

        /* [X,Y] avec X=left,Y=right */
        if (block_is_normal(Lx, low1_ei, w_num) && Lx->representative != -1) {
            const int X = uf_find(dfs->same_parent, Lx->representative);
            const int Y = (Ly->representative == -1) ? -1
                     : uf_find(dfs->same_parent, Ly->representative);
            const int U = get_partner(dfs->same_parent, dfs->same_partner, Bi_root);
            merge_pairs(dfs->same_parent, dfs->same_partner, Bi_root, U, X, Y);
            Bi_root = uf_find(dfs->same_parent, Bi_root);
        }
        /* [Y,X] avec X=right,Y=left */
        if (block_is_normal(Ly, low1_ei, w_num) && Ly->representative != -1) {
            const int X = uf_find(dfs->same_parent, Ly->representative);
            const int Y = (Lx->representative == -1) ? -1
                     : uf_find(dfs->same_parent, Lx->representative);
            const int U = get_partner(dfs->same_parent, dfs->same_partner, Bi_root);
            merge_pairs(dfs->same_parent, dfs->same_partner, Bi_root, U, X, Y);
            Bi_root = uf_find(dfs->same_parent, Bi_root);
        }
    }

    /* concatenation de tous les blocs de att(ei) en un seul bloc Bi */
    Block Bi;
    block_init(&Bi);
    Bi.representative = Bi_root;
    for (int p = 0; p < li->size; p++) {
        block_concat(&Bi, &li->pairs[p].left);
        block_concat(&Bi, &li->pairs[p].right);
    }

    attlist_free(li);
    return Bi;
}

/* Etape 3 : fusionne les blocs "hauts" (last > low1_e2) au sommet
   de att[e] en un seul bloc, en propageant SAME/DIFF. */
static void step3_merge_high_blocks(const DfsGraph *dfs, AttList *att, const int e,
    const int low1_e2, const int n) {

    Block mergedLeft, mergedRight;
    block_init(&mergedLeft);
    block_init(&mergedRight);
    mergedLeft.representative = -1;
    mergedRight.representative = -1;

    int has_merge = 0;

    while (att[e].size > 0) {
        BlockPair *top = &att[e].pairs[att[e].size - 1];
        if (block_last(&top->left, n) <= low1_e2 &&
            block_last(&top->right, n) <= low1_e2) {
            break; /* plus rien a fusionner */
        }

        if (has_merge) {
            const int X = mergedLeft.representative;
            const int Y = mergedRight.representative;
            const int X1 = top->left.representative;
            const int Y1 = top->right.representative;
            merge_pairs(dfs->same_parent, dfs->same_partner, X, Y, X1, Y1);
            mergedLeft.representative =
                (X != -1) ? uf_find(dfs->same_parent, X)
                          : ((X1 != -1) ? uf_find(dfs->same_parent, X1) : -1);
            mergedRight.representative =
                (Y != -1) ? uf_find(dfs->same_parent, Y)
                          : ((Y1 != -1) ? uf_find(dfs->same_parent, Y1) : -1);
        } else {
            mergedLeft.representative = top->left.representative;
            mergedRight.representative = top->right.representative;
            has_merge = 1;
        }

        block_concat(&mergedLeft, &top->left);
        block_concat(&mergedRight, &top->right);
        block_free(&top->left);
        block_free(&top->right);
        att[e].size--;
    }

    if (has_merge) {
        BlockPair newp;
        newp.left = mergedLeft;
        newp.right = mergedRight;
        attlist_push_back(&att[e], newp);
    } else {
        block_free(&mergedLeft);
        block_free(&mergedRight);
    }
}

/* Etape 4 : insere le bloc Bi dans att[e]. */
static void step4_insert_block(const DfsGraph *dfs, AttList *att, const int e, Block *Bi,
    const int low1_ei, const int n) {

    if (att[e].size == 0) {
        BlockPair p;
        block_init(&p.left);
        p.left.representative = Bi->representative;
        block_concat(&p.left, Bi);
        block_init(&p.right);
        attlist_push_back(&att[e], p);
        block_free(Bi);
        return;
    }

    BlockPair *top = &att[e].pairs[att[e].size - 1];
    const int interlaceX = (low1_ei < block_last(&top->left, n));
    const int interlaceY = (low1_ei < block_last(&top->right, n));

    const int U = get_partner(dfs->same_parent, dfs->same_partner,
                         Bi->representative);

    if (interlaceX && interlaceY) {
        /* non planaire : on ignore (suppose deja planaire en entree) */
        block_free(Bi);
        return;
    }

    if (interlaceX) {
        /* fusion Bi dans Y, en croix : (U,X) <- (Bi,Y) */
        block_concat(&top->right, Bi);
        merge_pairs_cross(dfs->same_parent, dfs->same_partner,
                           Bi->representative, U,
                           top->left.representative,
                           top->right.representative);
        if (block_last(&top->left, n) < block_last(&top->right, n)) {
            Block tmp = top->left;
            top->left = top->right;
            top->right = tmp;
        }
        block_free(Bi);
        return;
    }

    /* Bi n'entrelace ni X ni Y : nouvelle paire au sommet */
    BlockPair newp;
    newp.left = *Bi;
    block_init(&newp.right);
    newp.right.representative = U;
    attlist_push_back(&att[e], newp);
}

/* Construit att(e) pour un arc d'arbre e=[v,w] */
static void build_tree_edge_att(const DfsGraph *dfs, AttList *att, const int e) {

    if (att == NULL) return;
    const int w = dfs->edges[e].to;
    const int w_num = dfs->vertices[w].dfs_num;
    const int n = dfs->vertices_count;

    const DfsVertex *vw = &dfs->vertices[w];
    const int k = vw->nb_out_edges;
    if (k == 0) { attlist_init(&att[e]); return; }

    const int e1 = vw->out_edges[0];

    /* Etape 1 : transfert de propriete -- e1 ne doit plus posseder ce pointeur */
    att[e] = att[e1];
    att[e1].pairs = NULL;
    att[e1].size = 0;
    att[e1].capacity = 0;

    attlist_remove_value(&att[e], w_num);

    if (k == 1) return;

    /* Etape 2 */
    Block *Bs = malloc((k - 1) * sizeof(Block));
    for (int i = 1; i < k; i++) {
        const int ei = vw->out_edges[i];
        Bs[i - 1] = step2_build_Bi(dfs, att, ei, w_num);
    }

    /* Etape 3 (low1 du tout premier successeur non-reference, e2) */
    const int e2 = vw->out_edges[1];
    step3_merge_high_blocks(dfs, att, e, dfs->edges[e2].low1, n);

    /* Etape 4 */
    for (int i = 1; i < k; i++) {
        const int ei = vw->out_edges[i];
        step4_insert_block(dfs, att, e, &Bs[i - 1],
                            dfs->edges[ei].low1, n);
    }
    free(Bs);
}

/* parcours recursif simple (suffisant pour des graphes raisonnables) */
static void recurse(const int v, DfsGraph *dfs, AttList *att) {
    const DfsVertex *vv = &dfs->vertices[v];
    for (int i = 0; i < vv->nb_children; i++) {
        recurse(vv->children[i], dfs, att);
    }
    if (vv->parent_edge != -1) {
        build_tree_edge_att(dfs, att, vv->parent_edge);
    }
}

void build_SAME_DIFF(DfsGraph *dfs) {

    int m = dfs->edge_count;
    int *root_of = malloc(m * sizeof(int));
    int *count   = calloc(m, sizeof(int));

    for (int i = 0; i < dfs->E0_size; i++) {
        int e = dfs->E0[i];
        int r = uf_find(dfs->same_parent, e);
        root_of[e] = r;
        count[r]++;
    }

    int *group_index = malloc(m * sizeof(int));
    for (int i = 0; i < m; i++) group_index[i] = -1;

    int nb_groups = 0;
    for (int i = 0; i < dfs->E0_size; i++) {
        int e = dfs->E0[i];
        int r = root_of[e];
        if (group_index[r] == -1) group_index[r] = nb_groups++;
    }

    dfs->SAME.nb_groups = nb_groups;
    dfs->SAME.groups = malloc(nb_groups * sizeof(SameGroup));
    for (int i = 0; i < m; i++) {
        if (group_index[i] != -1) {
            int gi = group_index[i];
            dfs->SAME.groups[gi].edges = malloc(count[i] * sizeof(int));
            dfs->SAME.groups[gi].size  = 0;
        }
    }

    dfs->edge_to_group = malloc(m * sizeof(int));
    for (int i = 0; i < m; i++) dfs->edge_to_group[i] = -1;

    for (int i = 0; i < dfs->E0_size; i++) {
        int e  = dfs->E0[i];
        int gi = group_index[root_of[e]];
        SameGroup *g = &dfs->SAME.groups[gi];
        g->edges[g->size++] = e;
        dfs->edge_to_group[e] = gi;
    }

    dfs->DIFF.nb_groups = nb_groups;
    dfs->DIFF.partner   = malloc(nb_groups * sizeof(int));
    for (int i = 0; i < m; i++) {
        if (group_index[i] == -1) continue;
        const int gi = group_index[i];
        const int p_root = get_partner(dfs->same_parent, dfs->same_partner, i);
        dfs->DIFF.partner[gi] = (p_root == -1) ? -1 : group_index[p_root];
    }

    free(root_of); free(count); free(group_index);
}

static void print_same_sets(const DfsGraph *dfs)
{
    for(int g=0;g<dfs->SAME.nb_groups;g++)
    {
        printf("G%d : ",g);

        SameGroup *grp=&dfs->SAME.groups[g];

        for(int i=0;i<grp->size;i++)
        {
            int e=grp->edges[i];

            printf("%d(%d-%d) ",
                e,
                dfs->edges[e].from,
                dfs->edges[e].to);
        }

        printf("\n");
    }
}

static void print_diff_sets(const DfsGraph *dfs)
{
    for(int i=0;i<dfs->DIFF.nb_groups;i++)
    {
        int p=dfs->DIFF.partner[i];

        if(p!=-1 && i<p)
            printf("G%d <-> G%d\n",i,p);
    }
}

void compute_same_diff(DfsGraph *dfs)    {

    int m = dfs->edge_count;
    dfs->same_parent  = malloc(m * sizeof(int));
    dfs->same_partner = malloc(m * sizeof(int));
    for (int i = 0; i < m; i++) {
        dfs->same_parent[i]  = i;
        dfs->same_partner[i] = -1;
    }

    AttList *att = calloc(m, sizeof(AttList));

    for (int i = 0; i < m; i++)
        if (dfs->edges[i].type == BACK_EDGE)
            build_back_edge_att(dfs, att, i);

    recurse(dfs->root, dfs, att);

    for (int i = 0; i < m; i++)
        if (att[i].pairs != NULL) attlist_free(&att[i]);
    free(att);

    printf("\n================ SAME ================\n");
    print_same_sets(dfs);

    printf("\n================ DIFF ================\n");
    print_diff_sets(dfs);

    /* extraction finale dans dfs->SAME et dfs->DIFF */
    build_SAME_DIFF(dfs);
}

static void print_same_diff_prime(const DfsGraph *dfs)
{
    printf("SAME'\n");

    for(int g=0;g<dfs->SAME_prime.nb_groups;g++)
    {
        printf("G'%d : ",g);

        SameGroup *grp=&dfs->SAME_prime.groups[g];

        for(int i=0;i<grp->size;i++)
        {
            int e=grp->edges[i];

            printf("%d(%d-%d) ",
                e,
                dfs->edges[e].from,
                dfs->edges[e].to);
        }

        printf("\n");
    }

    printf("\nDIFF'\n");

    for(int i=0;i<dfs->DIFF_prime.nb_groups;i++)
    {
        int p=dfs->DIFF_prime.partner[i];

        if(p!=-1 && i<p)
            printf("G'%d <-> G'%d\n",i,p);
    }
}

void build_SAME_DIFF_prime(DfsGraph *dfs) {

    int m = dfs->edge_count;

    /* is_rep[e] = 1 si e est representant de son ensemble singulier
       (ou si e n'est pas singulier du tout) */
    int *is_rep = malloc(m * sizeof(int));
    for (int e = 0; e < m; e++) is_rep[e] = 1;

    for (int s = 0; s < dfs->nb_singular_sets; s++)
        for (int j = 1; j < dfs->singular_sets[s].size; j++)
            is_rep[dfs->singular_sets[s].edges[j]] = 0;

    const int ng = dfs->SAME.nb_groups;
    int *new_index = malloc(ng * sizeof(int));
    int *rep_count = calloc(ng, sizeof(int));

    for (int i = 0; i < ng; i++) {
        new_index[i] = -1;
        SameGroup *g = &dfs->SAME.groups[i];
        for (int j = 0; j < g->size; j++)
            if (is_rep[g->edges[j]]) rep_count[i]++;
    }

    /* un groupe d'origine ne survit dans SAME' que s'il contient
       au moins un representant */
    int nb_groups_prime = 0;
    for (int i = 0; i < ng; i++)
        if (rep_count[i] > 0) new_index[i] = nb_groups_prime++;

    dfs->SAME_prime.nb_groups = nb_groups_prime;
    dfs->SAME_prime.groups = malloc(nb_groups_prime * sizeof(SameGroup));

    for (int i = 0; i < ng; i++) {
        if (new_index[i] == -1) continue;
        SameGroup *dst = &dfs->SAME_prime.groups[new_index[i]];
        dst->edges = malloc(rep_count[i] * sizeof(int));
        dst->size = 0;

        SameGroup *src = &dfs->SAME.groups[i];
        for (int j = 0; j < src->size; j++)
            if (is_rep[src->edges[j]]) dst->edges[dst->size++] = src->edges[j];
    }

    /* DIFF' : remap des partenaires ; un partenaire dont le groupe
       a disparu (n'avait que des non-representants) devient "libre" */
    dfs->DIFF_prime.nb_groups = nb_groups_prime;
    dfs->DIFF_prime.partner = malloc(nb_groups_prime * sizeof(int));

    for (int i = 0; i < ng; i++) {
        if (new_index[i] == -1) continue;
        int p = dfs->DIFF.partner[i];
        int np = (p == -1 || new_index[p] == -1) ? -1 : new_index[p];
        dfs->DIFF_prime.partner[new_index[i]] = np;
    }

    /* FIX (bug #4) : on conserve ce mapping au lieu de le jeter -- il est
       indispensable pour traduire un indice SAME en indice SAME' plus tard
       (apply_partition_bits). Auparavant `free(new_index)` ici faisait
       perdre cette information, et apply_partition_bits melangeait les
       deux espaces d'indices. */
    dfs->same_to_prime = new_index;

    free(is_rep); free(rep_count);

    printf("\n================ SAME' / DIFF' ================\n");
    print_same_diff_prime(dfs);

}

static void print_reduced_system(const ReducedSystem *rs)
{
    printf("DIFF non triviales\n");

    for(int i=0;i<rs->npairs;i++)
        printf("G'%d <-> G'%d\n",
            rs->pair_g1[i],
            rs->pair_g2[i]);

    printf("\nGroupes libres\n");

    for(int i=0;i<rs->nfree;i++)
        printf("G'%d\n",
            rs->free_g[i]);
}

ReducedSystem build_reduced_system(DfsGraph *dfs) {

    ReducedSystem rs;
    int ng = dfs->SAME_prime.nb_groups;
    rs.pair_g1 = malloc(ng * sizeof(int));
    rs.pair_g2 = malloc(ng * sizeof(int));
    rs.free_g  = malloc(ng * sizeof(int));
    rs.npairs = 0;
    rs.nfree  = 0;

    int *done = calloc(ng, sizeof(int));

    for (int i = 0; i < ng; i++) {
        if (done[i]) continue;
        int p = dfs->DIFF_prime.partner[i];
        if (p == -1) {
            rs.free_g[rs.nfree++] = i;
            done[i] = 1;
        } else {
            rs.pair_g1[rs.npairs] = i;
            rs.pair_g2[rs.npairs] = p;
            rs.npairs++;
            done[i] = 1;
            done[p] = 1;
        }
    }

    free(done);

    printf("\n================ REDUCED SYSTEM ================\n");
    print_reduced_system(&rs);

    return rs;
}

/* FIX (bug #1) : fonction appelee dans enumerate_embeddings mais absente
   du fichier d'origine -- erreur de lien garantie sans elle. */
static void free_reduced_system(ReducedSystem *rs) {
    free(rs->pair_g1);
    free(rs->pair_g2);
    free(rs->free_g);
}

static void print_rotation_system(const PlanarEmbedding *emb)
{
    for(int v=0;v<emb->n;v++)
    {
        printf("%2d : ",v);

        for(int i=0;i<emb->rotation_lengths[v];i++)
            printf("%d ",
                emb->rotation_system[v][i]);

        printf("\n");
    }
}

static void apply_partition_bits(DfsGraph *dfs, ReducedSystem *rs, uint64_t k, DynPartition *part) {

    int ngp = dfs->SAME_prime.nb_groups;
    int *side_of_group = malloc(ngp * sizeof(int));
    for (int i = 0; i < ngp; i++) side_of_group[i] = -1;

    int bit = 0;
    for (int i = 0; i < rs->npairs; i++, bit++) {
        int b = (k >> bit) & 1;
        side_of_group[rs->pair_g1[i]] = b;
        side_of_group[rs->pair_g2[i]] = 1 - b;
    }
    for (int i = 0; i < rs->nfree; i++, bit++) {
        side_of_group[rs->free_g[i]] = (k >> bit) & 1;
    }

    int m = dfs->edge_count;
    part->side = malloc(m * sizeof(int));
    for (int e = 0; e < m; e++) part->side[e] = -1;

    for (int i = 0; i < dfs->E0_size; i++) {
        const int e = dfs->E0[i];
        const int gi = dfs->edge_to_group[e];      /* indice SAME (original) */
        const int gip = dfs->same_to_prime[gi];     /* -> indice SAME' */
        /* gip == -1 : e est un membre non-representant d'un ensemble
           singulier libre, dont le groupe a disparu de SAME'. Son cote
           n'est de toute facon jamais lu (cf build_Mv_prime, qui saute
           ces arcs et les gere via la variante choisie). */
        part->side[e] = (gip == -1) ? -1 : side_of_group[gip];
    }
    free(side_of_group);
}

void enumerate_partitions(DfsGraph *dfs, ReducedSystem *rs, void (*callback)(DfsGraph*,
    DynPartition*, uint64_t, void*), void *ctx) {
    int bits = rs->npairs + rs->nfree;
    uint64_t total = (uint64_t)1 << bits;

    for (uint64_t k = 0; k < total; k++) {
        DynPartition part;
        apply_partition_bits(dfs, rs, k, &part);
        callback(dfs, &part, k, ctx);
        free(part.side);
    }
}

static int popcount_u(unsigned x) {
    int c = 0;
    while (x) { c += x & 1; x >>= 1; }
    return c;
}

static void permute_rec(int *arr, const int l, const int r, void (*cb)(int*, int, void*), void *ctx) {

    if (arr==NULL) return;

    if (l == r) { cb(arr, r + 1, ctx); return; }
    for (int i = l; i <= r; i++) {
        int t = arr[l]; arr[l] = arr[i]; arr[i] = t;
        permute_rec(arr, l + 1, r, cb, ctx);
        t = arr[l]; arr[l] = arr[i]; arr[i] = t;
    }
}

static int singular_is_bound(const DfsGraph *dfs, const SingularSet *set) {
    const int r = uf_find(dfs->same_parent, set->edges[0]);
    for (int i = 1; i < set->size; i++)
        if (uf_find(dfs->same_parent, set->edges[i]) != r) return 0;
    return 1;
}

/* ---- cas Bound ---- */
typedef struct {
    void (*user_cb)(SingularVariant*, void*); void *user_ctx;
} BoundCtx;

static void emit_bound(int *arr, int n, void *ctx) {
    BoundCtx *c = ctx;
    SingularVariant v;
    v.left = arr; v.nleft = n;
    v.right = NULL; v.nright = 0;
    c->user_cb(&v, c->user_ctx);
}

/* ---- cas Free ---- */
typedef struct {
    int *s1; int s1n;
    void (*user_cb)(SingularVariant*, void*); void *user_ctx;
} S2Ctx;

static void emit_with_s2(int *s2perm, int n2, void *ctx) {
    S2Ctx *c = ctx;
    SingularVariant v;
    v.left = c->s1; v.nleft = c->s1n;
    v.right = s2perm; v.nright = n2;
    c->user_cb(&v, c->user_ctx);
}

typedef struct {
    int *s2; int s2n;
    void (*user_cb)(SingularVariant*, void*); void *user_ctx;
} S1Ctx;

static void emit_with_s1(int *s1perm, int n1, void *ctx) {
    S1Ctx *c = ctx;
    S2Ctx c2; c2.s1 = s1perm; c2.s1n = n1;
    c2.user_cb = c->user_cb; c2.user_ctx = c->user_ctx;

    int *s2copy = malloc(c->s2n * sizeof(int));
    memcpy(s2copy, c->s2, c->s2n * sizeof(int));
    permute_rec(s2copy, 0, c->s2n - 1, emit_with_s2, &c2);
    free(s2copy);
}

/* genere les h(x) variantes d'un ensemble singulier (set->edges[0] = x) */
void generate_singular_variants(DfsGraph *dfs, SingularSet *set, void (*cb)(SingularVariant*, void*),
    void *ctx) {
    int g = set->size;
    int rep = set->edges[0];
    int others_n = g - 1;
    int *others = malloc(others_n * sizeof(int));
    for (int i = 0; i < others_n; i++) others[i] = set->edges[i + 1];

    if (singular_is_bound(dfs, set)) {
        int *block = malloc(g * sizeof(int));
        block[0] = rep;
        for (int i = 0; i < others_n; i++) block[i + 1] = others[i];

        BoundCtx bc = { cb, ctx };
        permute_rec(block, 0, g - 1, emit_bound, &bc);
        free(block);
    } else {
        for (unsigned mask = 0; mask < (1u << others_n); mask++) {
            int k = popcount_u(mask);
            int s1n = g - k, s2n = k;
            int *s1 = malloc(s1n * sizeof(int));
            int *s2 = malloc(s2n * sizeof(int));
            int i1 = 0, i2 = 0;
            s1[i1++] = rep;
            for (int i = 0; i < others_n; i++) {
                if (mask & (1u << i)) s2[i2++] = others[i];
                else s1[i1++] = others[i];
            }
            S1Ctx c1 = { s2, s2n, cb, ctx };
            permute_rec(s1, 0, s1n - 1, emit_with_s1, &c1);
            free(s1); free(s2);
        }
    }
    free(others);
}

static DfsGraph *g_cmp_dfs;

static int cmp_phi_desc(const void *a, const void *b) {
    int A = *(const int*)a, B = *(const int*)b;
    int pa = g_cmp_dfs->edges[A].phi, pb = g_cmp_dfs->edges[B].phi;
    if (pa != pb) return pb - pa; /* decroissant */
    return A - B; /* tie-break stable */
}

static int cmp_phi_asc(const void *a, const void *b) {
    int A = *(const int*)a, B = *(const int*)b;
    int pa = g_cmp_dfs->edges[A].phi, pb = g_cmp_dfs->edges[B].phi;
    if (pa != pb) return pa - pb;
    return A - B;
}

/* chosen[s] = la SingularVariant choisie pour le set singulier s
   (NULL si non singulier) */
PartialMap build_Mv_prime(DfsGraph *dfs, DynPartition *part, SingularVariant **chosen,
    int *edge_set_id, int v) {
    DfsVertex *vv = &dfs->vertices[v];
    int ref = vv->out_edges[0];

    int cap = vv->nb_out_edges * 4 + 4;
    int *L = malloc(cap * sizeof(int)); int nL = 0;
    int *R = malloc(cap * sizeof(int)); int nR = 0;

    for (int i = 0; i < vv->nb_out_edges; i++) {
        int e = vv->out_edges[i];
        if (e == ref) continue;

        int s = edge_set_id[e];
        if (s != -1 && dfs->singular_sets[s].edges[0] != e) continue;
        /* (saute les membres non-representants : geres via le rep) */

        if (s != -1) {
            SingularVariant *var = chosen[s];
            int side_x = part->side[e]; /* 0 = LL, 1 = RR */
            int *same_side  = (side_x == 0) ? L : R;
            int *opp_side   = (side_x == 0) ? R : L;
            int *n_same = (side_x == 0) ? &nL : &nR;
            int *n_opp  = (side_x == 0) ? &nR : &nL;

            for (int j = 0; j < var->nleft; j++)
                same_side[(*n_same)++] = var->left[j];
            for (int j = 0; j < var->nright; j++)
                opp_side[(*n_opp)++] = var->right[j];
        } else {
            if (part->side[e] == 0) L[nL++] = e;
            else                    R[nR++] = e;
        }
    }

    g_cmp_dfs = dfs;
    qsort(L, nL, sizeof(int), cmp_phi_desc);
    qsort(R, nR, sizeof(int), cmp_phi_asc);

    PartialMap m;
    m.n = nL + 1 + nR;
    m.edges = malloc(m.n * sizeof(int));
    int idx = 0;
    for (int i = 0; i < nL; i++) m.edges[idx++] = L[i];
    m.edges[idx++] = ref;
    for (int i = 0; i < nR; i++) m.edges[idx++] = R[i];

    free(L); free(R);
    return m;
}

static void assign_labels(DfsGraph *dfs, PartialMap *mprime, int v, int *next_label) {
    PartialMap *mv = &mprime[v];
    for (int i = 0; i < mv->n; i++) {
        int e = mv->edges[i];
        if (dfs->edges[e].type == BACK_EDGE) {
            dfs->edges[e].label = (*next_label)++;
        } else {
            assign_labels(dfs, mprime, dfs->edges[e].to, next_label);
        }
    }
}

void compute_labels(DfsGraph *dfs, PartialMap *mprime) {
    int next_label = 1;
    assign_labels(dfs, mprime, dfs->root, &next_label);
}

static int compute_back(DfsGraph *dfs, int *back_cache, int e) {
    if (back_cache[e] != -1) return back_cache[e];
    if (dfs->edges[e].type == BACK_EDGE) { back_cache[e] = e; return e; }
    int w = dfs->edges[e].to;
    int eref = dfs->vertices[w].out_edges[0];
    int r = compute_back(dfs, back_cache, eref);
    back_cache[e] = r;
    return r;
}

static int compute_e_in(DfsGraph *dfs, int *back_cache, int v) {
    if (v != dfs->root) return dfs->vertices[v].parent_edge;
    int eref = dfs->vertices[v].out_edges[0];
    return compute_back(dfs, back_cache, eref);
}

static int cmp_label_desc(const void *a, const void *b) {
    int A = *(const int*)a, B = *(const int*)b;
    return g_cmp_dfs->edges[B].label - g_cmp_dfs->edges[A].label;
}

VertexRotation *build_M_all(DfsGraph *dfs, PartialMap *mprime) {
    int n = dfs->vertices_count;
    int m = dfs->edge_count;

    int *back_cache = malloc(m * sizeof(int));
    for (int i = 0; i < m; i++) back_cache[i] = -1;

    /* back_into[v] = arcs de retour entrant en v */
    int *back_count = calloc(n, sizeof(int));
    for (int i = 0; i < m; i++)
        if (dfs->edges[i].type == BACK_EDGE)
            back_count[dfs->edges[i].to]++;

    int **back_into = malloc(n * sizeof(int*));
    int *fill = calloc(n, sizeof(int));
    for (int v = 0; v < n; v++)
        back_into[v] = malloc(back_count[v] * sizeof(int));
    for (int i = 0; i < m; i++)
        if (dfs->edges[i].type == BACK_EDGE) {
            int v = dfs->edges[i].to;
            back_into[v][fill[v]++] = i;
        }

    VertexRotation *rot = malloc(n * sizeof(VertexRotation));
    g_cmp_dfs = dfs;

    for (int v = 0; v < n; v++) {
        int ein = compute_e_in(dfs, back_cache, v);
        /* DEBUG: ein doit être un arc incident à v */
        {
            int ef = dfs->edges[ein].from, et = dfs->edges[ein].to;
            if (ef != v && et != v)
                printf( "*** BUG build_M_all: v=%d ein=arc%d (%d→%d) NOT incident!\n",
                        v, ein, ef, et);
        }
        int cap = 8, cnt = 0;
        int *buf = malloc(cap * sizeof(int));
        buf[cnt++] = ein;

        PartialMap *mv = &mprime[v];

        for (int i = 0; i < mv->n; i++) {
            int e = mv->edges[i];

            if (dfs->edges[e].type == BACK_EDGE) {
                if (cnt >= cap) { cap *= 2; buf = realloc(buf, cap * sizeof(int)); }
                buf[cnt++] = e;
                continue;
            }

            int back_e = compute_back(dfs, back_cache, e);
            int lbl_back = dfs->edges[back_e].label;

            int *Lc = malloc(back_count[v] * sizeof(int)); int nl = 0;
            int *Rc = malloc(back_count[v] * sizeof(int)); int nr = 0;

            for (int j = 0; j < back_count[v]; j++) {
                int b = back_into[v][j];

                /* FIX (bug #3) : il fallait exclure e_{v,in} (calcule via
                   compute_e_in, deja disponible dans `ein` plus haut), pas
                   e_{v,ref} (eref_v). Pour tout sommet non-racine, e_{v,in}
                   est un arc d'ARBRE, donc il n'apparaissait de toute facon
                   jamais parmi les arcs de retour `b` -- le bug etait
                   invisible. Pour la RACINE en revanche, e_{racine,in} EST
                   un arc de retour (celui qui ferme cycle(e_{racine,ref})),
                   et l'ancien test `b == eref_v` (un arc sortant, jamais un
                   arc de retour entrant) ne l'excluait jamais : il se
                   retrouvait insere deux fois dans M(racine). Cf. l'exemple
                   du document : entrant(a) = {f}, pas {f,i}, car i =
                   e_{1,in} doit etre exclu. */
                if (b == ein) continue;

                if (!is_descendant_edge(dfs, e, b)) continue;
                if (dfs->edges[b].label < lbl_back) Lc[nl++] = b;
                else                                Rc[nr++] = b;
            }
            qsort(Lc, nl, sizeof(int), cmp_label_desc);
            qsort(Rc, nr, sizeof(int), cmp_label_desc);

            int needed = cnt + nl + 1 + nr;
            while (needed > cap) { cap *= 2; buf = realloc(buf, cap * sizeof(int)); }

            for (int j = 0; j < nl; j++) buf[cnt++] = Lc[j];
            buf[cnt++] = e;
            for (int j = 0; j < nr; j++) buf[cnt++] = Rc[j];

            free(Lc); free(Rc);
        }

        rot[v].order = buf;
        rot[v].n = cnt;
    }

    for (int v = 0; v < n; v++) free(back_into[v]);
    free(back_into); free(back_count); free(fill); free(back_cache);

    return rot;
}


/* une fois TOUTES les variantes choisies pour TOUS les ensembles singuliers :
   on construit la carte planaire complete et on la stocke */
static void emit_embedding(const VariantCtx *ctx) {
    DfsGraph *dfs = ctx->dfs;
    const int n = dfs->vertices_count;

    PartialMap *mprime = malloc(n * sizeof(PartialMap));
    for (int v = 0; v < n; v++)
        mprime[v] = build_Mv_prime(dfs, ctx->part, ctx->chosen,
                                    ctx->edge_set_id, v);

    compute_labels(dfs, mprime);
    VertexRotation *rot = build_M_all(dfs, mprime);

    /* stockage PERMANENT dans result */
    EmbeddingSet *res = ctx->result;
    res->embeddings = realloc(res->embeddings,
                               (res->count + 1) * sizeof(PlanarEmbedding));
    PlanarEmbedding *emb = &res->embeddings[res->count++];
    emb->n = n;
    emb->rotation_system = malloc(n * sizeof(int*));
    emb->rotation_lengths = malloc(n * sizeof(int));  /* FIX (bug #7) */
    for (int v = 0; v < n; v++) {
        emb->rotation_system[v] = malloc(rot[v].n * sizeof(int));
        memcpy(emb->rotation_system[v], rot[v].order, rot[v].n * sizeof(int));
        emb->rotation_lengths[v] = rot[v].n;          /* FIX (bug #7) */
    }

    printf("\n================ ROTATION ================\n");
    print_rotation_system(emb);

    for (int v = 0; v < n; v++) { free(mprime[v].edges); free(rot[v].order); }
    free(mprime); free(rot);
}

/* produit cartesien des variantes de chaque ensemble singulier */

/* FIX (bug #2) : `variant_cb` n'avait qu'un prototype (jamais de corps), et
   s'appuyait sur une variable globale `ctx_set_idx_holder` jamais declaree
   -- erreur de lien garantie. On la remplace par un petit contexte local
   explicite (plus robuste, pas de risque si l'ordre d'appel change un jour). */
static void variant_recurse(VariantCtx *ctx, int set_idx);

typedef struct {
    VariantCtx *ctx;
    int set_idx;
} VariantCbCtx;

static void variant_cb(SingularVariant *v, void *vctx) {
    VariantCbCtx *vc = vctx;
    vc->ctx->chosen[vc->set_idx] = v;
    variant_recurse(vc->ctx, vc->set_idx + 1);
}

static void variant_recurse(VariantCtx *ctx, int set_idx) {
    DfsGraph *dfs = ctx->dfs;
    if (set_idx == dfs->nb_singular_sets) {
        emit_embedding(ctx);   /* tout est choisi -> on construit la carte */
        return;
    }
    VariantCbCtx vc = { ctx, set_idx };
    generate_singular_variants(dfs, &dfs->singular_sets[set_idx],
                                variant_cb, &vc);
}

static void print_partition(const DfsGraph *dfs,
                            const DynPartition *part)
{
    for(int g=0;g<dfs->SAME_prime.nb_groups;g++)
    {
        printf("G'%d = %d\n",
            g,
            part->side[g]);
    }
}

/* callback appele a chaque partition reduite (ephemere, voir plus haut) */
static void on_partition(DfsGraph *dfs, DynPartition *part, uint64_t k, void *ctx) {
    (void)k;
    EmbedCtx *ec = ctx;

    VariantCtx vc;
    vc.dfs = dfs;
    vc.part = part;                 /* part vit le temps de cet appel */
    vc.chosen = calloc(dfs->nb_singular_sets, sizeof(SingularVariant*));
    vc.edge_set_id = ec->edge_set_id;
    vc.result = ec->result;

    printf("\n================ h ================\n");
    print_partition(dfs,part);

    variant_recurse(&vc, 0);

    free(vc.chosen);
}

EmbeddingSet enumerate_embeddings(DfsGraph *dfs) {

    EmbeddingSet result = { NULL, 0 };

    int m = dfs->edge_count;
    int *edge_set_id = malloc(m * sizeof(int));
    for (int i = 0; i < m; i++) edge_set_id[i] = -1;
    for (int s = 0; s < dfs->nb_singular_sets; s++)
        for (int j = 0; j < dfs->singular_sets[s].size; j++)
            edge_set_id[dfs->singular_sets[s].edges[j]] = s;

    EmbedCtx ec = { dfs, edge_set_id, &result };

    ReducedSystem rs = build_reduced_system(dfs);
    print_reduced_system(&rs);
    enumerate_partitions(dfs, &rs, on_partition, &ec);
    free_reduced_system(&rs);

    free(edge_set_id);
    return result;
}

/* FIX (bug #8) : aucune fonction ne liberait un EmbeddingSet -- fuite
   memoire garantie a chaque appel de generate_all_embeddings. */
void free_embedding_set(EmbeddingSet *es) {
    if (!es) return;
    for (int i = 0; i < es->count; i++) {
        PlanarEmbedding *emb = &es->embeddings[i];
        for (int v = 0; v < emb->n; v++) free(emb->rotation_system[v]);
        free(emb->rotation_system);
        free(emb->rotation_lengths);
    }
    free(es->embeddings);
    es->embeddings = NULL;
    es->count = 0;
}

EmbeddingSet generate_all_embeddings(Graph *g)
{
    DfsGraph *dfs = build_dfs_graph(g);

    compute_low_values(dfs);

    build_phi_lists(dfs);

    build_singular_sets(dfs);

    compute_same_diff(dfs);

    build_SAME_DIFF_prime(dfs);

    const EmbeddingSet result = enumerate_embeddings(dfs);

    free_dfs_graph(dfs);

    return result;
}

/* planar_map_generator.c */

long long factorial(int n) {
    long long r = 1;
    for (int i = 2; i <= n; i++) r *= i;
    return r;
}

/*
 * Calcule le nombre theorique de cartes planaires selon le Theoreme 1 de Cai :
 *
 *   2^(d+s) x produit_{x in RS} h(x)
 *
 * avec :
 *   d   = nombre de paires non-triviales dans DIFF'
 *   s   = nombre de groupes libres dans SAME'
 *   h(x) = g(x)!         si x est "bound"  (tout son ensemble singulier
 *                                            est dans le meme groupe SAME)
 *   h(x) = (g(x)+1)! / 2 si x est "free"
 *   g(x) = taille de l'ensemble singulier de x
 */
long long theoretical_embedding_count(const DfsGraph *dfs) {

    /* --- 1. Calculer d et s depuis DIFF' et SAME' --- */
    int d = 0, s = 0;
    int ng = dfs->SAME_prime.nb_groups;
    int *done = calloc(ng, sizeof(int));

    for (int i = 0; i < ng; i++) {
        if (done[i]) continue;
        int p = dfs->DIFF_prime.partner[i];
        if (p == -1) {
            s++;
            done[i] = 1;
        } else {
            d++;
            done[i] = 1;
            done[p] = 1;
        }
    }
    free(done);

    /* --- 2. Calculer le produit des h(x) sur RS --- */
    long long product = 1;

    for (int si = 0; si < dfs->nb_singular_sets; si++) {
        const SingularSet *set = &dfs->singular_sets[si];
        int g_x = set->size;

        /* Le representant est set->edges[0] par convention */
        int rep = set->edges[0];
        int gi_rep = dfs->edge_to_group[rep];

        /* x est "bound" si tous les arcs du set sont dans le meme groupe SAME */
        int is_bound = 1;
        for (int j = 1; j < set->size; j++) {
            if (dfs->edge_to_group[set->edges[j]] != gi_rep) {
                is_bound = 0;
                break;
            }
        }

        long long hx = is_bound ? factorial(g_x)
                                : factorial(g_x + 1) / 2;
        product *= hx;
    }

    /* --- 3. Resultat final --- */
    long long power = (long long)1 << (d + s);
    return power * product;
}

/*
 * Verifie que actual_count == nombre theorique.
 * Affiche un rapport et renvoie 1 si OK, 0 sinon.
 */
int verify_embedding_count(const DfsGraph *dfs, const int actual_count) {

    long long expected = theoretical_embedding_count(dfs);

    /* Recalculer d et s pour l'affichage */
    int d = 0, s = 0;
    int ng = dfs->SAME_prime.nb_groups;
    int *done = calloc(ng, sizeof(int));
    for (int i = 0; i < ng; i++) {
        if (done[i]) continue;
        int p = dfs->DIFF_prime.partner[i];
        if (p == -1) { s++; done[i] = 1; }
        else         { d++; done[i] = 1; done[p] = 1; }
    }
    free(done);

    long long prod = expected / ((long long)1 << (d + s));

    printf("=== Verification du nombre de plongements ===\n");
    printf("  d (paires DIFF non-triviales) = %d\n", d);
    printf("  s (groupes SAME libres)        = %d\n", s);
    printf("  2^(d+s)                        = %lld\n", (long long)1 << (d + s));
    printf("  prod h(x)                      = %lld\n", prod);
    printf("  Attendu  : %lld\n", expected);
    printf("  Obtenu   : %d\n",   actual_count);

    if ((long long)actual_count == expected) {
        printf("  [OK] Compte correct.\n");
        return 1;
    }
    printf("  [ERREUR] Difference de %lld plongements.\n",
           expected - (long long)actual_count);
    return 0;
}
