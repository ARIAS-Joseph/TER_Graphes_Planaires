#ifndef TER_GRAPHES_PLANAIRES_PLANAR_GRAPH_CREATOR_H
#define TER_GRAPHES_PLANAIRES_PLANAR_GRAPH_CREATOR_H
#include "graph.h"

Graph* create_outer_planar_graph(int nb_vertex, int nb_edges_target);

Graph* create_planar_graph(int nb_vertex, int nb_edges_target);

void create_tree(Graph* g, const int* perm);

Graph* create_three_connex_planar_graph(int nb_edges, int nb_edges_target);

void reduce_graph(Graph* g);

#endif //TER_GRAPHES_PLANAIRES_PLANAR_GRAPH_CREATOR_H
