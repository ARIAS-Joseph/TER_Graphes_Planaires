#ifndef TER_GRAPHES_PLANAIRES_PLANAR_GRAPH_CREATOR_H
#define TER_GRAPHES_PLANAIRES_PLANAR_GRAPH_CREATOR_H
#include "graph.h"

Graph* create_outer_planar_graph(int nb_vertex, int nb_edges_target);

Graph* create_planar_graph(int nb_vertex, int nb_edges_target);

void create_tree(Graph* g);

Graph* test_circle(int nb_vertex);

#endif //TER_GRAPHES_PLANAIRES_PLANAR_GRAPH_CREATOR_H
