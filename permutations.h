#ifndef TER_GRAPHES_PLANAIRES_PERMUTATIONS_H
#define TER_GRAPHES_PLANAIRES_PERMUTATIONS_H

int* inversion_to_permutation(const int *inv_tab, int* perm_tab, int n);

int next_inversion_table(int *inv, int n);

void create_all_edges(int nb_vertex, int* perm);

void fisher_yates_shuffle(int* array, int n);

void generate_random_inversion_table(int *inv, int n);

#endif //TER_GRAPHES_PLANAIRES_PERMUTATIONS_H