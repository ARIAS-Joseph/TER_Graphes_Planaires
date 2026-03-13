#include "permutations.h"

#include <stdlib.h>

/**
 * @brief Give the permutation corresponding to a given inversion table.
 *
 * @param inv_tab The inversion table, an array of size n where inv_tab[i] is the number of elements
 * before i that are greater than i in the permutation.
 * @param perm_tab The array where the resulting permutation will be stored. It should be an array
 * of size n full of 0's.
 * @param n The size of the permutation (number of elements)
 * @return The resulting permutation, an array of size n where perm_tab[i] is the element at
 * position i in the permutation. The elements of the permutation are the integers from 1 to n,
 * each appearing exactly once.
 */
int* inversion_to_permutation(const int *inv_tab, int* perm_tab, const int n) {
    for (int i = 0; i < n; i++) {
        int k = 0; // We look for the k-th position in the permutation that is still empty
        // (perm_tab[k] == 0) and such that there are inv_tab[i] empty positions before it
        // (count < inv_tab[i])
        int count = 0; // Count the number of empty positions before the k-th position
        while (perm_tab[k] != 0 || count < inv_tab[i]) {
            if (perm_tab[k] == 0) {
                count++;
            }
            k++;
        }
        perm_tab[k] = i+1;
    }
    return perm_tab;
}

/**
 * @brief Generate the next inversion table in lexicographic order.
 *
 * @param inv The current inversion table. Should be an array of size n.
 * @param n The size of the inversion table (number of elements in the permutation).
 * @return 1 if the next inversion table was generated successfully, 0 if all permutations have been
 * generated.
 */
int next_inversion_table(int *inv, const int n) {

    for (int i = 0; i < n-1; i++) {
        if (inv[i] < n-1 - i) {
            inv[i]++;
            return 1;
        }
        inv[i] = 0;
    }

    return 0;
}

void generate_random_inversion_table(int *inv, const int n) {
    for (int i = 0; i < n-1; i++) {
        inv[i] = rand() % (n - 1 - i);
    }
    inv[n-1] = 0;
}

/*
 * @brief Create an array containing all pairs of distinct vertices possibles of a graph with
 * nb_vertex vertices. he pairs are stored in the array perm, where perm[2*k] and perm[2*k+1] are
 * the two vertices of the k-th edge. The edges are ordered in lexicographic order, meaning that the
 * pairs are sorted first by the first vertex and then by the second vertex.
 *
 * @param nb_vertex The number of vertices in the complete graph.
 * @param perm The array where the pairs of vertices will be stored. It should be an array of size
 * nb_vertex*(nb_vertex-1).
 */
void create_all_edges(const int nb_vertex, int* perm) {
    int k = 0;
    for (int i = 0; i < nb_vertex; i++) {
        for (int j = i + 1; j < nb_vertex; j++) {
            perm[k++] = i;
            perm[k++] = j;
        }
    }
}

void fisher_yates_shuffle(int* array, const int n) {
    for (int i = n - 1; i > 0; i -= 2) {

        const int pair = rand() % ((i + 1) / 2);
        const int j = 2 * pair + 1;

        const int temp = array[i];
        array[i] = array[j];
        array[j] = temp;

        const int temp2 = array[i - 1];
        array[i - 1] = array[j - 1];
        array[j - 1] = temp2;
    }
}