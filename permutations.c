#include "permutations.h"

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