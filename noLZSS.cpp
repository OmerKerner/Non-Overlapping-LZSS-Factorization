/*
Compute the Non-Overlapping-LZSS-Factorization of a string
Copyright (c) 2024 Dr. Omer Kerner

Acknowledgements: This program is based on the pseudocode created by
Dominik Köppl in his paper <https://doi.org/10.3390/a14020044>.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <sdsl/suffix_trees.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace sdsl;

typedef cst_sada<> cst_t;

/**
 * @brief Report a new factor.
 *
 * @param j_sufnum The suffix number of the reference leaf. Fresh factors reference themselves.
 * @param witness The witness node of the reference leaf.
 * @param lambda_sufnum The suffix number of the beginning of the new factor.
 * @param l The length of the new factor.
 */
static void report(size_t j_sufnum, cst_t::node_type witness, size_t lambda_sufnum, size_t l) {
}

size_t lcp(cst_t& cst, size_t i, size_t j) {
    // use the lcp array of cst to calc the lcp of the suffixes i and j
    if (i == j) return cst.csa.size() - cst.csa[i];
    auto lca = cst.lca(cst.select_leaf(cst.csa.isa[i]+1), cst.select_leaf(cst.csa.isa[j]+1));
    return cst.depth(lca);
}

/**
 * @brief Get the next leaf of the suffix tree in text order.
 *
 * @param cst The suffix tree to search in.
 * @param lambda The current leaf.
 * @param iterations The number of iterations to apply to the Psi function (default is 1).
 * @return The next leaf in text order.
 */
static cst_t::node_type next_leaf(cst_t& cst, cst_t::node_type lambda, size_t iterations = 1)
{
    assert(cst.is_leaf(lambda)); // asert lambda is a leaf
    auto lambda_rank = cst.lb(lambda); // left boundary of a leaf is the rank of the leaf itself
    // get the Psi of lambda rank. If iterations is > 1, we calc psi^iterations(lambda_rank)
    for (size_t i = 0; i < iterations; i++) {
        lambda_rank = cst.csa.psi[lambda_rank];
    }
    auto next_leaf = cst.select_leaf(lambda_rank + 1);
    return next_leaf;
}


/**
 * @brief Compute the number of Non-Overlapping-LZSS factors of a string.
 *
 * @param cst The compressed suffix tree of the string.
 * @return The number of factors minus one.
 */
size_t lzss(cst_t& cst) {
    rmq_succinct_sct<> rmq(&cst.csa);
    size_t str_len = cst.size() - 1; // the length of the string is the size of the CST minus the '$' character

    // get the leaf of the first suffix
    auto lambda = cst.select_leaf(cst.csa.isa[0] + 1);
    size_t lambda_node_depth = cst.node_depth(lambda);
    size_t lambda_sufnum = 0;

    cst_t::node_type v;
    size_t v_min_leaf_sufnum = 0;
    size_t u_min_leaf_sufnum = 0;
    size_t num_of_factors = 0;

    while (lambda_sufnum < str_len) { // while the whole text is not processed, compute the next factor. The whole text is processed when the last suffix remaining is '$' 
        size_t d = 1;
        size_t l = 1;
        while (true) {
            v = cst.bp_support.level_anc(lambda, lambda_node_depth - d);
            v_min_leaf_sufnum = cst.csa[rmq(cst.lb(v), cst.rb(v))];
            l = cst.depth(v);
            if (v_min_leaf_sufnum + l - 1 < lambda_sufnum) {
                u_min_leaf_sufnum = v_min_leaf_sufnum;
                d++;
                continue;
            }
            auto u = cst.parent(v);
            auto u_depth = cst.depth(u);
            if (v_min_leaf_sufnum == lambda_sufnum) {
                if (u == cst.root()) {
                    l = 1;
                    report(lambda_sufnum, u, lambda_sufnum, l); // report a fresh factor
                    break;
                }
                else {
                    l = u_depth;
                    report(u_min_leaf_sufnum, u, lambda_sufnum, l); // report a type 1 factor
                    break;
                }
            }
            l = std::min(lcp(cst, lambda_sufnum, v_min_leaf_sufnum), (lambda_sufnum - v_min_leaf_sufnum));
            if (l <= u_depth) {
                l = u_depth;
                report(u_min_leaf_sufnum, u, lambda_sufnum, l); // report a type 2 factor
                break;
            }
            else {
                report(v_min_leaf_sufnum, v, lambda_sufnum, l); // report a type 3 factor
                break;
            }
        }
        num_of_factors++;
        lambda = next_leaf(cst, lambda, l);
        lambda_node_depth = cst.node_depth(lambda);
        lambda_sufnum = cst.sn(lambda);
    }
    return num_of_factors - 1;
}

/**
 * @brief Main function that computes the Non-Overlapping-LZSS-Factorization factorization.
 *
 * @param argc The number of command line arguments.
 * @param argv The array of command line arguments. argv[1] should be the path to the file containing the string.
 * @return Returns -1 for indicating an error. Otherwise, returns the number of factors minus one.
 */
int main(int argc, char* argv[])
{
    if (argc < 3) {
        cout << "usage: " << argv[0] << " -f filepath | -s string" << std::endl;
        return -1;
    }

    int opt;
    std::string input;
    bool is_file = false;
    while ((opt = getopt(argc, argv, "f:s:")) != -1) {
        switch (opt) {
        case 'f':
            input = std::string(optarg);
            is_file = true;
            break;
        case 's':
            input = std::string(optarg);
            is_file = false;
            break;
        default:
            cout << "usage: " << argv[0] << " -f filepath | -s string" << std::endl;
            return -1;
        }
    }

    // construct the CST
    cst_t cst;
    if (is_file) {
        construct(cst, input, 1);
    } else {
        construct_im(cst, input, 1);
    }
    size_t num_of_factors = lzss(cst);
    cout << num_of_factors << endl; // should return the number of factors minus 1
    return 0; 
}

extern "C" int nolzss(const char* str) {
    cst_t cst;
    construct_im(cst, str, 1);
    return lzss(cst);
}
