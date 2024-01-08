/*
LZK - Compute the Lempel–Ziv–Kerner (LZK) factorization of a string
Copyright (c) 2024 Omer Kerner

Acknowledgements: This program is based on the pseudocode created by
Dominik Köppl in his paper <https://doi.org/10.3390/a14020044>.
It was changed to compute the LZK factorization instead of the 
non-overlapping LZSS factorization.

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
#include <iostream>
#include <string>

using namespace std;
using namespace sdsl;

typedef cst_sada<> cst_t;


static void report(cst_sada<>::node_type node, int j, int l) {
    // Implement your report function here
    //print the node, j and l
    cout << "node " << node << " reffered sufnum " << j << " len " << l << endl;
}

//!  Get the next leaf of the suffix tree in text order.
/*! \param cst The suffix tree to search in.
    \param lambda The current leaf.
    \param iterations The number of iterations to apply to the Psi function (default is 1).
    \return The next leaf in text order.*/
static cst_t::node_type next_leaf(cst_t& cst, cst_t::node_type lambda, size_t iterations = 1)
{
    // asert lambda is a leaf
    assert(cst.is_leaf(lambda));

    // get the rank of lambda (left boundary of a leaf is the rank of the leaf itself)
    auto lambda_rank = cst.lb(lambda);

    // get the Psi of lambda rank. If iterations is > 1, we calc psi^iterations(lambda_rank)
    auto psi_of_lambda_rank = lambda_rank;
    for (size_t i = 0; i < iterations; i++) {
        // get the Psi of lambda rank
        psi_of_lambda_rank = cst.csa.psi[psi_of_lambda_rank];
    }
    auto next_leaf = cst.select_leaf(psi_of_lambda_rank + 1);
    // return the next leaf
    return next_leaf;
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "usage: " << argv[0] << " file" << std::endl;
        return 1;
    }

    // construct the CST with 1 as alphabet size
    cst_t cst;
    construct(cst, argv[1], 1);
    rmq_succinct_sct<> rmq(&cst.csa);
    auto str_len = cst.size() - 1; // the length of the string is the size of the CST minus the '$' character

    // get the leaf of the first suffix
    auto lambda = cst.select_leaf(cst.csa.isa[0] + 1);
    size_t lambda_node_depth = cst.node_depth(lambda);
    size_t lambda_sufnum = cst.sn(lambda);

    int ai = -1; // assembly index (number of factors reported - 1)

    while (lambda_sufnum < str_len) { // while the whole text is not processed, compute the next factor. The whole text is processed when the last suffix remaining is '$' 
        size_t d = 1;
        size_t l = 1;
        size_t u_min_leaf_sufnum = 0;
        while (true) {
            auto v = cst.bp_support.level_anc(lambda, lambda_node_depth - d);
            auto v_min_leaf_sufnum = cst.csa[rmq(cst.lb(v), cst.rb(v))];
            l = cst.depth(v);
            if (v_min_leaf_sufnum + l - 1 < lambda_sufnum) {
                if (lambda_sufnum + l == str_len) {
                    report(v, v_min_leaf_sufnum, l); // reached the end of the string, report last factor
                    break;
                }
                u_min_leaf_sufnum = v_min_leaf_sufnum;
                d++;
                continue;
            }
            auto u = cst.parent(v);
            auto u_depth = cst.node_depth(u);
            if (v_min_leaf_sufnum == lambda_sufnum) {
                if (u == cst.root()) {
                    l = 1;
                    report(u, u_min_leaf_sufnum, l); // report a fresh factor, u & u_min_leaf_sufnum are 0
                    break;
                }
                l = u_depth;
                report(u, u_min_leaf_sufnum, l); // report type 1 factor
                break;
            }
            l = std::min(cst.depth(cst.lca(lambda, (cst.select_leaf(cst.csa.isa[v_min_leaf_sufnum] + 1)))), (lambda_sufnum - v_min_leaf_sufnum));
            if (l <= u_depth) { // report type 2 factor
                l = u_depth;
                report(u, u_min_leaf_sufnum, l);
            }
            else {
                report(v, v_min_leaf_sufnum, l); // report type 3 factor
            }
            break;
        }
        ai++;
        lambda = next_leaf(cst, lambda, l);
        lambda_node_depth = cst.node_depth(lambda);
        lambda_sufnum = cst.sn(lambda);
    } 
    cout << "ai " << ai << endl;
    return ai;
}
