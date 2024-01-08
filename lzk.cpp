/*
LZSSK - Compute the Lempel–Ziv–Storer–Szymanski-Kerner factorization of a string
Copyright (c) 2024 Omer Kerner

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
    cout << "node " << node << " j " << j << " l " << l << endl;
}

//! Next leaf of the suffix tree in text-order
/*! \param cst The suffix tree
    \param lambda The current leaf
    \param iterations The number of iterations to apply to the Psi function
    \return The next leaf  */
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

    size_t ai = -1; // assembly index (number of factors reported - 1)

    while (cst.sn(lambda) < str_len) { // while the whole text is not processed, compute the next factor. The whole text is processed when the last suffix remaining is '$' {
        size_t d = 1;
        size_t l = 1;
        size_t j_u = 0;
        do {
            auto v = cst.bp_support.level_anc(lambda, lambda_node_depth - d);
            auto j_v = cst.csa[rmq(cst.lb(v), cst.rb(v))];
            l = cst.depth(v);
            if (j_v + l - 1 < cst.sn(lambda)) {
                j_u = j_v;
                d++;
                continue;
            }
            auto u = cst.parent(v);
            if (j_v == cst.sn(lambda)) {
                if (u == cst.root()) {
                    l = 1;
                    report(u, j_u, l); // report a fresh factor, u & j_u are 0
                    break;
                }
                l = cst.depth(u);
                report(u, j_u, l); // report type 1 factor
                break;
            }
            l = std::min(cst.depth(cst.lca(lambda, (cst.select_leaf(cst.csa.isa[j_v] + 1)))), (cst.sn(lambda) - j_v));
            cout << "lca " << cst.depth(cst.lca(lambda, (cst.select_leaf(cst.csa.isa[j_v] + 1)))) << " sn-j_v " << (cst.sn(lambda) - j_v) << endl;
            if (l <= cst.depth(u)) { // report type 2 factor
                l = cst.depth(u);
                report(u, j_u, l);
            }
            else {
                report(v, j_v, l); // report type 3 factor
            }
            break;
        } while (d != lambda_node_depth);
        ai++;
        lambda = next_leaf(cst, lambda, l);
        size_t lambda_node_depth = cst.node_depth(lambda);
        cout << "lambda " << lambda << " lambda_node_depth " << lambda_node_depth << " sn " << cst.sn(lambda) << endl;
    } 
    cout << "ai " << ai << endl;
    return ai;
}
