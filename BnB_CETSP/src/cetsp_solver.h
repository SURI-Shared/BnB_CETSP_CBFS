#ifndef CETSP_SOLVER_H
#define CETSP_SOLVER_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <cstdlib>
#include <stdio.h>
#include <cfloat>
#include <unordered_map>
#include <map>
#include <string>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>
#include <assert.h>

#include "util.h"
#include "PrintFunctions.h"
#include "Data.h"
#include "BranchNBound.h"
#include "SolveSocpCplex.h"
#include "structs.h"
#include "CBFS.h"
#include "subgraph.h"
#include "tsp_lb.h"
#include "local_search.h"

using namespace std;

// Stand alone CETSP BnB solver class that allows custom root node assignment
// Intersection tolerance can be set separately as well.
// Additionally, method for upadting lower bound for each node is implemented so that the lower bounds
// is accurate for the search tree after the BnB algo is terminated.
class CETSP_solver
{
public:
    CETSP_solver(Data* d, CETSP_Options& option);
    ~CETSP_solver();

    bool init_root_node();
    void set_root_node(node* root);
    int solve();
    int solve_keep_node();
    int solve_keep_history();
    void set_intersect_tol(double tol);
    void copy_node(node* org_node);
    void clean_up();
    void update_all_lb();
    void reset_root_node();
    double update_all_lb_helper(string current_seq);
    void update_all_lb_by_node_info();
    double update_all_lb_by_node_info_helper(string current_seq);
    bool reduce_uncovered_at_root();

    Data* m_data;
    node* m_root_node;
    string m_root_seq;
    BranchNBound *m_bnbPtr;
    int m_size_inst = 0;
    // Basic CETSP parameters
    double m_overlap_ratio;
    int m_branching_rule;
    int m_branching_strategy;
    double m_time_limit;
    int m_root_select_rule;
    int m_strong_branching_size;
    // CBFS related parameters
    int m_contour_option = 1;
    int m_tie_brak_option = 3;
    int m_uncov_cont_option = 1;
    int m_uncov_cont_param = 1;
    int m_measure_best_mode = 1;
    bool m_use_reverse_prune = false;
    bool m_use_uncov_lb = false;
    // Solution related parameters
    double m_best_sol = DBL_MAX;
    double m_best_ub = DBL_MAX;
    double m_best_lb = 0;
    double m_best_known = DBL_MAX;
    vector<int> m_solution_sequence;
    vector<vector<double>> m_solution_coords;
    // Search tree related parameters
    double m_root_lb;
    double m_gap_root;
    double m_gap_real;
    double m_gap_lb_bnb;
    unsigned long int m_count_SOCP_solved = 0;
    // Time related parameters
    double m_init_total_bnb_time = 0;
    double m_computation_time = 0;
    double m_init_socp_time = 0;
    double m_sb_computation_time = 0;
    double m_total_socp_time = 0;
    // Iteration related paramters
    int m_num_solves = 0;
    int m_iter_count = 0;
    int m_iter_to_incumbent = 0;
    int m_num_lnodes = 0;
    int m_num_nodes = 0;
    long double m_sum_violation = 0;
    int m_num_fea_sol_found = 0;
    // Aux
    bool m_print_results_on = false;
    bool m_delete_root_node = true;
    bool m_use_fsi = true;

    unordered_map<string, double> m_sequence_lb;
    unordered_map<string, lah_node_record> m_sequence_nodes;
    unordered_map<string, vector<string>> m_sequence_children;
    unordered_map<string, bool> m_sequence_checked;
    vector<node*> m_unexplored_nodes;
};

#endif