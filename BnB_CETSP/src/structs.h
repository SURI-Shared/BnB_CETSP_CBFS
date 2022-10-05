#ifndef structs_H
#define structs_H

#include <unordered_map>
#include <gmp.h>
#include <gmpxx.h>

#include "Data.h"
//#include "BranchNBound.h"

struct node {
   vector< int > pts;
   // vector< vector < double > > solXYZ;
   vector< int > notCovered;
   double lb;
   int s_lb;                              // 1 if representing current global lower bound
   int feasible = 0;
   //======= WENDA CHANGE =======
   int contour;
   int depth;
   long int id;
   double uncov_est = 0;
   bool add_uncovered_lb = true;
   double intersect_tol = 1;
   // vector<int> full_sequence;
   unordered_map<string, vector<vector<double>>> reversible_seg_infos;
   vector<pair<int,int>> insert_spans;
   unordered_map<int, vector<double>> uncovered_node_min_dist_to_edges;
   unordered_map<int, int> not_in_sequence_covered_node_locations;         // the position in sequence of the node after which the indexed node is covered
   //======= END WENDA CHANGE =======
};

struct sbAuxStruct {
   int index;
   double sum;
   vector< node* > candidates;

};

struct branching {
   Data * dataptr;
   node * current;
   mpz_class quantity;
   mpf_class temp;
   mpz_t sizeOftree;
   vector< mpz_class > * levels;
   vector< int > * solucao;
   vector< vector< double > > * solutionXYZ;
   list< node* > * open;
   int * k;
   int * eliminatedNodesCounter;
   long unsigned int * count_SOCP_solved;
   double * ub;
   double * best_ub;
   double * best_lb;
   double * best;
};

struct CETSP_Options
{
   double overlap_ratio = 1.0;
   int branching_rule = 1;
   int branching_strategy = 4;
   int strong_branching_size = 1;
   int time_limit = 1000000;
   int root_select_rule = 1;
   bool print_on = true;

   int contour_option = 1;
   int tie_break_option = 3;
   int uncov_cont_option = 1;
   int uncov_cont_param = 1;
   int measure_best_mode = 1;
};

struct myCompareStruct
{
   vector< int > vectorOfNotCov;
   vector< pair< int, double > > sortedNotCov;

   myCompareStruct( vector< int > a, vector< pair< int, double > > s )
      : vectorOfNotCov( a ), sortedNotCov( s ){
      }

   bool operator() ( pair< int, double > i, pair< int, double > j ){ 
      return i.second > j.second;
   }
};

struct point {
   double x;
   double y;
   double z;
};

struct reverse_node_info
{
   double lb;
   vector<int> notCovered;
   vector<pair<int,int>> insert_spans;
   unordered_map<int, int> not_in_sequence_covered_node_locations;
   unordered_map<string, vector<vector<double>>> reversible_seg_infos;

   reverse_node_info() {}
   reverse_node_info(double b, vector<int> nc, vector<pair<int,int>> insert, 
                     unordered_map<int, int> not_in_seq) : 
      lb(b), notCovered(nc), insert_spans(insert), 
      not_in_sequence_covered_node_locations(not_in_seq)
   {}
   reverse_node_info(double b, vector<int> nc, vector<pair<int,int>> insert, 
                     unordered_map<int, int> not_in_seq, unordered_map<string, vector<vector<double>>> rs) : 
      lb(b), notCovered(nc), insert_spans(insert), 
      not_in_sequence_covered_node_locations(not_in_seq),
      reversible_seg_infos(rs)
   {}
};

struct lah_node_record
{
   double lb;
   int depth;
   int feasible;
   vector<int> not_covered;
   unordered_map<int, int> not_in_sequence_covered_node_locations;

   lah_node_record() {}
   lah_node_record(double l, int d, int fea, vector<int> nc, 
                  unordered_map<int,int> nisc) : 
      lb(l), depth(d), feasible(fea), not_covered(nc),
      not_in_sequence_covered_node_locations(nisc)
   {}
};

#endif
