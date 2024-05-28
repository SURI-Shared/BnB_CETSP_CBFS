#ifndef BranchNBound_H
#define BranchNBound_H

#include <cmath>
#include <string>
#include <iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include <set>
#include<list>
#include<cstdlib>
#include<ctime>
#include<climits>
#include <cfloat>
#include<algorithm>
#include <unordered_map>
#include <numeric>
#include <map>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>

#include"SolveSocpCplex.h"
#include"local_search.h"
#include"Data.h"
#include"structs.h"
#include"util.h"
#include "clarabel_interface/SolveSocpClarabelWithRecycling.h"

using namespace std;

class BranchNBound{

   public:

      BranchNBound( Data * ); //constructor
      ~BranchNBound();

      vector< int > selectRoot ();
      vector< int > selectRootClarabel ();
      vector< int > selectRootClarabelWithRecycling (SolveSocpClarabelWithRecycling** solver_out_ptr, bool reduced_first_correction);
      vector< int > selectRoot2 ();
      vector< int > selectRoot3 ();
      bool check_feasibility_Q(vector<double>& solX, vector<double>& solY, vector<double>& solZ);
      bool check_feasibility_Q( node* cur_node, vector<double>& solX, vector<double>& solY, vector<double>& solZ );
      bool check_feasibility_Q(node* parent_node, node* child_node, int insert_pos, int insert_node_index, vector<double>& solX, vector<double>& solY, vector<double>& solZ);
      bool check_feasibility_Q_accrt_cvrd(node* parent_node, node* child_node, int insert_pos, int insert_node_index, vector<double>& solX, vector<double>& solY, vector<double>& solZ);
      vector< int > insert ( vector< int >, int, int );
      void computeSizeTree( int sizeInst, mpz_t sizeOfTree, vector< mpz_class > &levels );
      void computeLowerBounds( list< node* > *, node *, double * );
      void printLog( mpf_class, mpz_t, node*, list< node* >, unsigned long int, mpz_class, double, double, int, int * );
      void printLog( mpf_class, mpz_t, node*, unsigned long int, unsigned long int, mpz_class, double, double, int, int * );

      bool checkFeasibility( vector< vector< double > >, vector< int > ); // not using
      void bnb_algorithm(); //not using
      int getNotCoveredBalls( int i ); // not using		
      bool crossRoads( vector< int > &, vector< vector< double > > & ); // not using	
      void setBranchingRuleList2(); // not using
      int strongBranching( list< int >, vector< int > ); // not using
      void makeBranching( branching * ); // not using
      
      void set_branch_rule(int r);

      int check_reverse_prune(node* n, vector<double>& tempX, vector<double>& tempY, vector<double>& tempZ, int t1, int t2);
      void update_insert_spans(vector<pair<int, int>>& new_spans, int t1);

      void add_lb(double lb);
      void remove_lb(double lb);
      double get_current_best_lb() { return m_LB_Count.begin()->first; }
      void clear_lb() { m_LB_Count.clear(); };
      int get_branch_node_id(vector<int>& uncovered, vector<int>& sequence);
      void compute_intersects(int node_id, double& intersect_in, 
               double& intersect_out, vector<double>& line_from, vector<double>& line_to);
      void construct_full_sequence(vector<double>& solX, vector<double>& solY, vector<double>& solZ, 
               unordered_map<int, vector<int>>& covered_node_ids_by_location, 
               vector<int>& org_sequence, vector<int>& full_sequence);
      void construct_full_sequence(vector<double>& solX, vector<double>& solY, vector<double>& solZ, 
               unordered_map<int, vector<int>>& covered_node_ids_by_location, 
               vector<int>& org_sequence, vector<int>& full_sequence,
               vector<vector<double>>& sequence_coords);
      double find_new_solution_from_curr_sequence(vector<double>& solX, vector<double>& solY, vector<double>& solZ,
               node* curr_node);

      vector< int > notCoveredBalls;	
      list < int > branchingRule2;
      unordered_map<int,unordered_map<int,vector<double>>> m_uncov_nodes_dists_to_edges;
      unordered_map<int, vector<int>> m_not_in_seq_cvrd_node_locations;
      double m_intersect_tol = 0.0001;

   private:
      Data *objectOfData;
      SolveSocpCplex *objectOfSolveCplex;

      int variable_select_rule;
      double m_uncov_tour_estimate = 0;
      map<double, int> m_LB_Count;

      //	SET FUNCTIONS

      int sizeOfInstance;	
      vector< int > root;	
      void sortNotCovered( vector< int >& notCovered, vector<double>& solX, vector<double>& solY, vector<double>& solZ );
      void sortNotCovered2( vector< int >& notCovered, vector<double>& solX, vector<double>& solY, vector<double>& solZ );
      void sortNotCovered_new(node* cur_node);
      bool compute_dist_helper(int node_id, int& insert_pos, vector<int>& edge_index_to_check, vector<double>& results, 
                                 vector<double>& solX, vector<double>& solY, vector<double>& solZ);
      bool compute_dist_helper_with_update(int node_id, int& insert_pos, vector<int>& edge_index_to_check, vector<double>& results, vector<int>& sequence,
                                 int start_special_range, int end_special_range, vector<double>& solX, vector<double>& solY, vector<double>& solZ);
};

#endif
