/*
 *    BnB_CETSP
 *    main.cpp
 *    Purpose: Solves the CETSP by branch-and-bound
 *
 *    @author Walton Pereira Coutinho
 *    @version 1.0 02/02/2014
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *    
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *    GNU General Public License for more details.
 *    
 *    You should have received a copy of the GNU General Public License 
 *    along with this program. If not, see:
 *    <https://github.com/waltonpcoutinho/BnB_CETSP>.
*/

#include "util.h"
#include "PrintFunctions.h"
#include "Data.h"
#include "BranchNBound.h"
#include "clarabel_interface/SolveRedundantSocpClarabel.h"
#include "structs.h"
#include "CBFS.h"
#include "subgraph.h"
#include "tsp_lb.h"
#include "local_search.h"
#include "cetsp_solver.h"

using namespace std;

long double somaTeste = 0;

int main(int argc, char** argv) 
{
   // Check input 
   if ( argc < 9 || argc > 16) {
      cout << "Wrong calling command\n";
      cout << "./exeCETSP [path da Instancia] [OPTIONS] [OVERLAP FACTOR] [TIME LIMIT] [BRANCHING RULE] [BRANCHING STRATEGY] [ROOT SELECTION] [S.B SIZE]" << endl;	
      exit( 1 );
   }

   char *arqInstancia, *option;
   arqInstancia = argv[ 1 ];
   option = argv[ 2 ];
   double overlap = 0;
   overlap = atof( argv[ 3 ] );
   int timeLimit = atoi( argv[ 4 ] );
   int num_unexpd_nodes_limit = 500000;
   long memory_limit = 8388608;
   long memory_hard_limit = memory_limit;
   bool is_soft_limit_reached = true;
   int branchingRule = 0;
   int selectingRoot = atoi( argv[ 7 ] );
   int strong_branching_size = atoi( argv[ 8 ] );

   if( strcmp( argv[ 5 ], "V1" ) == 0 ) branchingRule = 1;
   if( strcmp( argv[ 5 ], "SB" ) == 0 ) branchingRule = 2;

   int branchingStrategy = 0;
   if( strcmp( argv[ 6 ], "DFS" ) == 0 ) branchingStrategy = 1;
   if( strcmp( argv[ 6 ], "BFS" ) == 0 ) branchingStrategy = 2;
   if( strcmp( argv[ 6 ], "BeFS" ) == 0 ) branchingStrategy = 3;
   if( strcmp( argv[ 6 ], "CBFS" ) == 0 ) branchingStrategy = 4;

   // CBFS variables
   int contour_option = 1;
   int tie_break_option = 3;
   int uncov_cont_option = 1;                   // if 1, the next argv is the bin size, if 2, the next argv is the number of contour
   int uncov_cont_param = 1;                    // size variable related to uncov_cont_option
   int measure_best_mode = 1;                   // if 1, use lower bound, if 2, use best estimate
   int force_branch_rule = 4;                   // represent the variable selection for branching options
   if (argc >= 10) contour_option = atoi(argv[9]);
   // if (argc >= 11) uncov_cont_option = atoi(argv[10]);
   // if (argc >= 12) uncov_cont_param = atoi(argv[11]);
   if (branchingStrategy == 3)
   {
      branchingStrategy = 4;
      contour_option = 3;
   }
   bool new_incumbent_found = false;

   // LAH parameters
   bool use_lah = false;
   if (argc >= 11) use_lah = (atoi(argv[10]) > 0);
   double uncov_lb_time_lim = 5;
   double lah_intersect_tol_init_factor = 1;
   double lah_intersect_tol_min_factor = 100;
   bool curr_lah_not_improving = false;
   // Additional pruning options
   bool use_reverse_prune = false;
   int prev_reverse_prune_flag = 0;
   // Concorde feasible solution improvement search variables
   bool use_fsi = false;
   if (use_lah || argc > 13)
   {
      if (argc >= 12) uncov_lb_time_lim = atof(argv[11]);
      if (argc >= 13) lah_intersect_tol_init_factor = atof(argv[12]);
      if (argc >= 14) lah_intersect_tol_min_factor = atof(argv[13]);
      if (argc >= 15) use_reverse_prune = (atoi(argv[14]) > 0);
      if (argc >= 16) use_fsi = (atoi(argv[15]) > 0);
   }
   else
   {
      if (argc >= 12) use_reverse_prune = (atoi(argv[11]) > 0);
      if (argc >= 13) use_fsi = (atoi(argv[12]) > 0);
   }
   cout<<"Use LAH: "<<use_lah<<endl;
   cout<<"Use Concorde Feasible Solution improvement"<<use_fsi<<endl;
   Data *dataptr = new Data( arqInstancia, option, overlap, argc, argv );

   int sizeInst = dataptr->getSizeInst();
   cout << "Instance Size: " << sizeInst << endl;
   cout << fixed << setiosflags ( ios::showpoint ) << setprecision( 8 );
   double bestKnown = dataptr->getUbFromFile();
   cout << "Initial Upper Bound: " << bestKnown << endl;
   cout << "Average radius size: " << dataptr->m_avg_radius << endl;

   // Additional LAH parameters
   double lah_intersect_tol_init = dataptr->m_avg_radius / lah_intersect_tol_init_factor;
   double lah_intersect_tol_min = dataptr->m_avg_radius / lah_intersect_tol_min_factor;
   lah_intersect_tol_min = max(lah_intersect_tol_min, 0.0001);
   double lah_intersect_tol;
   unordered_map<string, bool> sequence_visit_record;

   // Duplication detect variables
   unordered_map<string, reverse_node_info> reversed_sequence_infos;
   node* rp_temp_node;
   double rp_temp_node_lb;
   bool rp_prev_node_checked = false;
   vector<double> rp_prev_line1_from;
   vector<double> rp_prev_line2_to;

   // Auxilary variables
   int precision = 4;
   int printHeader = 0;
   int optimalFound = 0;
   int print_counter = floor( sizeInst/2 );
   list < node* > open;
   list < node* >:: iterator itOpen;
   // double ub = bestKnown + 1;
   double best = DBL_MAX;
   double best_lb = 0;
   double best_ub = bestKnown;
   unsigned long int count_SOCP_solved = 0;
   double rootLB = 0;
   
   vector< int > solucao;
   vector< vector< double > > solucaoXYZ;
   long int itCount = 0;
   unsigned long int iterCount = 0;
   int itToIncum = 0;
   int insert_id = 0;
   double sbComputationTime = 0;
   // double dbl_compare_constant = 0.000001;

   double initialTotalTimeBnB = cpuTime();
   // Branch and bound function class initialization
   BranchNBound *bnbPtr = new BranchNBound( dataptr );
   // bnbPtr->set_branch_rule(force_branch_rule);
   // ofstream node_info_file ("nodeinfo.out");
   node * root = new node;
   
   //####################################################################
   // Search tree ratio variables
   vector< mpz_class > levels;
   levels.resize( sizeInst );
   mpz_class quantity = 0;
   mpz_t sizeOfTree;
   bnbPtr->computeSizeTree( sizeInst, sizeOfTree, levels );
   mpf_class temp = 0;
   //####################################################################	

   //### root node selection strategy ###
   SolveSocpClarabelStatistics info_struct;
   if( selectingRoot == 1 ) root->pts = bnbPtr->selectRootRedundantClarabel(info_struct);
   if( selectingRoot == 2 ) root->pts = bnbPtr->selectRoot2();
   if( selectingRoot == 3 ) root->pts = bnbPtr->selectRoot3();

   cout << "Initial root: ";
   for ( int i = 0; i <root->pts.size() ; i++ ){
      cout << root->pts[ i ] << " ";
   }
   cout << endl;
   //####################################

   SolveRedundantSocpClarabel *solveSocpPtr = new SolveRedundantSocpClarabel( dataptr, root->pts.size() );
   //solve model
   double totalSocpCompTime = 0;
   double initialSocpCompTime = cpuTime();
   solveSocpPtr->solveSOCP( root->pts );
   totalSocpCompTime += ( cpuTime() - initialSocpCompTime );
   count_SOCP_solved++;

   root->lb = solveSocpPtr->getF_value();
   rootLB = root->lb;
   root->s_lb = 1;
   root->id = 0; root->depth = 1;
   root->intersect_tol = lah_intersect_tol_init;
   solveSocpPtr->printF_value();
   solveSocpPtr->finishSOCP();
   solveSocpPtr->printSolution( root->pts );
   somaTeste += solveSocpPtr->violation;
   
   // check feasibility
   bool feasibilityTest = false;
   vector< double > tempX;
   vector< double > tempY;
   vector< double > tempZ;
   for ( int i = 0; i < root->pts.size(); i++ ){
      tempX.push_back( solveSocpPtr-> getSolutionX( i ) );
      tempY.push_back( solveSocpPtr-> getSolutionY( i ) );
      tempZ.push_back( solveSocpPtr-> getSolutionZ( i ) );
   }
   feasibilityTest = bnbPtr->check_feasibility_Q( root, tempX, tempY, tempZ );
   solveSocpPtr->accumulate_info(info_struct);
   delete solveSocpPtr;

   cout << endl;
   cout << "Not covered in root: ";
   for ( int i = 0; i < bnbPtr->notCoveredBalls.size(); i++ ){
      cout << bnbPtr->notCoveredBalls[ i ] << " ";
   }
   cout << endl;
   //check not covered clients
   for ( int j = 0; j < bnbPtr->notCoveredBalls.size(); j++ ){
      root->notCovered.push_back( bnbPtr->notCoveredBalls[ j ] );
   }
   root->insert_spans.push_back(make_pair(0, root->pts.size()));

   if ( feasibilityTest == true ){
      cout << "FEASIBLE ROOT" << endl;
      best_lb = rootLB;
      best = rootLB;
      double computationTime = cpuTime() - initialTotalTimeBnB;
      double gap_root = ( ( best_ub - rootLB )/ best_ub )*100;
      printDataToFile( dataptr, option, overlap, sizeInst, bestKnown, best, best_lb, 0, count_SOCP_solved, 0, computationTime, 0, computationTime, 0, 0, 0, branchingStrategy );
      
      delete bnbPtr;
      delete dataptr;
      delete root;
      return 0;
   }
   else{
      cout << "INFEASIBLE ROOT" << endl;
   }
   cout << endl;

   tempX.clear();
   tempY.clear();
   tempZ.clear();

   // Additional LB method initialization
   int max_search_depth = -1;
   double temp_uncovered_lb = 0;
   CETSP_solver* aux_cetsp_solve;
   // if (use_lah)
   // {
      // CETSP solver option initialization:
      //    By default, BFS is used for best LB improvement. 
      //    Time limit is set by parameter uncov_lb_time_lim. CBFS options set to the same default values
      //    as the main function. Although with BFS as default search strategy, the only option that matters
      //    is arbitrary tie breaking rule.
      CETSP_Options ctop;
      ctop.overlap_ratio = overlap;
      ctop.branching_rule = 1;
      ctop.branching_strategy = 3;
      ctop.strong_branching_size = 1;
      ctop.time_limit = uncov_lb_time_lim;
      ctop.root_select_rule = selectingRoot;
      ctop.print_on = false;

      ctop.contour_option = 1;
      ctop.tie_break_option = 3;
      ctop.uncov_cont_option = 1;
      ctop.uncov_cont_param = 1;
      ctop.measure_best_mode = 1;

      aux_cetsp_solve = new CETSP_solver(dataptr, ctop);
   // }

   open.push_back( root );
   bnbPtr->add_lb(root->lb);
   best_lb = root->lb;
   
   CBFS* cbfs = new CBFS(dataptr, contour_option, tie_break_option);
   if (contour_option == 2)
   {
      switch (uncov_cont_option)
      {
      case 1:
         cbfs->setContourBinSize(uncov_cont_param);
         break;
      case 2:
         cbfs->setContourNumCont(uncov_cont_param);
         break;
      }
   }
   cbfs->setMeasureBestMode(measure_best_mode);
   cbfs->addNode(root);
   cbfs->resetCurContour();

   //set the size of the strong branching
   int strongBranchingSize = 0;
   if ( branchingRule == 1 ){
      strongBranchingSize = 1;
   }
   if ( branchingRule == 2 ){
      strongBranchingSize = strong_branching_size;
   }

   int sum = 0;
   vector<sbAuxStruct> vectorOfChildren;
   vector< int >::iterator stBrchit;

   // while( !open.empty() && cpuTime() - initialTotalTimeBnB <= timeLimit )
   while(true)
   {
      long temp_memory = get_curr_memory_consumption();
      // Termination conditions check:
      if (!is_soft_limit_reached && temp_memory > memory_limit)
      {
         is_soft_limit_reached = true;
         use_lah = true;
         cout << "Turn on LAH to save space." << endl;
      }
      if (cbfs->m_num_unexplrd_nodes == 0 || best_ub - best_lb < dbl_compare_constant)
      {
         optimalFound = 0;
         break;
      }
      else if (cpuTime() - initialTotalTimeBnB > timeLimit)
      {
         optimalFound = 1;
         break;
      }
      else if (is_soft_limit_reached && temp_memory > memory_hard_limit)
      {
         // terminate if more than memory_hard_limit memory is used
         optimalFound = 1;
         break;
      }

      node * current;
      //######### Depth First Search ###########
      // DFS is currently disabled	
      //######### Depth First Search ###########

      //######### Breadth First Search ###########
      // BrFS is currently disabled
      //######### Breadth First Search ###########

      //######### Best First Search ###########
      // BFS procedure is included in CBFS
      //######### Best First Search ###########

      //######### Cyclic Best First Search ###########
      if (branchingStrategy == 4)
      {
         current = cbfs->getNextNode();
         // open.remove(current);
         bnbPtr->remove_lb(current->lb);
      }
      new_incumbent_found = false;
      //######### Cyclic Best First Search ###########
      
      // print_node_info_to_file(current, node_info_file);
      // print_turn_points_to_file(current, node_info_file);

      if ( current->notCovered.empty() || current->lb >= best_ub - dbl_compare_constant ){
         // FATHOMED BY BOUND
         quantity += mpz_class( sizeOfTree )/levels[ current->pts.size() - 3 ] - 1;
         // cbfs->keep_curr_contour = true;
         if (use_lah)
         {
            string parent_seq = convert_seq_to_string(current->pts);
            if (sequence_visit_record.count(parent_seq) != 0 && sequence_visit_record[parent_seq])
            {
               sequence_visit_record.erase(parent_seq);
            }
         }
      }
      else{
         stBrchit = current->notCovered.begin();
         //control strong branching size
         if ( branchingRule == 2 ){
            strongBranchingSize = strong_branching_size;
            strongBranchingSize -= current->pts.size() - 3;
            strongBranchingSize = max(1, strongBranchingSize);
         }
         
         // Alternative LB calculation
         //============================================================
         // Use high intersection tolerance to look ahead into the search tree rooted
         // at the current node to update LB.
         temp_uncovered_lb = 0;
         curr_lah_not_improving = false;
         if (use_lah && current->add_uncovered_lb)
         // if (use_lah && cbfs->m_from_top_level && current->add_uncovered_lb)
         // if (use_lah && cbfs->m_from_top_level)
         {
            string parent_seq = convert_seq_to_string(current->pts);
            // Prevent the size of the lb dictionary from getting too large
            // if (aux_cetsp_solve->m_sequence_lb.size() > 50000)
            // {
            //    int count = 0;
            //    for (auto node_iter = aux_cetsp_solve->m_sequence_lb.begin(); count < 10000; count++)
            //    {
            //       node_iter = aux_cetsp_solve->m_sequence_lb.erase(node_iter);
            //    }
            //    // aux_cetsp_solve->m_sequence_lb.clear();
            // }
            if (aux_cetsp_solve->m_sequence_nodes.size() > 50000)
            {
               int count = 0;
               for (auto node_iter = aux_cetsp_solve->m_sequence_nodes.begin(); count < 10000; count++)
               {
                  node_iter = aux_cetsp_solve->m_sequence_nodes.erase(node_iter);
               }
               // aux_cetsp_solve->m_sequence_nodes.clear();
            }

            // double temp_time_lim = uncov_lb_time_lim / 
                     // (static_cast<double>(pow(static_cast<double>(current->depth),2)));
            // double temp_time_lim = uncov_lb_time_lim / 
            //          (pow(1.2, static_cast<double>(current->depth-1)));
            double temp_time_lim = uncov_lb_time_lim;

            if (sequence_visit_record.count(parent_seq) != 0 && sequence_visit_record[parent_seq])
            {
               sequence_visit_record.erase(parent_seq);
            }
            else if (temp_time_lim >= 0.05)
            {
               if (best_ub < bestKnown - dbl_compare_constant)
               {
                  aux_cetsp_solve->m_best_known = best_ub;
                  aux_cetsp_solve->m_best_ub = best_ub;
               }   
               lah_intersect_tol = current->intersect_tol;
               aux_cetsp_solve->m_time_limit = temp_time_lim;
               aux_cetsp_solve->set_root_node(current);
               aux_cetsp_solve->set_intersect_tol(0.0001);
               int aux_cetsp_flag = aux_cetsp_solve->solve_keep_history();
               // aux_cetsp_solve->reset_root_node();

               if (aux_cetsp_solve->m_best_ub < best_ub - dbl_compare_constant)
               {
                  best = aux_cetsp_solve->m_best_ub;
                  best_ub = best;
                  solucao = aux_cetsp_solve->m_solution_sequence;
                  itToIncum = itCount;
                  new_incumbent_found = true;
                  // cbfs->clean_up(best_ub);
               }

               if (aux_cetsp_flag == 0)
               {                  
                  iterCount++;
                  quantity += mpz_class(sizeOfTree)/levels[current->pts.size() - 3] - 1;
                  temp = mpz_class( sizeOfTree );
                  temp = ( quantity/temp )*100;
                  cbfs->delNode(current);
                  delete current;
                  continue;
               }
               else if (aux_cetsp_flag == 1)
               {
                  temp_uncovered_lb = aux_cetsp_solve->m_best_lb;
               }
               // cout << "===========================" << endl;
               // cout << "Original LB: " << current->lb << endl;
               // cout << "CETSP LB Algo test. Uncov size: " << current->notCovered.size() 
               //    << "; Intersect tol: " << lah_intersect_tol 
               //    << "; Contour level: " << current->contour 
               //    << "; LAH time limit: " << uncov_lb_time_lim << endl;
               // cout << "flag: " << aux_cetsp_flag 
               //    << "; Cycle: " << cbfs->m_num_cycle
               //    << "; LAH run time: " << aux_cetsp_solve->m_computation_time
               //    << endl;
               // cout << "===========================" << endl;
               aux_cetsp_solve->clean_up();
            }
         }
         //============================================================

         //begin strong branching
         double initialTimeSB = cpuTime();

         if ( strongBranchingSize > 1 )
         {
            cout << "Doing Strong Branching: " << strongBranchingSize << endl;
         }
         for ( int t = 0; t < strongBranchingSize && t < current->notCovered.size(); t++ )
         {
            insert_id = ( *stBrchit );
            stBrchit++;
            // insert_id = bnbPtr->get_branch_node_id(current->notCovered, current->pts);
            sbAuxStruct tempSbStr;		
            tempSbStr.index = insert_id;
            tempSbStr.sum = 0;

            SolveRedundantSocpClarabel *solveSocpPtr2 = new SolveRedundantSocpClarabel(dataptr, current->pts.size()+1);
            solveSocpPtr2->initialize_model();

            for (int span_id = 0; span_id < current->insert_spans.size(); span_id++)
            {
               int start_pos = current->insert_spans[span_id].first;
               int end_pos = current->insert_spans[span_id].second;
               prev_reverse_prune_flag = 0;

               int prev_insert_pos = start_pos + 1;
               int curr_insert_pos = start_pos + 1;
               
               rp_prev_node_checked = false;
               for (int pos = start_pos; pos < end_pos; pos++)
               {
                  node * child = new node;
                  curr_insert_pos = pos + 1;

                  itCount++;
                  child->add_uncovered_lb = current->add_uncovered_lb;
                  // child->add_uncovered_lb = true;
                  child->s_lb = 0;
                  child->depth = current->depth + 1;
                  child->id = itCount;
                  child->intersect_tol = current->intersect_tol;
                  child->pts = current->pts;
                  child->pts.insert(child->pts.begin() + (pos+1), insert_id);

                  double curr_child_socp_lb = 0;
                  double child_temp_lb = 0;
                  bool has_lah_node_hit = false;
                  if (use_lah)
                  {
                     string child_seq = convert_seq_to_string(child->pts);
                     // if (aux_cetsp_solve->m_sequence_lb.count(child_seq) != 0)
                     // {
                     //    // get lah lb and remove corresponding entry from lb dictionary
                     //    child_temp_lb = aux_cetsp_solve->m_sequence_lb[child_seq];
                     //    aux_cetsp_solve->m_sequence_lb.erase(child_seq);
                     //    if (aux_cetsp_solve->m_sequence_nodes.count(child_seq) > 0)
                     //    {
                     //       // if it is in node dictionary, remove
                     //       // corresponding entry
                     //       aux_cetsp_solve->m_sequence_nodes.erase(child_seq);
                     //    }
                     //    sequence_visit_record[child_seq] = true;
                     //    // if the lb is no better than incumbent, the node can be pruned
                     //    if (child_temp_lb >= best_ub - dbl_compare_constant)
                     //    {
                     //       rp_prev_node_checked = false;
                     //       tempSbStr.sum += best_ub;
                     //       delete child;
                     //       continue;
                     //    }
                     // }
                     // if lah node dictionary has corresponding entry, replace current child node with
                     // the one stored (that has socp and feasibility checked)
                     if (aux_cetsp_solve->m_sequence_nodes.count(child_seq) > 0)
                     {
                        // Sanity check to make sure the hit is not a phantom hit
                        // if (aux_cetsp_solve->m_sequence_nodes[child_seq].depth == child->depth)
                        // {
                           // copy info
                           child->lb = aux_cetsp_solve->m_sequence_nodes[child_seq].lb;
                           child->notCovered = aux_cetsp_solve->m_sequence_nodes[child_seq].not_covered;
                           child->feasible = aux_cetsp_solve->m_sequence_nodes[child_seq].feasible;
                           child->not_in_sequence_covered_node_locations = 
                                 aux_cetsp_solve->m_sequence_nodes[child_seq].not_in_sequence_covered_node_locations;

                           sequence_visit_record[child_seq] = true;
                           curr_child_socp_lb = child->lb;
                           // if (child_temp_lb > 0 && child_temp_lb > child->lb)
                           child->lb = max(temp_uncovered_lb, max(child_temp_lb, max(current->lb, child->lb)));

                           aux_cetsp_solve->m_sequence_nodes.erase(child_seq);
                           has_lah_node_hit = true;
                           count_SOCP_solved++;
                        // }
                     }
                  }

                  if (!has_lah_node_hit)
                  {
                     if (use_reverse_prune)
                     {
                        string child_seq = convert_seq_to_string(child->pts);
                        if (reversed_sequence_infos.count(child_seq) > 0)
                        {
                           rp_prev_node_checked = false;
                           reverse_node_info temp_reverse_node_info = reversed_sequence_infos[child_seq];
                           reversed_sequence_infos.erase(child_seq);
                           child->lb = max(temp_reverse_node_info.lb, max(child_temp_lb, max(temp_uncovered_lb, current->lb)));
                           if (child->lb >= best_ub - dbl_compare_constant)
                           {
                              tempSbStr.sum += best_ub;
                              delete child;
                              continue;
                           }
                           else
                           {
                              tempSbStr.sum += child->lb;
                           }
                           child->notCovered = temp_reverse_node_info.notCovered;
                           child->insert_spans = temp_reverse_node_info.insert_spans;
                           child->not_in_sequence_covered_node_locations = 
                              temp_reverse_node_info.not_in_sequence_covered_node_locations;
                           // child->reversible_seg_infos = current->reversible_seg_infos;
                           child->reversible_seg_infos = temp_reverse_node_info.reversible_seg_infos;
                           tempSbStr.candidates.push_back( child );
                           continue;
                        }
                     }
                     
                     initialSocpCompTime = cpuTime();
                     if (solveSocpPtr2->m_num_solves == 0)
                     {
                        solveSocpPtr2->populate_removable_constraints(child->pts);
                     }
                     else
                     {
                        solveSocpPtr2->clear_removable_constraints(prev_insert_pos, curr_insert_pos);
                        solveSocpPtr2->populate_removable_constraints(child->pts, prev_insert_pos, curr_insert_pos);
                     }
                     solveSocpPtr2->solveSOCP();
                     solveSocpPtr2->accumulate_info(info_struct);
                     prev_insert_pos = curr_insert_pos;
                     totalSocpCompTime += ( cpuTime() - initialSocpCompTime );
                     count_SOCP_solved++;
                     curr_child_socp_lb = solveSocpPtr2->getF_value();

                     child->lb = max(curr_child_socp_lb, max(child_temp_lb, max(temp_uncovered_lb, current->lb)));
                     
                     if (child->lb >= best_ub - dbl_compare_constant)
                     {
                        rp_prev_node_checked = false;
                        tempSbStr.sum += best_ub;
                        delete child;
                        continue;
                     }
                     else
                     {
                        tempSbStr.sum += child->lb;
                     }

                     somaTeste += solveSocpPtr2->violation;
                     int temp_size_curr_seq = child->pts.size();
                     tempX.resize(temp_size_curr_seq);
                     tempY.resize(temp_size_curr_seq);
                     tempZ.resize(temp_size_curr_seq);
                     for ( int i2 = 0; i2 < child->pts.size(); i2++ ){
                        tempX[i2] = solveSocpPtr2-> getSolutionX(i2);
                        tempY[i2] = solveSocpPtr2-> getSolutionY(i2);
                        tempZ[i2] = solveSocpPtr2-> getSolutionZ(i2);
                     }

                     //check feasibility
                     bool feasibilityTest;
                     feasibilityTest = bnbPtr->check_feasibility_Q( current, child, pos, insert_id, tempX, tempY, tempZ );
                     child->feasible = (feasibilityTest) ? 1 : 0;
                  }

                  if (child->feasible == 1)
                  {
                     rp_prev_node_checked = false;
                     if (use_fsi && !has_lah_node_hit)
                     {
                        double alter_lb = bnbPtr->find_new_solution_from_curr_sequence(tempX, tempY, tempZ, child);
                        if (alter_lb > 0)
                           child->lb = alter_lb;
                     }

                     if (child->lb < best - dbl_compare_constant)
                     {
                        best = child->lb;
                        solucao = child->pts;
                        solucaoXYZ.clear();
                        solucaoXYZ.push_back(tempX); solucaoXYZ.push_back(tempY); solucaoXYZ.push_back(tempZ);
                        itToIncum = itCount;
                        quantity += mpz_class(sizeOfTree)/levels[child->pts.size() - 3] - 1;
                        if (best < bestKnown - dbl_compare_constant)
                        {
                           best_ub = best;
                           new_incumbent_found = true;
                        }
                     }

                     if (use_lah)
                     {
                        string child_seq = convert_seq_to_string(child->pts);
                        if (sequence_visit_record.count(child_seq) != 0 && sequence_visit_record[child_seq])
                        {
                           sequence_visit_record.erase(child_seq);
                        }
                     }
                     delete child;
                  }
                  else
                  {
                     if (!has_lah_node_hit)
                     {
                        child->notCovered = bnbPtr->notCoveredBalls;
                     }
                     // prev_reverse_prune_flag = 0;
                     child->insert_spans.push_back(make_pair(start_pos, end_pos+1));

                     // RP procedure
                     if (use_reverse_prune && !has_lah_node_hit)
                     {
                        string temp_key;
                        int begin_of_seg_pos = (!rp_prev_node_checked)? pos : pos-1;
                        int end_of_seg_pos = 0;
                        
                        if (!rp_prev_node_checked && pos == end_pos-1)
                        {}
                        else
                        {
                           if (!rp_prev_node_checked)
                           {
                              end_of_seg_pos = (pos == end_pos-2) ? 0 : pos+3;
                           }
                           else
                           {
                              end_of_seg_pos = (pos == end_pos-1) ? 0 : pos+2;
                           }
                           vector<double> line1_from = {solveSocpPtr2->getSolutionX(begin_of_seg_pos), 
                                                         solveSocpPtr2->getSolutionY(begin_of_seg_pos), 
                                                         solveSocpPtr2->getSolutionZ(begin_of_seg_pos)};
                           vector<double> line2_to = {solveSocpPtr2->getSolutionX(end_of_seg_pos), 
                                                      solveSocpPtr2->getSolutionY(end_of_seg_pos), 
                                                      solveSocpPtr2->getSolutionZ(end_of_seg_pos)};
                           if (!rp_prev_node_checked)
                           {
                              rp_temp_node = child;
                              rp_temp_node_lb = curr_child_socp_lb;
                              rp_prev_line1_from = line1_from;
                              rp_prev_line2_to = line2_to;
                              rp_prev_node_checked = true;
                           }
                           else
                           {
                              double diff_sum = Norm_2(difference(rp_prev_line1_from, line1_from)) + Norm_2(difference(rp_prev_line2_to, line2_to));
                              if (diff_sum < 0.01 && rp_temp_node->notCovered == child->notCovered)
                              {
                                 temp_key = to_string(child->pts[begin_of_seg_pos]) + '-' + to_string(child->pts[end_of_seg_pos]);
                                 vector<double> temp_pair_lb_info = {static_cast<double>(child->pts[begin_of_seg_pos+1]), 
                                                                     static_cast<double>(child->pts[begin_of_seg_pos+2]), 
                                                                     curr_child_socp_lb, rp_temp_node_lb};
                                 rp_temp_node->reversible_seg_infos[temp_key] = {line1_from, line2_to, temp_pair_lb_info};
                                 child->reversible_seg_infos = rp_temp_node->reversible_seg_infos;
                              }
                              // Since end_of_seg_pos is 0 only when pos = end_pos-1, after the current iteration, no more child
                              // node will be generated.
                              if(end_of_seg_pos != 0)
                              {
                                 rp_temp_node = child;
                                 rp_temp_node_lb = curr_child_socp_lb;
                                 rp_prev_line1_from = {solveSocpPtr2->getSolutionX(begin_of_seg_pos+1), 
                                                         solveSocpPtr2->getSolutionY(begin_of_seg_pos+1), 
                                                         solveSocpPtr2->getSolutionZ(begin_of_seg_pos+1)};
                                 end_of_seg_pos = (pos == end_pos-2) ? 0 : pos+3;
                                 rp_prev_line2_to = {solveSocpPtr2->getSolutionX(end_of_seg_pos), 
                                                      solveSocpPtr2->getSolutionY(end_of_seg_pos), 
                                                      solveSocpPtr2->getSolutionZ(end_of_seg_pos)};
                              }
                           }
                        }

                        double rp_temp_rev_lb;

                        if (!current->reversible_seg_infos.empty())
                        {
                           vector<pair<int, int>> candidate_reversible_pairs;
                           vector<pair<int, int>> candidate_reversible_endids;
                           vector<double> candidate_reversible_lbs;
                           int num_reverse_count = 0; int begin_of_seg_pos = 0; int end_of_seg_pos = 0;
                           for (int pos_index = 0; pos_index < child->pts.size()-2; pos_index++)
                           {
                              begin_of_seg_pos = pos_index;
                              end_of_seg_pos = (pos_index + 3 == child->pts.size()) ? 0 : pos_index + 3;
                              int begin_of_seg_id = child->pts[begin_of_seg_pos];
                              int end_of_seg_id = child->pts[end_of_seg_pos];
                              temp_key = to_string(begin_of_seg_id) + '-' + to_string(end_of_seg_id);
                              // bool keep_seg = false;
                              if (current->reversible_seg_infos.count(temp_key) > 0 && !current->reversible_seg_infos[temp_key].empty())
                              {
                                 vector<double> prev_line1_from = current->reversible_seg_infos[temp_key][0];
                                 vector<double> prev_line2_to = current->reversible_seg_infos[temp_key][1];
                                 vector<double> curr_line1_from = {solveSocpPtr2->getSolutionX(begin_of_seg_pos), 
                                                                   solveSocpPtr2->getSolutionY(begin_of_seg_pos), 
                                                                   solveSocpPtr2->getSolutionZ(begin_of_seg_pos)};
                                 vector<double> curr_line2_to = {solveSocpPtr2->getSolutionX(end_of_seg_pos), 
                                                                 solveSocpPtr2->getSolutionY(end_of_seg_pos), 
                                                                 solveSocpPtr2->getSolutionZ(end_of_seg_pos)};
                                 double diff_sum = Norm_2(difference(prev_line1_from, curr_line1_from)) + Norm_2(difference(prev_line2_to, curr_line2_to));
                                 if (diff_sum < 0.01)
                                 {
                                    num_reverse_count++;
                                    vector<double> curr_reversible_seq_lb_info = current->reversible_seg_infos[temp_key][2];
                                    if (child->pts[begin_of_seg_pos+1] == static_cast<int>(curr_reversible_seq_lb_info[0]))
                                    {
                                       rp_temp_rev_lb = (child->lb - curr_reversible_seq_lb_info[2]) + curr_reversible_seq_lb_info[3];
                                       curr_reversible_seq_lb_info[2] = child->lb;
                                       curr_reversible_seq_lb_info[3] = rp_temp_rev_lb;
                                    }
                                    else
                                    {
                                       rp_temp_rev_lb = (child->lb - curr_reversible_seq_lb_info[3]) + curr_reversible_seq_lb_info[2]; 
                                       curr_reversible_seq_lb_info[3] = child->lb;
                                       curr_reversible_seq_lb_info[2] = rp_temp_rev_lb;
                                    }
                                    
                                    child->reversible_seg_infos[temp_key] = {prev_line1_from, prev_line2_to, curr_reversible_seq_lb_info};
                                    candidate_reversible_pairs.push_back(make_pair(begin_of_seg_pos, end_of_seg_pos));
                                    // candidate_reversible_endids.push_back(make_pair(begin_of_seg_id, end_of_seg_id));
                                    candidate_reversible_lbs.push_back(rp_temp_rev_lb);
                                 }
                              }
                           }
                           if (num_reverse_count > 0)
                           {
                              for (int rev_index = 0; rev_index < num_reverse_count; rev_index++)
                              {
                                 // Add reverse tour to the hash table so that to avoid an SOCP and feasibility check run.
                                 begin_of_seg_pos = candidate_reversible_pairs[rev_index].first;
                                 vector<int> reversed_sequence = child->pts;
                                 double temp_for_swap = reversed_sequence[begin_of_seg_pos+1];
                                 reversed_sequence[begin_of_seg_pos+1] = reversed_sequence[begin_of_seg_pos+2];
                                 reversed_sequence[begin_of_seg_pos+2] = temp_for_swap;
                                 string reversed_sequence_string = convert_seq_to_string(reversed_sequence);
                                 reversed_sequence_infos[reversed_sequence_string] = 
                                       reverse_node_info(candidate_reversible_lbs[rev_index], child->notCovered, child->insert_spans, 
                                                         child->not_in_sequence_covered_node_locations,
                                                         child->reversible_seg_infos);
                              }
                           }
                        }
                     }
                     child->uncovered_node_min_dist_to_edges.clear();
                     tempSbStr.candidates.push_back( child );     
                  }                  
               }
            }
            bnbPtr->m_uncov_nodes_dists_to_edges.clear();
            vectorOfChildren.push_back( tempSbStr );
            solveSocpPtr2->finishSOCP();
            delete solveSocpPtr2;
         }
         sbComputationTime += cpuTime() - initialTimeSB;

         int pos = 0;
         //take the best sum when there are more than one set of candidates
         if (strongBranchingSize > 1)
         {
            double bestSum = 0;
            for ( int s0 = 0; s0 < vectorOfChildren.size(); s0++ )
            {
               if( vectorOfChildren[ s0 ].sum > bestSum )
               {
                  bestSum = vectorOfChildren[ s0 ].sum;
                  pos = s0;
               }
            }
            // clean unused nodes
            for ( int s0 = 0; s0 < vectorOfChildren.size(); s0++ )
            {
               if( s0 != pos )
               {
                  for ( int s00 = 0; s00 < vectorOfChildren[ s0 ].candidates.size(); s00++ )
                     delete vectorOfChildren[ s0 ].candidates[ s00 ];
                  // delete vectorOfChildren[ s0 ];
               }
            }
         }

         // do branching for the selected vertex
         for ( int s = 0; s < vectorOfChildren[ pos ].candidates.size(); s++ )
         {
            if( vectorOfChildren[ pos ].candidates[ s ]->lb < best_ub - dbl_compare_constant )
            {
               // open.push_back( vectorOfChildren[ pos ].candidates[ s ] );
               bnbPtr->add_lb(vectorOfChildren[pos].candidates[s]->lb);
               cbfs->addNode(vectorOfChildren[pos].candidates[s]);
               //print log
               if ( print_counter == floor( sizeInst/2 ) )
               {
                  double temp_ub = (best_ub < bestKnown) ? best_ub : DBL_MAX;
                  // bnbPtr->printLog(temp, sizeOfTree, vectorOfChildren[ pos ].candidates[ s ], open, count_SOCP_solved, quantity, temp_ub, best_lb, vectorOfChildren[ pos ].candidates[ s ]->notCovered.size(), &printHeader );
                  bnbPtr->printLog(temp, sizeOfTree, vectorOfChildren[ pos ].candidates[ s ], cbfs->m_num_unexplrd_nodes, count_SOCP_solved, quantity, temp_ub, best_lb, vectorOfChildren[ pos ].candidates[ s ]->notCovered.size(), &printHeader );
                  print_counter = 0;
               }
               print_counter++;
            }
            else
            {
               if (use_lah)
               {
                  string temp_seq_str = convert_seq_to_string(vectorOfChildren[pos].candidates[s]->pts);
                  if (sequence_visit_record.count(temp_seq_str) != 0 && sequence_visit_record[temp_seq_str])
                  {
                     sequence_visit_record.erase(temp_seq_str);
                  }
               }
               quantity += mpz_class( sizeOfTree )/levels[ vectorOfChildren[ pos ].candidates[ s ]->pts.size() - 3 ] - 1;
               delete vectorOfChildren[ pos ].candidates[ s ];
            }
         }

         // update global lower bound
         //---------------------------------------------------------------------
         best_lb = bnbPtr->get_current_best_lb();
         //---------------------------------------------------------------------
         vectorOfChildren.clear();
      }
      temp = mpz_class( sizeOfTree );
      temp = ( quantity/temp )*100;

      iterCount++;
      cbfs->delNode(current);
      delete current;
      if(new_incumbent_found)
         cbfs->clean_up(best_ub);
   }

   mpf_class prctComp = mpz_class( sizeOfTree );
   prctComp = count_SOCP_solved/prctComp*100;

   double computationTime = cpuTime() - initialTotalTimeBnB;

   double gap_root = ( ( best_ub - rootLB )/ best_ub )*100;
   double gap_real = ( ( bestKnown - best_lb )/ bestKnown )*100;
   double gap_lb_bnb = ( ( best - best_lb )/ best )*100;

   // int numLNodes = cbfs->getNumNodes(best_ub);
   // int numNodes = itCount;
   int numLNodes = -1;
   int numNodes = cbfs->m_max_num_unexplrd_nodes;

   //Finish Branch and Bound
   best = min(best, bestKnown);
   cout << endl;
   cout << "### Final Log ###" << endl << endl;
   if( optimalFound == 0 ){
      cout << "OPTIMAL SOLUTION FOUND" << endl;
      cout << "Function objective value: " << setiosflags (ios::fixed | ios::showpoint) << setprecision( 15 ) << best << endl;
      cout << "Sum of infeasibilities: " << somaTeste << endl;
      cout << "Sequence: ";
      for ( int i = 0; i < solucao.size(); i++ ){
         cout << solucao[ i ] << " ";
      }
      cout << endl;
      cout << "Solution: \n";
      // for ( int j = 0; j < solucaoXYZ[ 0 ].size(); j++ ){
      // }
      cout << endl;		

      // printDataToFile ( dataptr, option, overlap, sizeInst, bestKnown, best, best_lb, gap_root, 
      //                   count_SOCP_solved, itCount, computationTime, sbComputationTime, totalSocpCompTime, 
      //                   itToIncum, numLNodes, numNodes, branchingStrategy );

      // printDataToMatlab( dataptr, sizeInst, overlap, best, solucao, solucaoXYZ );
   }
   else{

      cout << "NO OPTIMAL SOLUTION FOUND" << endl;
      cout << setiosflags ( ios::showpoint ) << setprecision( 30 );
      cout << "Lower Bound: " << best_lb << endl;
      cout << "Upper Bound: " << best_ub << endl;
      cout << "GAP(LB): " << ( (best_ub - best_lb)/best_ub )*100 << "% " << endl;

      // printDataToFile ( dataptr, option, overlap, sizeInst, bestKnown, best_ub, best_lb, gap_real, gap_lb_bnb, 
      //                   gap_root, count_SOCP_solved, itCount, computationTime, sbComputationTime, totalSocpCompTime, 
      //                   itToIncum, numLNodes, numNodes, branchingStrategy );
   }

   cout << "Number of Resolved Nodes: " << count_SOCP_solved << endl;
   cout << "Pruned tree percentage: " << temp << "%" << endl;
   cout << "Computed Tree Percentage: " << prctComp <<  "%" << endl;
   cout << "Total: " << temp + prctComp << endl;
   cout << "Time total: " << computationTime << endl;
   cout << "Time S.B: " << sbComputationTime << endl;
   cout << "Time SOCP: "<< totalSocpCompTime << endl;
   cout << "Iterations to incumbent: " << itToIncum << endl;

   cout << endl << "#################" << endl;	
   cout << info_struct<<endl;
   
   //  for (auto it = open.begin(); it != open.end(); it++)
   //  {
   //      delete (*it);
   //  }
   //  open.clear();
   cbfs->clean_up();

   // node_info_file.close();
   mpz_clear ( sizeOfTree );
   // delete root;

   delete aux_cetsp_solve;
   delete bnbPtr;
   delete cbfs;
   delete dataptr;
   
   return 0;
}
