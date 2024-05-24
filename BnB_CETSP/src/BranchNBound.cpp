#include"BranchNBound.h"

   BranchNBound::BranchNBound( Data *dataObject )
: objectOfData( dataObject )
{			
   sizeOfInstance = objectOfData->getSizeInst();

   //select variable selection rule
   variable_select_rule = 0;
   for ( int i = 2; i < sizeOfInstance; i++ ){
      if ( objectOfData->getRadius( i ) != objectOfData->getRadius( i - 1 ) ){
         variable_select_rule = 2;
         cout << "Using variable selection rule #2!" << endl;
         break;
      }
      else{
         variable_select_rule = 1;
      }
   }
   if ( variable_select_rule == 1)	{
      cout << "Using variable selection rule #1!" << endl;
   }
}	
// destructor
BranchNBound::~BranchNBound()
{
}

int BranchNBound::getNotCoveredBalls( int i )
{
   return notCoveredBalls[ i ];
}

double R3_2norm( vector< double > element )
{	
   double normOfElement = 0;
   normOfElement = sqrt( pow( element[ 0 ], 2 ) + pow( element[ 1 ], 2 ) + pow( element[ 2 ], 2 ) );
   return normOfElement;
}

bool sortVectorOfPairs( pair< int, double > i, pair< int, double > j )
{
   return ( i.second > j.second );
}

bool sortVectorOfPairsMin(pair<int, double> i, pair<int, double> j)
{
   return i.second < j.second;
}

vector< int > BranchNBound::selectRoot()
{
   //	escolher os elementos que entram na raiz
   double greatestSolution = 0;
   double temp = 0;
   vector < int > tempSequence;
   vector < int > sequence;
   tempSequence.resize( 3 );
   sequence.resize( 3 );

   int sizeInst = objectOfData->getSizeInst();

   sequence[ 0 ] = 0;
   sequence[ 1 ] = objectOfData->getDepotFarthest( 0 );

   tempSequence[ 0 ] = 0;
   tempSequence[ 1 ] = objectOfData->getDepotFarthest( 0 );

   for ( int i = 1; i < sizeInst; i++ ){
      tempSequence[ 2 ] = i;
      SolveSocpCplex *solveCplexSGR = new SolveSocpCplex( objectOfData, tempSequence );
      solveCplexSGR->solveSOCP( tempSequence );
      temp = solveCplexSGR->getF_value();
      if ( temp > greatestSolution ){
         greatestSolution = temp;
         sequence[ 2 ] = i;
      }
      delete solveCplexSGR;
   }

   cout << "Raiz: ";
   for ( int i = 0; i < 3; i++ ){
      cout << sequence[ i ] << " ";
   }
   cout << endl;

   return sequence;

}

//select root as node starting at depot, going to farthest neighborhood, then inserting the node that maximizes the cost
//solves the SOCP to pick the last node using a SolveSocpClarabelWithRecycling that is output via solver_out_ptr, which should have been passed as nullptr
vector< int > BranchNBound::selectRootClarabelWithRecycling(SolveSocpClarabelWithRecycling** solver_out_ptr)
{
   //	escolher os elementos que entram na raiz
   double greatestSolution = 0;
   double temp = 0;
   vector < int > tempSequence;
   vector < int > sequence;
   tempSequence.resize( 3 );
   sequence.resize( 3 );

   int sizeInst = objectOfData->getSizeInst();

   sequence[ 0 ] = 0;
   sequence[ 1 ] = objectOfData->getDepotFarthest( 0 );

   tempSequence[ 0 ] = 0;
   tempSequence[ 1 ] = objectOfData->getDepotFarthest( 0 );

   *solver_out_ptr=new SolveSocpClarabelWithRecycling(objectOfData,3);

   for ( int i = 1; i < sizeInst; i++ ){
      tempSequence[ 2 ] = i;
      (*solver_out_ptr)->solveSOCP(tempSequence);
      temp = (*solver_out_ptr)->getF_value();
      if ( temp > greatestSolution ){
         greatestSolution = temp;
         sequence[ 2 ] = i;
      }
   }

   cout << "Raiz: ";
   for ( int i = 0; i < 3; i++ ){
      cout << sequence[ i ] << " ";
   }
   cout << endl;

   return sequence;

}

vector< int > BranchNBound::selectRoot2()
{
   //escolher os elementos que entram na raiz
   double greatestSolution = 0;
   double temp = 0;
   vector < int > tempSequence;
   vector < int > sequence;
   tempSequence.resize( 3 );
   sequence.resize( 3 );

   int sizeInst = objectOfData->getSizeInst();

   for ( int i = 0; i < sizeInst; i++ ){
      for ( int j = i + 1; j < sizeInst; j++ ){
         for ( int k = j + 1; k < sizeInst; k++ ){	
            tempSequence[ 0 ] = i;
            tempSequence[ 1 ] = j;
            tempSequence[ 2 ] = k;
            SolveSocpCplex *solveCplexSGR = new SolveSocpCplex( objectOfData, tempSequence );
            solveCplexSGR->solveSOCP( tempSequence );
            temp = solveCplexSGR->getF_value();
            if ( temp > greatestSolution ){
               greatestSolution = temp;
               sequence = tempSequence;
            }

            solveCplexSGR->finishSOCP();
            delete solveCplexSGR;

            if( cpuTime() > 1800 ){
               return sequence;
            }	
         }					
      }
   }

   cout << "Raiz: ";
   for ( int i = 0; i < 3; i++ ){
      cout << sequence[ i ] << " ";
   }
   cout << endl;

   return sequence;

}

vector< int > BranchNBound::selectRoot3()
{
   //escolher os elementos que entram na raiz
   double greatestSolution = 0;
   double temp = 0;
   vector < int > tempSequence;
   vector < int > sequence;
   tempSequence.resize( 3 );
   sequence.resize( 3 );
   vector< vector< double > > solutionXYZ;	
   vector< double > tempX;
   vector< double > tempY;
   vector< double > tempZ;

   int notCoveredNumber = INT_MAX;	
   int sizeInst = objectOfData->getSizeInst();

   for ( int i = 0; i < sizeInst; i++ ){
      for ( int j = i + 1; j < sizeInst; j++ ){
         for ( int k = j + 1; k < sizeInst; k++ ){				
            tempSequence[ 0 ] = i;
            tempSequence[ 1 ] = j;
            tempSequence[ 2 ] = k;
            SolveSocpCplex *solveCplexSGR3 = new SolveSocpCplex( objectOfData, tempSequence );
            solveCplexSGR3->solveSOCP( tempSequence );
            //temp = solveCplexSGR->getF_value();

            //get solution	
            for ( int i = 0; i < tempSequence.size(); i++ ){
               tempX.push_back( solveCplexSGR3-> getSolutionX( i ) );
               tempY.push_back( solveCplexSGR3-> getSolutionY( i ) );
               tempZ.push_back( solveCplexSGR3-> getSolutionZ( i ) );
            }		
            // solutionXYZ.push_back( tempX );
            // solutionXYZ.push_back( tempY );
            // solutionXYZ.push_back( tempZ );
            //end get solution

            BranchNBound * BnB = new BranchNBound( objectOfData );
            bool feasibilityTest;
            feasibilityTest = BnB->check_feasibility_Q( tempX, tempY, tempZ );
            //check not covered clients
            if ( BnB->notCoveredBalls.size() < notCoveredNumber ){
               notCoveredNumber = BnB->notCoveredBalls.size();
               sequence = tempSequence;
            }

            solveCplexSGR3->finishSOCP();
            // solutionXYZ.clear();
            tempX.clear();
            tempY.clear();
            tempZ.clear();

            delete solveCplexSGR3;
            delete BnB;

            if( cpuTime() > 1800 ){
               return sequence;
            }	
         }					
      }
   }

   cout << "Raiz: ";
   for ( int i = 0; i < 3; i++ ){
      cout << sequence[ i ] << " ";
   }
   cout << endl;

   return sequence;

}

bool BranchNBound::checkFeasibility( vector< vector< double > > solution, vector< int > sequence )
{
   int sizeSequence = 0;
   long double dist = 0;
   int feasibility = 1;

   sizeSequence = solution[ 0 ].size() - 1;

   vector< int > coveredBalls;		
   coveredBalls.resize( sizeOfInstance );

   vector< double > r;
   r.resize( 3 );
   vector< double > p1v;
   p1v.resize( 3 );
   vector< double > p2v;
   p2v.resize( 3 );
   vector< double > v;
   v.resize( 3 );
   vector< double > cv;
   cv.resize( 3 );
   vector< double > cp1;
   cp1.resize( 3 );
   vector< double > cp2;
   cp2.resize( 3 );

   int index = 0;		
   double mr = 0;
   double norm_cv = 0;
   double test = 0;

   notCoveredBalls.clear();

   for ( int j = 0; j < sizeOfInstance; j++ ){
      for ( int i = 0; i < sizeSequence; i++ ){

         cp1[ 0 ] = solution[ 0 ][ i ] - objectOfData->getCoordx( j );
         cp1[ 1 ] = solution[ 1 ][ i ] - objectOfData->getCoordy( j );

         cp2[ 0 ] = solution[ 0 ][ i + 1 ] - objectOfData->getCoordx( j );
         cp2[ 1 ] = solution[ 1 ][ i + 1 ] - objectOfData->getCoordy( j );

         if ( R3_2norm( cp1 ) <= objectOfData->getRadius( j ) + 0.005*objectOfData->getRadius( j ) || R3_2norm( cp2 ) <= objectOfData->getRadius( j ) + 0.005*objectOfData->getRadius( j ) ){
            coveredBalls[ j ] = 1;
            i = sizeSequence;
         }
         else{				
            r[ 0 ] = solution[ 0 ][ i + 1 ] - solution[ 0 ][ i ];
            r[ 1 ] = solution[ 1 ][ i + 1 ] - solution[ 1 ][ i ];

            mr = r[ 1 ]/r[ 0 ];

            v[ 0 ]	= ( mr*( mr*solution[ 0 ][ i ] + solution[ 1 ][ i ] - objectOfData->getCoordy( j ) + objectOfData->getCoordx( j ) ) )/( 1 - mr*mr );
            v[ 1 ] = objectOfData->getCoordy( j ) + ( 1/mr )*( v[ 0 ] - objectOfData->getCoordx( j ) );

            cv[ 0 ] = v[ 0 ] - objectOfData->getCoordx( j );
            cv[ 1 ] = v[ 1 ] - objectOfData->getCoordy( j );

            norm_cv = 0;
            norm_cv = R3_2norm( cv );

            p1v[ 0 ] = v[ 0 ] - solution[ 0 ][ i ];
            p1v[ 1 ] = v[ 1 ] - solution[ 1 ][ i ];
            p1v[ 2 ] = v[ 2 ] - solution[ 2 ][ i ];

            p2v[ 0 ] = v[ 0 ] - solution[ 0 ][ i + 1 ];
            p2v[ 1 ] = v[ 1 ] - solution[ 1 ][ i + 1 ];
            p2v[ 2 ] = v[ 2 ] - solution[ 2 ][ i + 1 ];

            test = R3_2norm( p1v );

            if ( R3_2norm( p2v ) > test ){
               test = R3_2norm( p2v );
            }

            if ( norm_cv <= objectOfData->getRadius( j ) ){
               if ( test <= R3_2norm( r ) ){
                  coveredBalls[ j ] = 1;
                  i = sizeSequence;
               }
               else{
               }
            }
            else{
            }
         }
      }

      cp1[ 0 ] = solution[ 0 ][ sizeSequence ] - objectOfData->getCoordx( j );
      cp1[ 1 ] = solution[ 1 ][ sizeSequence ] - objectOfData->getCoordy( j );

      cp2[ 0 ] = solution[ 0 ][ 0 ] - objectOfData->getCoordx( j );
      cp2[ 1 ] = solution[ 1 ][ 0 ] - objectOfData->getCoordy( j );

      if ( R3_2norm( cp1 ) <= objectOfData->getRadius( j ) + 0.1*objectOfData->getRadius( j ) || R3_2norm( cp2 ) <= objectOfData->getRadius( j ) + 0.1*objectOfData->getRadius( j ) ){
         coveredBalls[ j ] = 1;
      }				
      else{

         r[ 0 ] = solution[ 0 ][ 0 ] - solution[ 0 ][ sizeSequence ];
         r[ 1 ] = solution[ 1 ][ 0 ] - solution[ 1 ][ sizeSequence ];
         r[ 2 ] = solution[ 2 ][ 0 ] - solution[ 2 ][ sizeSequence ];

         mr = r[ 1 ]/r[ 0 ];

         v[ 0 ]	= ( mr*( mr*solution[ 0 ][ sizeSequence ] + solution[ 1 ][ sizeSequence ] - objectOfData->getCoordy( j ) + objectOfData->getCoordx( j ) ) )/( 1 - mr*mr );
         v[ 1 ] = objectOfData->getCoordy( j ) + ( 1/mr )*( v[ 0 ] - objectOfData->getCoordx( j ) );

         cv[ 0 ] = v[ 0 ] - objectOfData->getCoordx( j );
         cv[ 1 ] = v[ 1 ] - objectOfData->getCoordy( j );

         norm_cv = R3_2norm( cv );

         p1v[ 0 ] = v[ 0 ] - solution[ 0 ][ sizeSequence ];
         p1v[ 1 ] = v[ 1 ] - solution[ 1 ][ sizeSequence ];
         p1v[ 2 ] = v[ 2 ] - solution[ 2 ][ sizeSequence ];

         p2v[ 0 ] = v[ 0 ] - solution[ 0 ][ 0 ];
         p2v[ 1 ] = v[ 1 ] - solution[ 1 ][ 0 ];
         p2v[ 2 ] = v[ 2 ] - solution[ 2 ][ 0 ];

         test = R3_2norm( p1v );

         if ( R3_2norm( p2v ) > test ){
            test = R3_2norm( p2v );
         }

         if ( norm_cv <= objectOfData->getRadius( j ) ){
            if ( test <= R3_2norm( r ) ){
               coveredBalls[ j ] = 1;
            }
            else{
            }
         }
         else{
         }
      }
   }

   for ( int i = 0; i < sizeOfInstance; i++ ){
      feasibility *= coveredBalls[ i ];
   }

   for ( int i = 0; i < coveredBalls.size(); i++ ){
      if ( coveredBalls[ i ] != 1 ){
         notCoveredBalls.push_back( i );
      }
   }

   if ( feasibility == 1 ){
      return true;
   }
   else{
      return false;
   }
}

bool BranchNBound::check_feasibility_Q( vector<double>& solX, vector<double>& solY, vector<double>& solZ )
{
   vector< double > c( 3 );		
   vector< double > p1( 3 );		
   vector< double > p2( 3 );
   vector< double > point( 3 );

   double constant = 0.0001;
   long double theta = 0;
   double test;

   solX.push_back(solX.front());
   solY.push_back(solY.front());
   solZ.push_back(solZ.front());

   vector< int > coveredBalls( sizeOfInstance );

   for ( int i = 0; i < sizeOfInstance; i++ ){
      for ( int j = 0; j < solX.size() - 1; j++ ){

         c[ 0 ] = objectOfData->getCoordx( i );
         c[ 1 ] = objectOfData->getCoordy( i );
         c[ 2 ]= objectOfData->getCoordz( i );

         p1[ 0 ] = solX[ j ];
         p1[ 1 ] = solY[ j ];
         p1[ 2 ] = solZ[ j ];

         p2[ 0 ] = solX[ j + 1 ];
         p2[ 1 ] = solY[ j + 1 ];
         p2[ 2 ] = solZ[ j + 1 ];

         theta = - ( dot_product( difference( p2, p1 ), difference( c, p2 ) ) )/( pow ( Norm_2( difference( p2, p1 ) ), 2 ) );
         if ( theta >= 1 )
         {
            theta = 1;
         }
         else if (theta < 0)
         {
            theta = 0;
         }
         point = sum_vector ( scalar_product( theta, p1 ), scalar_product( 1 - theta, p2 ) );
         test = Norm_2( difference( point, c ) );
         if ( test <  objectOfData->getRadius( i ) + constant )
         {
            coveredBalls[ i ] = 1;
         }
      }
   }

   double feasibility = 1;

   for ( int i = 0; i < sizeOfInstance; i++ ){
      feasibility *= coveredBalls[ i ];
   }

   notCoveredBalls.clear();

   if( variable_select_rule == 1 ){
      sortNotCovered( notCoveredBalls, solX, solY, solZ );
   }
   if( variable_select_rule == 2 ){
      sortNotCovered2( notCoveredBalls, solX, solY, solZ );
   }

   if ( feasibility == 1 ){
      return true;
   }
   else{
      return false;
   }
}

// Should only used at root node. Initialize the uncovered node distances to edges from scratch.
bool BranchNBound::check_feasibility_Q( node* cur_node, vector<double>& solX, vector<double>& solY, vector<double>& solZ )
{
   vector< double > c( 3 );		
   vector< double > p1( 3 );		
   vector< double > p2( 3 );
   vector< double > point( 3 );

   double constant = 0.0001;
   long double theta = 0;
   double test;
   int temp_insert_pos;
   int node_id;
   bool is_covered_flag = false;
   bool is_feasible_flag = false;
   int sequence_length = cur_node->pts.size();

   solX.push_back(solX.front());
   solY.push_back(solY.front());
   solZ.push_back(solZ.front());

   vector<int> nodes_to_compute_dists;
   vector<int> edges_to_compute_dists_from;
   vector<double> temp_results(2);

   cur_node->uncovered_node_min_dist_to_edges.clear();
   cur_node->not_in_sequence_covered_node_locations.clear();

   // Populate the list of nodes to find dists and the edges to check
   set<int> nodes_in_sequence;
   for (int i = 0; i < sequence_length; i++)
   {
      nodes_in_sequence.insert(cur_node->pts[i]);
      edges_to_compute_dists_from.push_back(i);
   }
   for (int i = 0; i < sizeOfInstance; i++)
   {
      if (nodes_in_sequence.count(i) == 0)
         nodes_to_compute_dists.push_back(i);
   }

   for ( int i = 0; i < nodes_to_compute_dists.size(); i++ )
   {
      node_id = nodes_to_compute_dists[i];
      double min = DBL_MAX;
      cur_node->uncovered_node_min_dist_to_edges[node_id].resize(2, 0);
      is_covered_flag = compute_dist_helper(node_id, temp_insert_pos, edges_to_compute_dists_from, temp_results, solX, solY, solZ);
      if (is_covered_flag)
      {
         cur_node->not_in_sequence_covered_node_locations[node_id] = temp_insert_pos;
      }
      else
      {
         cur_node->uncovered_node_min_dist_to_edges[node_id] = temp_results;
      }
   }

   if ( cur_node->uncovered_node_min_dist_to_edges.empty() )
      is_feasible_flag = true;
   else
      sortNotCovered_new(cur_node);
   return is_feasible_flag;
}

bool BranchNBound::check_feasibility_Q(node* parent_node, node* child_node, int insert_pos, int insert_node_index, vector<double>& solX, vector<double>& solY, vector<double>& solZ)
{
   vector< double > c( 3 );		
   vector< double > p1( 3 );		
   vector< double > p2( 3 );
   vector< double > point( 3 );

   double constant = 0.0001;
   long double theta = 0;
   double test;
   int temp_insert_pos = 0;
   int temp_insert_pos_2 = 0;
   int sequence_length = child_node->pts.size();
   bool is_covered_flag = false;
   int node_id;
   bool is_feasible_flag = false;
   int start_check_pos = max(0, insert_pos-1);
   int end_check_pos = min(sequence_length, insert_pos+3);

   solX.push_back(solX.front());
   solY.push_back(solY.front());
   solZ.push_back(solZ.front());

   vector<int> edges_to_compute_dists_from;
   vector<double> temp_results(2);
   // DIfferent size of the edges to check may affect performance: 
   //    the ones that are most likely to have substantial changes are the ones
   //    affected by the inserted node
   // for (int i = max(0, insert_pos - max(static_cast<int>(sequence_length*0.8),2)) ; i < min(sequence_length,insert_pos + max(static_cast<int>(sequence_length*0.8),2)+2); i++)
   for (int i = 0; i < sequence_length; i++)
   {
      edges_to_compute_dists_from.push_back(i);
   }

   // Update lookup table for min distance to edges for uncovered nodes
   for ( auto it = parent_node->notCovered.begin(); it != parent_node->notCovered.end(); it++ )
   {
      node_id = (*it);
      if (node_id == insert_node_index) continue;

      is_covered_flag = compute_dist_helper(node_id, temp_insert_pos, edges_to_compute_dists_from, temp_results, solX, solY, solZ);
      // is_covered_flag = compute_dist_helper_with_update(node_id, temp_insert_pos, edges_to_compute_dists_from, temp_results, child_node->pts,
      //                                           start_check_pos, end_check_pos, solX, solY, solZ);

      if (is_covered_flag)
      {
         child_node->not_in_sequence_covered_node_locations[node_id] = temp_insert_pos;
         continue;
      }

      // Vector parent_node->uncovered_node_min_dist_to_edges[node_id] only holds the min value to edge value for sequence in parent node
      // therefore, this code does not work directly.
      // TODO: for this method to work, either change the Vector to hold dists to all edges so that we can compare them individually in current iteration.
      //       Or fix the range of edges to check in each iteration and then keep track of the minimum value outside the range as well as the overall minimum
      // if (temp_results[0] < parent_node->uncovered_node_min_dist_to_edges[node_id][0])
         child_node->uncovered_node_min_dist_to_edges[node_id] = temp_results;
      // else
      //    child_node->uncovered_node_min_dist_to_edges[node_id] = parent_node->uncovered_node_min_dist_to_edges[node_id];
   }

   // Check if covered nodes that are not in sequence are no longer covered
   vector<int> sequence_edges(sequence_length, 0);
   iota(sequence_edges.begin(), sequence_edges.end(), 0);

   for (auto it = parent_node->not_in_sequence_covered_node_locations.begin();
            it != parent_node->not_in_sequence_covered_node_locations.end(); it++)
   {
      node_id = it->first;
      if (it->second >= start_check_pos && it->second < end_check_pos)
      {
         edges_to_compute_dists_from.clear();
         edges_to_compute_dists_from.resize(end_check_pos - start_check_pos);
         iota(edges_to_compute_dists_from.begin(), edges_to_compute_dists_from.end(), start_check_pos);
         // edges_to_compute_dists_from = {insert_pos, insert_pos+1};
      }
      else
      {
         if (it->second < insert_pos)
            child_node->not_in_sequence_covered_node_locations[node_id] = it->second;
         else
            child_node->not_in_sequence_covered_node_locations[node_id] = it->second+1;
         continue;
      }

      is_covered_flag = compute_dist_helper(node_id, temp_insert_pos, edges_to_compute_dists_from, temp_results, solX, solY, solZ);
      if (!is_covered_flag)
      {
         is_covered_flag = compute_dist_helper(node_id, temp_insert_pos_2, sequence_edges, temp_results, solX, solY, solZ);
         if (is_covered_flag)
         {
            child_node->not_in_sequence_covered_node_locations[node_id] = temp_insert_pos_2;
         }
         else
         {
            child_node->uncovered_node_min_dist_to_edges[node_id] = temp_results;
         }
      }
      else
      {
         child_node->not_in_sequence_covered_node_locations[node_id] = temp_insert_pos;
      }
   }

   // When there are no uncovered node indicated, check needs to be performed for each node
   // covered but not in sequence.
   if ( child_node->uncovered_node_min_dist_to_edges.empty() )
   {
      is_feasible_flag = true;
      vector<int> all_nodes;
      m_not_in_seq_cvrd_node_locations.clear();
      for (auto it = child_node->not_in_sequence_covered_node_locations.begin();
            it != child_node->not_in_sequence_covered_node_locations.end(); it++)
      {
         all_nodes.push_back(it->first);
      }
      for (int i = 0; i < all_nodes.size(); i++)
      {
         node_id = all_nodes[i];
         is_covered_flag = compute_dist_helper(node_id, temp_insert_pos, sequence_edges, temp_results, solX, solY, solZ);
         if (!is_covered_flag)
         {
            child_node->not_in_sequence_covered_node_locations.erase(node_id);
            child_node->uncovered_node_min_dist_to_edges[node_id] = temp_results;
            is_feasible_flag = false;
         }
         else
         {
            child_node->not_in_sequence_covered_node_locations[node_id] = temp_insert_pos;
            m_not_in_seq_cvrd_node_locations[temp_insert_pos].push_back(node_id);
         }
      }
   }

   if (!is_feasible_flag)
      sortNotCovered_new(child_node);

   return is_feasible_flag;
}

double BranchNBound::find_new_solution_from_curr_sequence(vector<double>& solX, vector<double>& solY, vector<double>& solZ,
               node* curr_node)
{
   vector<int> full_sequence;
   vector<vector<double>> full_sequence_coords;
   construct_full_sequence(solX, solY, solZ, m_not_in_seq_cvrd_node_locations, 
                           curr_node->pts, full_sequence, full_sequence_coords);
   Local_search_cetsp* alter_solver = new Local_search_cetsp(objectOfData);
   alter_solver->initialize_sequence(full_sequence, full_sequence_coords[0],
                                       full_sequence_coords[1], full_sequence_coords[2]);
   alter_solver->m_init_solution = curr_node->lb;
   double new_sol = alter_solver->solve_new_solution();
   
   delete alter_solver;
   return new_sol;
}

bool BranchNBound::check_feasibility_Q_accrt_cvrd(node* parent_node, node* child_node, int insert_pos, int insert_node_index, vector<double>& solX, vector<double>& solY, vector<double>& solZ)
{
   vector< double > c( 3 );		
   vector< double > p1( 3 );		
   vector< double > p2( 3 );
   vector< double > point( 3 );

   double constant = 0.0001;
   long double theta = 0;
   double test;
   int temp_insert_pos = 0;
   int temp_insert_pos_2 = 0;
   int sequence_length = child_node->pts.size();
   bool is_covered_flag = false;
   int node_id;
   bool is_feasible_flag = false;
   int start_check_pos = max(0, insert_pos-1);
   int end_check_pos = min(sequence_length, insert_pos+3);

   solX.push_back(solX.front());
   solY.push_back(solY.front());
   solZ.push_back(solZ.front());

   vector<int> edges_to_compute_dists_from;
   vector<double> temp_results(2);

   m_not_in_seq_cvrd_node_locations.clear();
   // DIfferent size of the edges to check may affect performance: 
   //    the ones that are most likely to have substantial changes are the ones
   //    affected by the inserted node
   for (int i = 0; i < sequence_length; i++)
   {
      edges_to_compute_dists_from.push_back(i);
   }

   // Update lookup table for min distance to edges for uncovered nodes
   for ( auto it = parent_node->notCovered.begin(); it != parent_node->notCovered.end(); it++ )
   {
      node_id = (*it);
      if (node_id == insert_node_index) continue;

      is_covered_flag = compute_dist_helper(node_id, temp_insert_pos, edges_to_compute_dists_from, temp_results, solX, solY, solZ);

      if (is_covered_flag)
      {
         child_node->not_in_sequence_covered_node_locations[node_id] = temp_insert_pos;
         m_not_in_seq_cvrd_node_locations[temp_insert_pos].push_back(node_id);
         continue;
      }
      child_node->uncovered_node_min_dist_to_edges[node_id] = temp_results;
   }

   // Check if covered nodes that are not in sequence are no longer covered
   vector<int> sequence_edges(sequence_length, 0);
   iota(sequence_edges.begin(), sequence_edges.end(), 0);

   for (auto it = parent_node->not_in_sequence_covered_node_locations.begin();
            it != parent_node->not_in_sequence_covered_node_locations.end(); it++)
   {
      node_id = it->first;
      is_covered_flag = compute_dist_helper(node_id, temp_insert_pos_2, sequence_edges, temp_results, solX, solY, solZ);
      if (is_covered_flag)
      {
         child_node->not_in_sequence_covered_node_locations[node_id] = temp_insert_pos_2;
         m_not_in_seq_cvrd_node_locations[temp_insert_pos_2].push_back(node_id);
      }
      else
      {
         child_node->uncovered_node_min_dist_to_edges[node_id] = temp_results;
      }
   }

   // When there are no uncovered node indicated, check needs to be performed for each node
   // covered but not in sequence.
   if (child_node->uncovered_node_min_dist_to_edges.empty())
   {
      is_feasible_flag = true;
   }

   if (!is_feasible_flag)
   {
      sortNotCovered_new(child_node);
      // construct_full_sequence(solX, solY, solZ, m_not_in_seq_cvrd_node_locations, child_node->pts, child_node->full_sequence);
   }
   return is_feasible_flag;
}

void BranchNBound::sortNotCovered_new(node* cur_node)
{
   vector<pair<int, double>> aux;
   pair<int, double> temp;
   notCoveredBalls.clear();
   for (auto it = cur_node->uncovered_node_min_dist_to_edges.begin(); it != cur_node->uncovered_node_min_dist_to_edges.end(); it++)
   {
      if (variable_select_rule == 1)
         aux.push_back(make_pair(it->first, it->second[0]));
      else if (variable_select_rule == 2)
         aux.push_back(make_pair(it->first, it->second[1]));
   }
   sort( aux.begin(), aux.end(), sortVectorOfPairs );
   for ( int i = 0; i < aux.size(); i++ ){
      notCoveredBalls.push_back(aux[i].first);
   }
}

bool BranchNBound::compute_dist_helper(int node_id, int& insert_pos, vector<int>& edge_index_to_check, vector<double>& results,
                                          vector<double>& solX, vector<double>& solY, vector<double>& solZ)
{
   vector< double > c( 3 );		
   vector< double > p1( 3 );		
   vector< double > p2( 3 );
   vector< double > point( 3 );
   vector< double > point_in_the_border( 3 );

   // double constant = 0.0001;
   double constant = m_intersect_tol;
   long double theta = 0;
   double test;
   bool check_intersect = false;
   int edge_size = edge_index_to_check.size();
   double cur_min = DBL_MAX;
   double sum_dist = 0;

   results.clear();
   results.resize(2, 0);

   for (int i = 0; i < edge_size; i++)
   {
      int edge_id = edge_index_to_check[i];
      c[ 0 ] = objectOfData->getCoordx( node_id );
      c[ 1 ] = objectOfData->getCoordy( node_id );
      c[ 2 ]= objectOfData->getCoordz( node_id );

      p1[ 0 ] = solX[ edge_id ];
      p1[ 1 ] = solY[ edge_id ];
      p1[ 2 ] = solZ[ edge_id ];

      p2[ 0 ] = solX[ edge_id + 1 ];
      p2[ 1 ] = solY[ edge_id + 1 ];
      p2[ 2 ] = solZ[ edge_id + 1 ];

      // double length_p2p1 = Norm_2(difference(p1, p2));
      // double dot_product_p2p1_cp1 = dot_product(difference(p1, p2), difference(p1, c));
      // theta = (dot_product_p2p1_cp1/pow(length_p2p1,2));
      // theta = max(0.0, min(static_cast<double>(theta), 1.0));
      // point = sum_vector(scalar_product(1-theta, p1), scalar_product(theta, p2));
      theta = ( dot_product( difference( p2, p1 ), difference( p2, c ) ) )/( pow ( Norm_2( difference( p2, p1 ) ), 2 ) );
      theta = max(0.0, min(static_cast<double>(theta), 1.0));
      point = sum_vector ( scalar_product( theta, p1 ), scalar_product( 1 - theta, p2 ) );

      test = Norm_2( difference( point, c ) );
      if( test <  objectOfData->getRadius( node_id ) + constant )
      {
         check_intersect = true;
         insert_pos = edge_id;
         break;
      }

      // sum_dist += test - objectOfData->getRadius( node_id );
      if (test < cur_min)
      {
         cur_min = test;
         if (variable_select_rule == 1)
         {
            results[0] = test;
         }
         else if (variable_select_rule == 2)
         {
            // double lambda = objectOfData->getRadius( node_id )/test;
            // point_in_the_border = sum_vector ( scalar_product( lambda, point ), scalar_product( 1 - lambda, c ) );
            // double side_one = Norm_2( difference( point_in_the_border, p1 ) );
            // double side_two = Norm_2( difference( point_in_the_border, p2 ) );
            // double base = Norm_2( difference( p2, p1 ) );
            results[0] = test;
            // results[1] = side_one + side_two - base;
            results[1] = test - objectOfData->getRadius( node_id );
         }
      }
   }
   // if (variable_select_rule == 2)
   //    results[1] = sum_dist;
   return check_intersect;
}

bool BranchNBound::compute_dist_helper_with_update(int node_id, int& insert_pos, vector<int>& edge_index_to_check, vector<double>& results, vector<int>& sequence,
                                                int start_special_range, int end_special_range, vector<double>& solX, vector<double>& solY, vector<double>& solZ)
{
   vector< double > c( 3 );		
   vector< double > p1( 3 );		
   vector< double > p2( 3 );
   vector< double > point( 3 );
   vector< double > point_in_the_border( 3 );

   double constant = 0.0001;
   long double theta = 0;
   double test;
   bool check_intersect = false;
   int edge_size = edge_index_to_check.size();
   double cur_min = DBL_MAX;
   double special_branch_res = 0;

   results.clear();
   results.resize(2, 0);

   for (int i = 0; i < edge_size; i++)
   {
      special_branch_res = 0;
      int edge_id = edge_index_to_check[i];
      int edge_from_node_index = sequence[edge_id];
      if (m_uncov_nodes_dists_to_edges.count(node_id) == 0 || m_uncov_nodes_dists_to_edges[node_id].count(edge_from_node_index) == 0
                                                           || (edge_id >= start_special_range && edge_id < end_special_range))
      {
         c[ 0 ] = objectOfData->getCoordx( node_id );
         c[ 1 ] = objectOfData->getCoordy( node_id );
         c[ 2 ]= objectOfData->getCoordz( node_id );

         p1[ 0 ] = solX[ edge_id ];
         p1[ 1 ] = solY[ edge_id ];
         p1[ 2 ] = solZ[ edge_id ];

         p2[ 0 ] = solX[ edge_id + 1 ];
         p2[ 1 ] = solY[ edge_id + 1 ];
         p2[ 2 ] = solZ[ edge_id + 1 ];

         theta = ( dot_product( difference( p2, p1 ), difference( p2, c ) ) )/( pow ( Norm_2( difference( p2, p1 ) ), 2 ) );
         theta = max(0.0, min(static_cast<double>(theta), 1.0));

         point = sum_vector ( scalar_product( theta, p1 ), scalar_product( 1 - theta, p2 ) );
         test = Norm_2( difference( point, c ) );

         if (variable_select_rule == 2)
         {
            double lambda = objectOfData->getRadius( node_id )/test;
            point_in_the_border = sum_vector ( scalar_product( lambda, point ), scalar_product( 1 - lambda, c ) );
            double side_one = Norm_2( difference( point_in_the_border, p1 ) );
            double side_two = Norm_2( difference( point_in_the_border, p2 ) );
            double base = Norm_2( difference( p2, p1 ) );
            special_branch_res = side_one + side_two - base;
         }

         if (edge_id < start_special_range || edge_id >= end_special_range)
         {
            m_uncov_nodes_dists_to_edges[node_id][edge_from_node_index] = {test, special_branch_res};
         }
      }
      else
      {
         test = m_uncov_nodes_dists_to_edges[node_id][edge_from_node_index][0];
         special_branch_res = m_uncov_nodes_dists_to_edges[node_id][edge_from_node_index][1];
      }

      if( test <  objectOfData->getRadius( node_id ) + constant )
      {
         check_intersect = true;
         insert_pos = edge_id;
         break;
      }

      if (test < cur_min)
      {
         cur_min = test;
         results = {test, special_branch_res};
      }
   }
   return check_intersect;
}

void BranchNBound::sortNotCovered( vector< int >& notCovered, vector<double>& solX, vector<double>& solY, vector<double>& solZ )
{
   vector< pair< int, double > > aux;

   pair< int, double > temp;

   vector< double > c;		
   vector< double > p1;		
   vector< double > p2;
   vector< double > point;

   c.resize( 3 );
   p1.resize( 3 );
   p2.resize( 3 );

   double min = DBL_MAX;
   long double theta = 0;
   m_uncov_tour_estimate = 0;

   for ( int i = 0; i < notCovered.size(); i++ ){
      for ( int j = 0; j < solX.size() - 1; j++ ){

         c[ 0 ] = objectOfData->getCoordx( notCovered[ i ] );
         c[ 1 ] = objectOfData->getCoordy( notCovered[ i ] );
         c[ 2 ]= objectOfData->getCoordz( notCovered[ i ] );

         p1[ 0 ] = solX[ j ];
         p1[ 1 ] = solY[ j ];
         p1[ 2 ] = solZ[ j ];

         p2[ 0 ] = solX[ j + 1 ];
         p2[ 1 ] = solY[ j + 1 ];
         p2[ 2 ] = solZ[ j + 1 ];

         theta = - ( dot_product( difference( p2, p1 ), difference( c, p2 ) ) )/( pow ( Norm_2( difference( p2, p1 ) ), 2 ) );
         if ( theta >= 1 ){
            theta = 1;
            point = sum_vector ( scalar_product( theta, p1 ), scalar_product( 1 - theta, p2 ) );
            double test = Norm_2( difference( point, c ) );	
            if( test < min ){
               min = test;
               temp = make_pair( notCovered[ i ], min );
            }
         }
         else{
            if ( theta > 0 && theta < 1 ){
               point = sum_vector ( scalar_product( theta, p1 ), scalar_product( 1 - theta, p2 ) );
               double test = Norm_2( difference( point, c ) );
               if( test < min ){
                  min = test;
                  temp = make_pair( notCovered[ i ], min );
               }
            }
            else{
               if ( theta <= 0){
                  theta = 0;
                  point = sum_vector ( scalar_product( theta, p1 ), scalar_product( 1 - theta, p2 ) );
                  double test = Norm_2( difference( point, c ) );
                  if( test < min ){
                     min = test;
                     temp = make_pair( notCovered[ i ], min );
                  }
               }
            }
         }
      }

      min = DBL_MAX;
      aux.push_back( temp );
   }

   sort( aux.begin(), aux.end(), sortVectorOfPairs );

   for ( int i = 0; i < notCovered.size(); i++ ){
      notCoveredBalls[ i ] = aux[ i ].first;
   }
}

void BranchNBound::sortNotCovered2( vector< int >& notCovered, vector<double>& solX, vector<double>& solY, vector<double>& solZ )
{
   vector< pair< int, double > > aux;

   pair< int, double > temp;

   vector< double > c( 3 );		
   vector< double > p1( 3 );		
   vector< double > p2( 3 );
   vector< double > point( 3 );
   vector< double > point_in_the_border( 3 );

   double min = DBL_MAX;
   long double theta = 0;
   m_uncov_tour_estimate = 0;

   for ( int i = 0; i < notCovered.size(); i++ ){
      for ( int j = 0; j < solX.size() - 1; j++ ){

         c[ 0 ] = objectOfData->getCoordx( notCovered[ i ] );
         c[ 1 ] = objectOfData->getCoordy( notCovered[ i ] );
         c[ 2 ]= objectOfData->getCoordz( notCovered[ i ] );

         p1[ 0 ] = solX[ j ];
         p1[ 1 ] = solY[ j ];
         p1[ 2 ] = solZ[ j ];

         p2[ 0 ] = solX[ j + 1 ];
         p2[ 1 ] = solY[ j + 1 ];
         p2[ 2 ] = solZ[ j + 1 ];

         theta = - ( dot_product( difference( p2, p1 ), difference( c, p2 ) ) )/( pow ( Norm_2( difference( p2, p1 ) ), 2 ) );
         if ( theta >= 1 ){
            theta = 1;
            point = sum_vector ( scalar_product( theta, p1 ), scalar_product( 1 - theta, p2 ) );
            double test = Norm_2( difference( point, c ) );
            if( test < min ){
               double lambda = objectOfData->getRadius( notCovered[ i ] )/test;
               point_in_the_border = sum_vector ( scalar_product( lambda, point ), scalar_product( 1 - lambda, c ) );
               double side_one = Norm_2( difference( point_in_the_border, p1 ) );
               double side_two = Norm_2( difference( point_in_the_border, p2 ) );
               double base = Norm_2( difference( p2, p1 ) );
               double gamma = side_one + side_two - base;
               min = test;
               temp = make_pair( notCovered[ i ], gamma );
            }
         }
         else{
            if ( theta > 0 && theta < 1 ){
               point = sum_vector ( scalar_product( theta, p1 ), scalar_product( 1 - theta, p2 ) );
               double test = Norm_2( difference( point, c ) );
               if( test < min ){
                  double lambda = objectOfData->getRadius( notCovered[ i ] )/test;
                  point_in_the_border = sum_vector ( scalar_product( lambda, point ), scalar_product( 1 - lambda, c ) );
                  double side_one = Norm_2( difference( point_in_the_border, p1 ) );
                  double side_two = Norm_2( difference( point_in_the_border, p2 ) );
                  double base = Norm_2( difference( p2, p1 ) );
                  double gamma = side_one + side_two - base;
                  min = test;
                  temp = make_pair( notCovered[ i ], gamma );
               }
            }
            else{
               if ( theta <= 0){
                  theta = 0;
                  point = sum_vector ( scalar_product( theta, p1 ), scalar_product( 1 - theta, p2 ) );
                  double test = Norm_2( difference( point, c ) );
                  if( test < min ){
                     double lambda = objectOfData->getRadius( notCovered[ i ] )/test;
                     point_in_the_border = sum_vector ( scalar_product( lambda, point ), scalar_product( 1 - lambda, c ) );
                     double side_one = Norm_2( difference( point_in_the_border, p1 ) );
                     double side_two = Norm_2( difference( point_in_the_border, p2 ) );
                     double base = Norm_2( difference( p2, p1 ) );
                     double gamma = side_one + side_two - base;
                     min = test;
                     temp = make_pair( notCovered[ i ], gamma );
                  }
               }
            }
         }
      }

      min = DBL_MAX;
      aux.push_back( temp );
   }

   sort( aux.begin(), aux.end(), sortVectorOfPairs );

   for ( int i = 0; i < notCovered.size(); i++ ){
      notCoveredBalls[ i ] = aux[ i ].first;
   }
}

vector< int > BranchNBound::insert ( vector< int > aux, int i, int k )
{
   vector< int >::iterator it;
   it = aux.begin();
   for ( int j = 0; j < i + 1; j++ ){
      it++;
   }

   aux.insert ( it , k );

   return aux;

}

void BranchNBound::computeLowerBounds( list< node* > * open_nodes, node *current_node, double * best_lb_aux )
{
   list < node* >:: iterator itOpen_aux;
   list < node* >:: iterator aux;	
   current_node->s_lb = 0;						
   double min = DBL_MAX;

   for ( itOpen_aux = ( *open_nodes ).begin(); itOpen_aux != ( *open_nodes ).end(); itOpen_aux++ ){
      if( ( *itOpen_aux )->lb < min ){
         min = ( *itOpen_aux )->lb;
         aux = itOpen_aux;
      }
   }

   ( *aux )->s_lb = 1;
   *best_lb_aux = min;
}

void BranchNBound::printLog( mpf_class temp_aux, mpz_t treeSize, node* child_aux, list< node* > open_aux, unsigned long int count_socp, mpz_class quantity_aux, double b_ub, double b_lb, int notCovSize, int * printHeader2 )
{
   temp_aux = mpz_class( treeSize );
   temp_aux = ( quantity_aux/temp_aux )*100;

   cout.unsetf(ios_base::floatfield);

   if( *printHeader2 == 0 ){
      *printHeader2 = 1;

      cout 	<< "F. value\t"
         << "BestUB\t"
         << "        BestLB\t"
         << "        GAP%(I)\t"
         << "        Level\t"
         << "#Uncov\t"
         << "#Solved\t"
         << " open\t"
         << "        Tree (%)\t" << endl;
   }
   cout 	<< setfill( '0' ) << setw( 3 ) << child_aux->lb << "\t";
   if ( b_ub > pow( 10, 100 ) )	{
      cout << setfill( '\t' ) << setw( 3 ) << "- " << "\t";
   }
   else{
      cout << setfill( '0' ) << setw( 3 ) << b_ub << "\t";
   }
   if ( b_lb > pow( 10, 100 ) )	{
      cout << setfill( '\t' ) << setw( 3 ) << "- " << "\t";
   }
   else{
      cout << setfill( '0' ) << setw( 3 ) << b_lb << "\t";
   }
   double gap = ( ( b_ub - b_lb )/b_lb )*100;
   if ( gap > pow( 10, 100 ) || gap == 0 )	{
      cout << setfill( '\t' ) << setw( 3 ) << "- " << "\t";
   }
   else{
      cout << setfill( '0' ) << setw( 3 ) << ( ( b_ub - b_lb )/b_lb )*100 << "\t";
   }
   cout << setfill( '0' ) << setw( 3 ) << child_aux->pts.size() + 1 - 3 << "\t"
      << setfill( '0' ) << setw( 3 ) << notCovSize << "\t"
      << setfill( '0' ) << setw( 7 ) << count_socp << "\t "
      << setfill( '0' ) << setw( 7 ) << open_aux.size() << "\t"
      << setfill( '0' ) << setw( 3 ) << temp_aux << "%" << endl;
}

void BranchNBound::printLog( mpf_class temp_aux, mpz_t treeSize, node* child_aux, unsigned long int num_unexplrd, unsigned long int count_socp, mpz_class quantity_aux, double b_ub, double b_lb, int notCovSize, int * printHeader2 )
{
   temp_aux = mpz_class( treeSize );
   temp_aux = ( quantity_aux/temp_aux )*100;

   cout.unsetf(ios_base::floatfield);

   if( *printHeader2 == 0 ){
      *printHeader2 = 1;

      cout 	<< "F. value\t"
         << "BestUB\t"
         << "        BestLB\t"
         << "        GAP%(I)\t"
         << "        Level\t"
         << "#Uncov\t"
         << "#Solved\t"
         << " open\t"
         << "        Tree (%)\t" << endl;
   }
   cout 	<< setfill( '0' ) << setw( 3 ) << child_aux->lb << "\t";
   if ( b_ub > pow( 10, 100 ) )	{
      cout << setfill( '\t' ) << setw( 3 ) << "- " << "\t";
   }
   else{
      cout << setfill( '0' ) << setw( 3 ) << b_ub << "\t";
   }
   if ( b_lb > pow( 10, 100 ) )	{
      cout << setfill( '\t' ) << setw( 3 ) << "- " << "\t";
   }
   else{
      cout << setfill( '0' ) << setw( 3 ) << b_lb << "\t";
   }
   double gap = ( ( b_ub - b_lb )/b_lb )*100;
   if ( gap > pow( 10, 100 ) || gap == 0 )	{
      cout << setfill( '\t' ) << setw( 3 ) << "- " << "\t";
   }
   else{
      cout << setfill( '0' ) << setw( 3 ) << ( ( b_ub - b_lb )/b_lb )*100 << "\t";
   }
   cout << setfill( '0' ) << setw( 3 ) << child_aux->pts.size() + 1 - 3 << "\t"
      << setfill( '0' ) << setw( 3 ) << notCovSize << "\t"
      << setfill( '0' ) << setw( 7 ) << count_socp << "\t "
      << setfill( '0' ) << setw( 12 ) << num_unexplrd << "\t"
      << setfill( '0' ) << setw( 3 ) << temp_aux << "%" << endl;
}

void BranchNBound::computeSizeTree( int sizeInst, mpz_t sizeOfTree, vector< mpz_class > &levels )
{
   mpz_init ( sizeOfTree );

   int heightOfTree = sizeInst - 2;

   for ( int i = 0; i < heightOfTree; i++ ){
      mpz_t factorial;
      mpz_init ( factorial );
      unsigned long int n = i + 2;		
      mpz_fac_ui ( factorial, n );
      levels[ i ] = mpz_class( factorial )/2;
      mpz_add ( sizeOfTree, sizeOfTree, factorial );
      mpz_clear ( factorial );
   }
   unsigned long int div = 2;
   mpz_div_ui ( sizeOfTree, sizeOfTree , div );

   mpf_class print;
   print = mpz_class( sizeOfTree );

   cout << scientific;
   cout << "Size of the tree: " << print << endl;
}

void BranchNBound::makeBranching( branching *branchingS )
{

}

bool BranchNBound::crossRoads( vector< int > & sequence, vector< vector< double > > & solution )
{	
   return false;
}

void BranchNBound::setBranchingRuleList2()
{
   vector< pair< int, double > > dist;
   dist.resize( sizeOfInstance );

   for ( int i = 0; i < sizeOfInstance; i++ ){
      ( dist[ i ] ).first = i;
      ( dist[ i ] ).second = sqrt( pow( objectOfData->getCoordx( 0 ) - objectOfData->getCoordx( i ) ,2 ) + pow( objectOfData->getCoordy( 0 ) - objectOfData->getCoordy( i ),2 ) );

   }

   sort( dist.begin(), dist.end(), sortVectorOfPairs );

   for ( int i = 0; i < dist.size(); i++ ){
      branchingRule2.push_back( ( dist[ i ] ).first );
   }
}

int BranchNBound::strongBranching( list< int > nCovered, vector< int > sbSequence )
{
   vector< pair < int, double > > sb_aux;
   pair < int, double > sbPair;
   list< int >::iterator sbIt;
   sbIt = nCovered.begin();
   int sbRange = 4;
   double sum = 0;

   for ( int i = 0; i < sbRange; i++ ){
      int k = ( *sbIt );
      sbIt++;
      for ( int j = 0; j < sbSequence.size(); j++ ){
         node *sbChild = new node;
         sbChild->pts.push_back( 0 );
         sbChild->pts = insert( sbSequence, j, k );
         SolveSocpCplex *sbSOCP = new SolveSocpCplex( objectOfData, sbChild->pts );	
         sbSOCP->solveSOCP( sbChild->pts );
         sbChild->lb = sbSOCP->getF_value();
         sum += sbChild->lb;
         delete sbChild;
      }
      sbPair = make_pair( k, sum );
      sb_aux.push_back( sbPair );
      sum = 0;
   }

   sort( sb_aux.begin(), sb_aux.end(), sortVectorOfPairs );
   cout << "teste: ";
   for ( int i = 0; i < sb_aux.size(); i++ ){
      cout << sb_aux[ i ].first << " ";
   }
   cout << endl;
   return sb_aux[ 0 ].first;
}

void BranchNBound::set_branch_rule(int r)
{
   variable_select_rule = r;
}

int BranchNBound::check_reverse_prune(node* n, vector<double>& tempX, vector<double>& tempY, vector<double>& tempZ, int t1, int t2)
{
   int t1_index = n->pts[t1];
   int t2_index = n->pts[t2];
   if (objectOfData->is_intersect(t1_index, t2_index))
      return 0;
   // two turn points that determines the line 
   vector<double> tp1;
   vector<double> tp2;
   vector<double> ct1;
   vector<double> ct2;
   tp1 = {tempX[t1], tempY[t1], tempZ[t1]};
   tp2 = {tempX[t2], tempY[t2], tempZ[t2]};
   ct1 = {objectOfData->getCoordx(t1_index),
          objectOfData->getCoordy(t1_index),
          objectOfData->getCoordz(t1_index)};
   ct2 = {objectOfData->getCoordx(t2_index),
          objectOfData->getCoordy(t2_index),
          objectOfData->getCoordz(t2_index)};
   
   // check intersection to nodes t1_index and t2_index
   double aux_a = pow(tp2[0]-tp1[0],2) + pow(tp2[1]-tp1[1],2) + pow(tp2[2]-tp1[2],2);
   double aux_b1 = 2 * ((tp2[0]-tp1[0])*(tp1[0]-ct1[0])
                     +  (tp2[1]-tp1[1])*(tp1[1]-ct1[1])
                     +  (tp2[2]-tp1[2])*(tp1[2]-ct1[2]));
   double aux_b2 = 2 * ((tp2[0]-tp1[0])*(tp1[0]-ct2[0])
                     +  (tp2[1]-tp1[1])*(tp1[1]-ct2[0])
                     +  (tp2[2]-tp1[2])*(tp1[2]-ct2[0]));
   double aux_c1 = pow(ct1[0],2) 
                 + pow(ct1[1],2)
                 + pow(ct1[2],2)
                 + pow(tp1[0],2) + pow(tp1[1],2) + pow(tp1[2],2)
                 - 2 * (ct1[0]*tp1[0] + ct1[1]*tp1[1] + ct1[2]*tp1[2])
                 - pow(objectOfData->getRadius(t1_index),2);
   double aux_c2 = pow(ct2[0],2) 
                 + pow(ct2[1],2)
                 + pow(ct2[2],2)
                 + pow(tp1[0],2) + pow(tp1[1],2) + pow(tp1[2],2)
                 - 2 * (ct2[0]*tp1[0] + ct2[1]*tp1[1] + ct2[2]*tp1[2])
                 - pow(objectOfData->getRadius(t2_index),2);
   double aux_diff1 = pow(aux_b1,2) - 4*aux_a*aux_c1;
   double aux_diff2 = pow(aux_b2,2) - 4*aux_a*aux_c2;

   double aux_u1_in = (-aux_b1-sqrt(aux_diff1)) / (2*aux_a);
   double aux_u1_out = (-aux_b1+sqrt(aux_diff1)) / (2*aux_a);
   double aux_u2_in = (-aux_b2-sqrt(aux_diff2)) / (2*aux_a);
   double aux_u2_out = (-aux_b2+sqrt(aux_diff2)) / (2*aux_a);

   // TODO: perhaps we should introduce float number comparison tolerance.
   if (aux_u1_in <= aux_u2_in)
   {
      if (aux_u1_out < aux_u2_in)
         return 0;
      else
         return 1;
   }
   else if (aux_u2_in < aux_u1_in)
   {
      if (aux_u2_out < aux_u1_in)
         return 0;
      else 
         return 2;
   }
}

void BranchNBound::compute_intersects(int node_id, 
   double& intersect_in, double& intersect_out, 
   vector<double>& line_from, vector<double>& line_to)
{
   vector<double> node_center = {objectOfData->getCoordx(node_id), objectOfData->getCoordy(node_id), objectOfData->getCoordz(node_id)};

   double aux_a = pow(line_to[0]-line_from[0],2) + pow(line_to[1]-line_from[1],2) + pow(line_to[2]-line_from[2],2);
   double aux_b1 = 2 * ((line_to[0]-line_from[0])*(line_from[0]-node_center[0])
                     +  (line_to[1]-line_from[1])*(line_from[1]-node_center[1])
                     +  (line_to[2]-line_from[2])*(line_from[2]-node_center[2]));
   double aux_c1 = pow(node_center[0],2) + pow(node_center[1],2) + pow(node_center[2],2)
                 + pow(line_from[0],2) + pow(line_from[1],2) + pow(line_from[2],2)
                 - 2 * (node_center[0]*line_from[0] + node_center[1]*line_from[1] + node_center[2]*line_from[2])
                 - pow(objectOfData->getRadius(node_id),2);
   double aux_diff1 = pow(aux_b1,2) - 4*aux_a*aux_c1;

   intersect_in = (-aux_b1-sqrt(aux_diff1)) / (2*aux_a);
   intersect_out = (-aux_b1+sqrt(aux_diff1)) / (2*aux_a);
}

void BranchNBound::construct_full_sequence(
   vector<double>& solX, vector<double>& solY, vector<double>& solZ, 
   unordered_map<int, vector<int>>& covered_node_ids_by_location, vector<int>& org_sequence, 
   vector<int>& full_sequence)
{
   double intersect_in_value, intersect_out_value;
   vector<pair<int, double>> covered_node_dists;
   full_sequence.clear();
   for (int i = 0; i < org_sequence.size(); i++)
   {
      full_sequence.push_back(org_sequence[i]);
      if (covered_node_ids_by_location.count(i) > 0 && !covered_node_ids_by_location[i].empty())
      {
         covered_node_dists.clear();
         // solX, solY, solZ has one more item than in org_sequence as the last item is 
         // the same as the first to represent deopt location
         vector<double> line_from = {solX[i], solY[i], solZ[i]};
         vector<double> line_to = {solX[i+1], solY[i+1], solZ[i+1]};
         for (int node_id : covered_node_ids_by_location[i])
         {
            compute_intersects(node_id, intersect_in_value, intersect_out_value, 
               line_from, line_to);
            covered_node_dists.push_back(make_pair(node_id, intersect_in_value));
         }
         sort(covered_node_dists.begin(), covered_node_dists.end(), sortVectorOfPairsMin);
         for (int j = 0; j < covered_node_dists.size(); j++)
         {
            full_sequence.push_back(covered_node_dists[j].first);
         }
      }
   }
}

void BranchNBound::construct_full_sequence(vector<double>& solX, vector<double>& solY, vector<double>& solZ, 
               unordered_map<int, vector<int>>& covered_node_ids_by_location, 
               vector<int>& org_sequence, vector<int>& full_sequence,
               vector<vector<double>>& sequence_coords)
{
   double intersect_in_value, intersect_out_value;
   vector<pair<int, double>> covered_node_dists;
   full_sequence.clear();
   sequence_coords.clear();
   sequence_coords.resize(3);
   for (int i = 0; i < org_sequence.size(); i++)
   {
      full_sequence.push_back(org_sequence[i]);
      // sequence_coords.push_back({solX[i], solY[i], solZ[i]});
      sequence_coords[0].push_back(solX[i]); sequence_coords[1].push_back(solY[i]); sequence_coords[2].push_back(solZ[i]);
      if (covered_node_ids_by_location.count(i) > 0 && !covered_node_ids_by_location[i].empty())
      {
         covered_node_dists.clear();
         // solX, solY, solZ has one more item than in org_sequence as the last item is 
         // the same as the first to represent deopt location
         vector<double> line_from = {solX[i], solY[i], solZ[i]};
         vector<double> line_to = {solX[i+1], solY[i+1], solZ[i+1]};
         vector<double> coeff = {line_to[0]-line_from[0], line_to[1]-line_from[1], line_to[2]-line_from[2]};
         for (int node_id : covered_node_ids_by_location[i])
         {
            compute_intersects(node_id, intersect_in_value, intersect_out_value, 
               line_from, line_to);
            // intersect_in_value = min(1.0, max(0.0, intersect_in_value));
            if (intersect_in_value < 1 && intersect_in_value > 0)
            {
               covered_node_dists.push_back(make_pair(node_id, intersect_in_value));
            }
         }
         sort(covered_node_dists.begin(), covered_node_dists.end(), sortVectorOfPairsMin);
         for (int j = 0; j < covered_node_dists.size(); j++)
         {
            double temp_coeff_multi = covered_node_dists[j].second;
            full_sequence.push_back(covered_node_dists[j].first);
            sequence_coords[0].push_back(line_from[0] + temp_coeff_multi * coeff[0]);
            sequence_coords[1].push_back(line_from[1] + temp_coeff_multi * coeff[1]);
            sequence_coords[2].push_back(line_from[2] + temp_coeff_multi * coeff[2]);
         }
      }
   }
}

void BranchNBound::update_insert_spans(vector<pair<int, int>>& new_spans, int t1)
{
   for (int i = 0; i < new_spans.size(); i++)
   {
      if (t1 > new_spans[i].second)
      {
         continue;
      }
      else if (t1 >= new_spans[i].first && t1 < new_spans[i].second)
      {
         new_spans[i].second++;
      }
      else
      {
         new_spans[i].first++;
         new_spans[i].second++;
      }
   }
}

int BranchNBound::get_branch_node_id(vector<int>& uncovered, vector<int>& sequence)
{
   int branch_node_id = -1;
   for(auto temp_node_ptr = uncovered.begin(); temp_node_ptr != uncovered.end(); temp_node_ptr++)
   {
      int temp_id = (*temp_node_ptr);
      bool intersect_flag = false;
      for (int inseq_id : sequence)
      {
         if (objectOfData->is_intersect(temp_id, inseq_id))
         {
            intersect_flag = true;
            break;
         }
      }
      if (!intersect_flag)
      {
         branch_node_id = temp_id;
         break;
      }
   }
   if (branch_node_id == -1)
   {
      branch_node_id = (*uncovered.begin());
   }
   return branch_node_id;
}

void BranchNBound::add_lb(double lb)
{
   if (m_LB_Count.count(lb) == 0)
        m_LB_Count[lb] = 1;
    else
        m_LB_Count[lb]++;
}

void BranchNBound::remove_lb(double lb)
{
   if (m_LB_Count.count(lb) > 0)
   {
      m_LB_Count[lb]--;
      if (m_LB_Count[lb] == 0)
      {
         m_LB_Count.erase(lb);
      }
   }
}