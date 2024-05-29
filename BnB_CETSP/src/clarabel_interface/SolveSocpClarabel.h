#ifndef SOLVESOCPCLARABEL_H
#define SOLVESOCPCLARABEL_H

#include<cmath>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<list>
#include<cstdlib>
#include<ctime>
#include<climits>
#include<algorithm>
#include<stdio.h>
#include <pthread.h>
#include <Clarabel>
#include"Data.h"

using namespace std;

#define SOCP_NDIM 3//TODO: template the SOCP solvers on the dimension of the space

struct SolveSocpClarabelStatistics{
   public:
      size_t m_num_solves=0;
      size_t solvers_made=0;
      double solve_time=0;
      double setup_time=0;
      double equilibration_time=0;
      double kktinit_time=0;
      double initialization_time=0;
      double ip_iteration_time=0;
      double kkt_update_time=0;
      double kkt_solve_time=0;
      uint32_t iterations=0;
};
std::ostream& operator<<(std::ostream& os, const SolveSocpClarabelStatistics& info_struct);
class SolveSocpClarabel{

   public:

      //construtor	
      SolveSocpClarabel( Data *, vector< int >& ); //constructor
      SolveSocpClarabel( Data *, vector< int >& , bool); //constructor
      SolveSocpClarabel( Data *,int);
      SolveSocpClarabel( Data *,int, bool);
      //functions
      virtual void solveSOCP( vector< int >& );
      virtual void solveSOCP();
      void finishSOCP();
      double getF_value();
      double getSolutionX( int );
      double getSolutionY( int );
      double getSolutionZ( int );
      void printF_value();
      void printSolution(  vector<int>& );
      void clear_removable_constraints();
      void clear_removable_constraints(int prev_pos, int curr_pos);
      void populate_removable_constraints(vector<int>& sequence);
      void populate_removable_constraints(vector<int>& sequence, int prev_pos, int curr_pos);
      void initialize_model();
      void initialize_init_point();
      void set_init_point_values();
      double violation;
      int m_num_solves = 0;

      void accumulate_info(SolveSocpClarabelStatistics&);

      Eigen::Map<Eigen::VectorXd> primals(){return solution.x;}
      Eigen::Map<Eigen::VectorXd> duals(){return solution.z;}
      Eigen::Map<Eigen::VectorXd> slacks(){return solution.s;}
      int get_latest_iterations(){return solution.iterations;}
      double get_latest_solve_time(){return solution.solve_time;}
      double get_latest_setup_time(){return solution.setup_time;}
      double get_latest_equilibration_time(){return solution.equilibration_time;}
      double get_latest_kktinit_time(){return solution.kktinit_time;}
      double get_latest_initialization_time(){return solution.initialization_time;}
      double get_latest_ip_iteration_time(){return solution.ip_iteration_time;}
      double get_latest_kkt_update_time(){return solution.kkt_update_time;}
      double get_latest_kkt_solve_time(){return solution.kkt_solve_time;}

   protected:
      
      clarabel::DefaultSettings<double> settings;
      clarabel::DefaultSolver<double>* solver_ptr;
      clarabel::DefaultSolution<double> solution;

      Data *objectData;

      int sizeProblem;
      double f_value;
      vector<int> sequence;

      //functions for solving SOCP
      virtual void createModel( vector< int >& );
};

#endif
