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

class SolveSocpClarabel{

   public:

      //construtor	
      SolveSocpClarabel( Data *, vector< int >& ); //constructor
      SolveSocpClarabel( Data *,int);
      //functions
      void solveSOCP( vector< int >& );
      void solveSOCP();
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

      Eigen::Map<Eigen::VectorXd> primals(){solution.x;}
      Eigen::Map<Eigen::VectorXd> duals(){solution.z;}
      Eigen::Map<Eigen::VectorXd> slacks(){solution.s;}

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
