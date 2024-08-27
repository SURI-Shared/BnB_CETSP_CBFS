#ifndef SOLVESOCPSCS_H
#define SOLVESOCPSCS_H

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
#include <scs.h>
#include"Data.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>

using namespace std;

#define SOCP_NDIM 3//TODO: template the SOCP solvers on the dimension of the space

ScsMatrix* no_copy_eigen_SparseMatrix_to_ScsMatrix(Eigen::SparseMatrix<double>&);

struct SolveSocpSCSStatistics{
   public:
      size_t m_num_solves=0;
      size_t solvers_made=0;
      double solve_time=0;
      double setup_time=0;
      double kktinit_time=0;
      double kkt_solve_time=0;
      uint32_t iterations=0;
};
std::ostream& operator<<(std::ostream& os, const SolveSocpSCSStatistics& info_struct);
size_t compute_nvar(size_t);
size_t compute_ncon(size_t);
class SolveSocpSCS{

   public:

      //construtor	
      SolveSocpSCS( Data *, vector< int >& ); //constructor
      SolveSocpSCS( Data *, vector< int >& ,bool); //constructor
      SolveSocpSCS( Data *,int);
      SolveSocpSCS( Data *,int,bool);
      ~SolveSocpSCS();
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

      void accumulate_info(SolveSocpSCSStatistics&);

      Eigen::Map<Eigen::VectorXd> primals(){return Eigen::Map<Eigen::VectorXd>(solution.x,nvar);}
      Eigen::Map<Eigen::VectorXd> duals(){return Eigen::Map<Eigen::VectorXd>(solution.y,ncon);}
      Eigen::Map<Eigen::VectorXd> slacks(){return Eigen::Map<Eigen::VectorXd>(solution.s,ncon);}
      int get_latest_iterations(){return info.iter;}
      double get_latest_solve_time(){return info.solve_time;}
      double get_latest_setup_time(){return info.setup_time;}
      double get_latest_kktinit_time(){return info.init_lin_sys_time;}
      double get_latest_kkt_solve_time(){return info.lin_sys_time;}

   protected:
      
      ScsSettings settings;
      ScsWork* workspace_ptr;
      ScsSolution solution;
      ScsInfo info;

      size_t ncon;
      size_t nvar;

      Data *objectData;

      int sizeProblem;
      double f_value;
      vector<int> sequence;

      bool as_tridiagonal;

      //functions for solving SOCP
      virtual void createModel( vector< int >& );
      virtual void createRegularModel( vector< int >&);
      virtual void createTridiagonalModel( vector< int >& );
};

#endif
