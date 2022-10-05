#ifndef SOLVESOCPCPLEX_H
#define SOLVESOCPCPLEX_H
#define IL_STD 1

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
#include<ilcplex/ilocplex.h>
#include <pthread.h>

#include"Data.h"

using namespace std;

class SolveSocpCplex{

   public:

      //construtor	
      SolveSocpCplex( Data *, vector< int >& ); //constructor
      SolveSocpCplex(Data*, int size_seq);
      ~SolveSocpCplex();

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
      IloNum violation;
      int m_num_solves = 0;

   private:

      IloEnv env;
      IloModel model; 
      IloCplex SOCP;
      IloNumArray xCoord;
      IloNumArray yCoord;
      IloNumArray zCoord;
      IloNumVarArray x; 
      IloNumVarArray y; 
      IloNumVarArray z;
      IloNumVarArray m_allf;
      IloNumVarArray m_allw;
      IloNumVarArray m_allu;
      IloNumVarArray m_allv;
      IloNumVarArray m_alls;
      IloNumVarArray m_allt;
      IloNumVarArray m_allq;

      IloNumVarArray m_allvars;
      IloNumArray m_allvarvals;
      IloNumArray m_allvar_redcosts;

      IloExtractableArray m_dists_to_c;
      IloExtractableArray m_coord_x;
      IloExtractableArray m_coord_y;
      IloExtractableArray m_coord_z;

      Data *objectData;

      int sizeProblem;
      double f_value;

      //functions for solving SOCP
      int setSizeProblem( vector< int >& );		
      void createModel( vector< int >& );
      void setF_value();

};

#endif
