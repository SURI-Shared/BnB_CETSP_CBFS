#ifndef SOLVESOCPCLARABELRECYCLING_H
#define SOLVESOCPCLARABELRECYCLING_H

#include "SolveSocpClarabel.h"
#include <unordered_map>

class SolveSocpClarabelWithRecycling: public SolveSocpClarabel{
   public:
      //constructor
      SolveSocpClarabelWithRecycling( Data *, vector< int >& );
      SolveSocpClarabelWithRecycling( Data *,int);
      //functions
      void clear_removable_constraints();
      void clear_removable_constraints(int prev_pos, int curr_pos);
      void solve_warm(double*, double*, double*);
      void accumulate_info();
      using SolveSocpClarabel::solveSOCP;
      virtual void solveSOCP() override;

      double solve_time;
      double setup_time;
      double equilibration_time;
      double kktinit_time;
      double initialization_time;
      double ip_iteration_time;
      double kkt_update_time;
      double kkt_solve_time;
      uint32_t iterations;
   protected:
      unordered_map<size_t,clarabel::DefaultSolver<double>*> solvers;
      void createModel( vector< int >& ) override;
      void update_b();
};
#endif