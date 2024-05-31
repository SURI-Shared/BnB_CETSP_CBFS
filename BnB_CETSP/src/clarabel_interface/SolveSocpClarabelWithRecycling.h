#ifndef SOLVESOCPCLARABELRECYCLING_H
#define SOLVESOCPCLARABELRECYCLING_H

#include "SolveSocpClarabel.h"
#include <unordered_map>

class SolveSocpClarabelWithRecycling: public SolveSocpClarabel{
   public:
      //constructor
      SolveSocpClarabelWithRecycling( Data *, vector< int >&, bool);
      SolveSocpClarabelWithRecycling( Data *,int, bool);
      ~SolveSocpClarabelWithRecycling();
      //functions
      void clear_removable_constraints();
      void clear_removable_constraints(int prev_pos, int curr_pos);
      void solve_warm(double*, double*, double*, int, double);
      void accumulate_info();
      size_t solvers_made(){return solvers.size();}
      using SolveSocpClarabel::solveSOCP;
      virtual void solveSOCP() override;

      SolveSocpClarabelStatistics info_struct;
   protected:
      unordered_map<size_t,clarabel::DefaultSolver<double>*> solvers;
      void createModel( vector< int >& ) override;
      void update_b();
};
#endif