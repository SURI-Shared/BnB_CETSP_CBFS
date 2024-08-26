#ifndef SOLVESOCPSCSRECYCLING_H
#define SOLVESOCPSCSRECYCLING_H

#include "SolveSocpSCS.h"
#include <unordered_map>

class SolveSocpSCSWithRecycling: public SolveSocpSCS{
   public:
      //constructor
      SolveSocpSCSWithRecycling( Data *, vector< int >&);
      SolveSocpSCSWithRecycling( Data *,int);
      ~SolveSocpSCSWithRecycling();
      //functions
      void clear_removable_constraints();
      void clear_removable_constraints(int prev_pos, int curr_pos);
      void solve_warm(double*, double*, double*);
      void accumulate_info();
      size_t solvers_made(){return solvers.size();}
      using SolveSocpSCS::solveSOCP;
      virtual void solveSOCP() override;

      SolveSocpSCSStatistics info_struct;
   protected:
      unordered_map<size_t,ScsWork*> solvers;
      void createModel( vector< int >& ) override;
      void update_b();
};
#endif