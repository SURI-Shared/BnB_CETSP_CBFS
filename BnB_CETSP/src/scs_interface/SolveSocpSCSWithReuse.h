#ifndef SOLVESOCPSCSRECYCLING_H
#define SOLVESOCPSCSRECYCLING_H

#include "SolveSocpSCS.h"
#include <unordered_map>

class SolveSocpSCSWithReuse: public SolveSocpSCS{
   public:
      //constructor
      SolveSocpSCSWithReuse( Data *, vector< int >&);
      SolveSocpSCSWithReuse( Data *,int);
      SolveSocpSCSWithReuse( Data *, vector< int >&,bool);
      SolveSocpSCSWithReuse( Data *,int,bool);
      ~SolveSocpSCSWithReuse();
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