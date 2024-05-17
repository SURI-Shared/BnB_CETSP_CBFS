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
   protected:
      unordered_map<size_t,clarabel::DefaultSolver<double>*> solvers;
      void createModel( vector< int >& ) override;
      void update_b();
};
#endif