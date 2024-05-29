#ifndef SOLVEREDUNDANTSOCPCLARABEL_H
#define SOLVEREDUNDANTSOCPCLARABEL_H

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
#include "clarabel_interface/SolveSocpClarabel.h"

using namespace std;

class SolveRedundantSocpClarabel: public SolveSocpClarabel{
    public:

      //constructor
      SolveRedundantSocpClarabel( Data *,int);
      SolveRedundantSocpClarabel( Data *,int, bool);
      //constructors passed a sequence are not provided since we cannot call the overriden createModel during super class construction

   protected:
      //functions for solving SOCP
      virtual void createModel( vector< int >& ) override;
};

#endif