#include "clarabel_interface/SolveSocpClarabel.h"

int main(int argc, char** argv) 
{
    //Data constructor demands command line arguments
    //[path to instance] [OPTIONS] [OVERLAP FACTOR] [TIME LIMIT] [BRANCHING RULE] [BRANCHING STRATEGY] [ROOT SELECTION] [S.B SIZE]
    char exe[]="test_SolveSocpClarabel";
    char instanceName[]="Behdani/CETSP-50-1.txt";
    char option[]="2D";
    char overlapRato[]="1.0";
    char timeLimit[]="1000";
    char branchingRule[]="SB";
    char search_strategy[]="BeFS";
    char root_selection[]="1";
    char sb_size[]="1";
    char* load_args[9]={exe,instanceName,option,overlapRato,timeLimit,branchingRule,search_strategy,root_selection,sb_size};
    Data *dataptr = new Data(instanceName, option, 1, 9, load_args);
    vector<int> in_order;
    int size=3;
    for(int i=0;i<size;i++){
        in_order.push_back(i);
    }
    SolveSocpClarabel solver_handle=SolveSocpClarabel(dataptr,in_order);
    solver_handle.solveSOCP();
    solver_handle.printSolution(in_order);
    return 0;
}