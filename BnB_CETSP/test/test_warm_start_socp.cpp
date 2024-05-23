#include "warm_start_socp.h"
#include "clarabel_interface/SolveSocpClarabelWithRecycling.h"

int main(int argc, char** argv) 
{
    //Data constructor demands command line arguments
    //[path to instance] [OPTIONS] [OVERLAP FACTOR] [TIME LIMIT] [BRANCHING RULE] [BRANCHING STRATEGY] [ROOT SELECTION] [S.B SIZE]
    char exe[]="test_warm_start_socp";
    char instanceName[]="2D/rotatingDiamonds1.txt";
    char option[]="2D";
    char overlapRato[]="1.0";
    char timeLimit[]="1000";
    char branchingRule[]="SB";
    char search_strategy[]="BeFS";
    char root_selection[]="1";
    char sb_size[]="1";
    char* load_args[9]={exe,instanceName,option,overlapRato,timeLimit,branchingRule,search_strategy,root_selection,sb_size};
    Data *dataptr = new Data(instanceName, option, 1, 9, load_args);
    dataptr->printSizes();
    vector<int> in_order;
    int size=3;
    cout<<"Parent sequence: ";
    for(int i=0;i<size;i++){
        in_order.push_back(i);
        cout<<i<<" ";
    }
    cout<<endl;
    SolveSocpClarabelWithRecycling solver_handle=SolveSocpClarabelWithRecycling(dataptr,in_order);
    solver_handle.solveSOCP();

    BnBNodeForWarmStart parent(in_order,solver_handle.primals(),solver_handle.duals());

    vector<int> child_sequence;
    for(int i=0;i<size;i++){
        child_sequence.push_back(i);
    }
    size_t insertion_index=2;
    child_sequence.insert(child_sequence.begin()+insertion_index,size);
    cout<<"Child sequence: ";
    for(auto item:child_sequence){
        cout<<item<<" ";
    }
    cout<<endl;

    WarmStartHandler warm_start_handler;
    std::vector<double> primal_guess;
    std::vector<double> slack_guess;
    std::vector<double> dual_guess;
    warm_start_handler.construct_initial_guess(child_sequence,parent,dataptr,primal_guess,slack_guess,dual_guess);

    cout<<"Primal Guess length: "<<primal_guess.size()<<endl;
    cout<<"Slack Guess length: "<<slack_guess.size()<<endl;
    cout<<"Dual Guess length: "<<dual_guess.size()<<endl;

    solver_handle.clear_removable_constraints(0, 0);
    solver_handle.populate_removable_constraints(child_sequence, 0, 0);

    solver_handle.solve_warm(primal_guess.data(),slack_guess.data(),dual_guess.data());

    solver_handle.printSolution(child_sequence);
    return 0;
}