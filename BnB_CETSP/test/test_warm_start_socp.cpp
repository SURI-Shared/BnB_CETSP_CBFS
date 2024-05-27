#include "warm_start_socp.h"
#include "clarabel_interface/SolveSocpClarabelWithRecycling.h"

void print_vector(const double* vector, size_t length){
    for(size_t i=0;i<length-1;i++){
        cout<<vector[i]<<" ";
    }
    cout<<vector[length-1];
}

void print_vector(const Eigen::Map<Eigen::VectorXd>& vector, size_t length){
    for(size_t i=0;i<length-1;i++){
        cout<<vector[i]<<" ";
    }
    cout<<vector[length-1];
}

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

    solver_handle.clear_removable_constraints(0, 0);
    solver_handle.populate_removable_constraints(child_sequence, 0, 0);

    solver_handle.solve_warm(primal_guess.data(),slack_guess.data(),dual_guess.data(),0,1);

    // Eigen::Array2Xd primal_comp(primal_guess.size());
    // for(size_t i=0;i<primal_guess.size();i++){
    //     primal_comp.row(0)[i]=primal_guess.at(i);
    // }
    // primal_comp.row(1)=solver_handle.primals();
    // cout<<"Primal Guess vs Primal Sol"<<endl;
    // cout<<primal_comp<<endl;
    cout<<"Primal Guess: ";
    print_vector(primal_guess.data(),primal_guess.size());
    cout<<endl;
    cout<<"Primal Solut: ";
    print_vector(solver_handle.primals(),primal_guess.size());
    cout<<endl;

    cout<<"Slack Guess: ";
    print_vector(slack_guess.data(),slack_guess.size());
    cout<<endl;
    cout<<"Slack Solut: ";
    print_vector(solver_handle.slacks(),slack_guess.size());
    cout<<endl;

    cout<<"Parent Dual: ";
    print_vector(parent.duals.data(),parent.duals.size());
    cout<<endl;
    cout<<"Dual Guess: ";
    print_vector(dual_guess.data(),dual_guess.size());
    cout<<endl;
    cout<<"Dual Solut: ";
    print_vector(solver_handle.duals(),dual_guess.size());
    cout<<endl;

    solver_handle.printSolution(child_sequence);
    return 0;
}