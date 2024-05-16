#include "clarabel_interface/SolveSocpClarabel.h"

int main(int argc, char** argv) 
{
    //Data constructor demands command line arguments
    //[path to instance] [OPTIONS] [OVERLAP FACTOR] [TIME LIMIT] [BRANCHING RULE] [BRANCHING STRATEGY] [ROOT SELECTION] [S.B SIZE]
    vector<char*> load_args={"test_SolveSocpClarabel","/home/ggutow/eclipse-workspace/BnB_CETSP/Behdani/CETSP-50-1.txt","3D","1.0","1000","SB","BeFS","1","1","1"};
    Data *dataptr = new Data("/home/ggutow/eclipse-workspace/BnB_CETSP/Behdani/CETSP-50-1.txt", "3D", 1, 9, load_args.data());

    vector<int> in_order;
    int size=dataptr->getSizeInst();
    for(int i=0;i<size;i++){
        in_order.push_back(i);
    }
    SolveSocpClarabel solver_handle=SolveSocpClarabel(dataptr,in_order);

    solver_handle.solveSOCP();

    solver_handle.printSolution(in_order);
    return 0;
}