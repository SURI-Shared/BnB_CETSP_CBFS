#include "warm_start_socp.h"
#include "clarabel_interface/SolveSocpClarabelWithReuse.h"

int main(int argc, char** argv) 
{
    WarmStartHandler warm_start_handler;
    double dual_value;
    double turning_point[3];
    Eigen::Vector3d center={0,0,0};
    Eigen::Vector3d p1={0,2,0};
    Eigen::Vector3d p2={2,0,0};
    double r=1;
    warm_start_handler.solve_insertion_problem(p1,p2,center,r,turning_point,&dual_value);
    cout<<"Turning Point: "<<turning_point[0]<<", "<<turning_point[1]<<", "<<turning_point[2]<<endl;
    cout<<"Dual: "<<dual_value<<endl;
    return 0;
}