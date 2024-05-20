#include "warm_start_socp.h"

void WarmStartHandler::construct_initial_guess(const std::vector<int>& current_sequence, 
                             const std::vector<int>& parent_sequence, const Eigen::ArrayX3d & parent_turning_points, Eigen::Map<Eigen::VectorX<double>> parent_dual,
                             Data* instanceData,std::vector<double>& out_primals,std::vector<double>& out_slacks, std::vector<double>& out_duals){
    size_t ndim=3;
    size_t per_turn=ndim+1;
    size_t m=current_sequence.size()-1;
    size_t nf=m+1;
    size_t nx=ndim*m;
    size_t nvar=nf+nx;
    size_t nw=ndim*nf;
    size_t nslack=nx+m+nw+nf;
    size_t ndual=nslack;
    size_t xstart=nf;

    out_primals.reserve(nvar);
    out_slacks.reserve(nslack);
    out_duals.reserve(ndual);

    //figure out guesses for turning points (either copied from parent, or via a small SOCP)
    //for the new turning point, also computes the dual variable related to the constraint that the turning point stay in the neighborhood
    bool not_yet_inserted=true;
    bool not_in_parent=false;
    int cnid;
    int pnid;
    Eigen::Vector3d depot={instanceData->getCoordx(current_sequence[0]),instanceData->getCoordy(current_sequence[0]),instanceData->getCoordz(current_sequence[0])};
    Eigen::Vector3d last_tp=depot;
    Eigen::Vector3d next_tp;
    for(size_t i=0;i<current_sequence.size();i++){
        if (not_yet_inserted){
            cnid=current_sequence.at(i);
            pnid=parent_sequence.at(i);
            not_in_parent=pnid!=cnid;
            //once the neighborhood in current_sequence differs from that in parent, we have found the spot where a new neighborhood was inserted
            not_yet_inserted=not_in_parent;
        }else{
            cnid=current_sequence.at(i);
            pnid=parent_sequence.at(i-1);
            not_in_parent=pnid!=cnid;
        }
        if(not_in_parent){
            //guess the turning point to be at the optimal location if all other turning points are held fixed
            if(i+1>=current_sequence.size()){
                next_tp=depot;
            }else{
                next_tp=parent_turning_points.row(i);
            }
            Eigen::Vector3d center={instanceData->getCoordx(current_sequence[i]),instanceData->getCoordy(current_sequence[i]),instanceData->getCoordz(current_sequence[i])};
            double radius=instanceData->getRadius(current_sequence[i]);
            solve_insertion_problem(last_tp,next_tp,center,radius,out_primals.data()+xstart+ndim*i,out_duals.data()+i*per_turn);
        }
        else{
            //turning points in previously considered neighborhoods are guessed to not change
            last_tp=parent_turning_points.row(i);
            out_primals.at(xstart+ndim*i)=last_tp[0];
            out_primals.at(xstart+ndim*i+1)=last_tp[0];
            out_primals.at(xstart+ndim*i+2)=last_tp[0];
        }
    }

    //populate f guesses
    last_tp=depot;
    for(size_t i=0;i<nf;i++){
        next_tp=out_primals
        Eigen::Vector3d displacement_guess=
    }

}

/*
compute the point in a sphere that is closest to point1 and point2

Parameters: point1 : const Eigen::Vector3d& 
                first point
            point2 : const Eigen::Vector3d&
                second point
            center : const Eigen::Vector3d&
                center of the sphere
            radius : double
                radius of the sphere
            out_turning_point : double* MUTATED
                values are set to the coordinates of the point in the sphere
            out_dual : double* MUTATED
                value is set to the lagrange multiplier of the inequality constraint that says the point stays in the sphere

Mutated:    this->insertion_problem_solver (calls this->insertion_problem_solver.update_b(), this->insertion_problem_solver.solve(), this->insertion_problem_solver.solution())
*/
void WarmStartHandler::solve_insertion_problem(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2, const Eigen::Vector3d& center,double radius,double* out_turning_point,double* out_dual){
    Eigen::Vector<double,12> new_b;
    new_b[0]=radius;
    new_b.segment<3>(1)=center;
    new_b[4]=0;
    new_b.segment<3>(5)=-point1;
    new_b[8]=0;
    new_b.segment<3>(9)=-point2;
    insertion_problem_solver.update_b(new_b);

    insertion_problem_solver.solve();
    clarabel::DefaultSolution<double> solution=insertion_problem_solver.solution();

    for (size_t i=0;i<3;i++){
        out_turning_point[i]=solution.x[2+i];
    }
    *out_dual=solution.z[0];

}