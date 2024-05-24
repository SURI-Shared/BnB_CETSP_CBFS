#include "warm_start_socp.h"
#include <algorithm>

WarmStartHandler::WarmStartHandler(){
    //global matrices for the tiny insertion problem
    Eigen::SparseMatrix<double> _insertion_problem_P(5,5);
    Eigen::SparseMatrix<double> _insertion_problem_A(12,5);
    Eigen::Matrix<double,5,1> _insertion_problem_q{1,1,0,0,0};
    Eigen::Array<double,12,1> _insertion_problem_initial_b={0,0,0,0,0,0,0,0,0,0,0,0};//if left uninitialze valgrind complains about a conditional jump or move depending on unintiailized values
    const std::vector<clarabel::SupportedConeT<double>> _insertion_problem_cones={clarabel::SecondOrderConeT<double>(4),clarabel::SecondOrderConeT<double>(4),clarabel::SecondOrderConeT<double>(4)};
    clarabel::DefaultSettings<double> _insertion_problem_settings=clarabel::DefaultSettings<double>::rerunnable_settings();

    //populate constraint matrix using a triplet list
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(11);//2 for the fi, 3 for xi three times

    //||x-c||<=r
    //implemented as s in SOC, s=[r,c1-x1,c2-x2]
    //first a row of zeros in A and r in b
    //then x+s==c
    tripletList.emplace_back(1,2,1);
    tripletList.emplace_back(2,3,1);
    tripletList.emplace_back(3,4,1);

    //||x-p1||<=f1
    //first -f1+s=0
    tripletList.emplace_back(4,0,-1);
    //then -x+s=-p1
    tripletList.emplace_back(5,2,-1);
    tripletList.emplace_back(6,3,-1);
    tripletList.emplace_back(7,4,-1);

    //||x-p2||<=f2
    //first -f2+s=0
    tripletList.emplace_back(8,1,-1);
    //then -x+s=-p2
    tripletList.emplace_back(9,2,-1);
    tripletList.emplace_back(10,3,-1);
    tripletList.emplace_back(11,4,-1);

    _insertion_problem_A.setFromTriplets(tripletList.cbegin(),tripletList.cend());

    insertion_problem_solver_ptr=new clarabel::DefaultSolver<double>(_insertion_problem_P,Eigen::Ref<Eigen::VectorXd>(_insertion_problem_q),
                                                              _insertion_problem_A,Eigen::Ref<Eigen::VectorXd>(_insertion_problem_initial_b),
                                                              _insertion_problem_cones,_insertion_problem_settings);
}

WarmStartHandler::~WarmStartHandler(){
    delete insertion_problem_solver_ptr;
}
void WarmStartHandler::construct_initial_guess(const std::vector<int>& current_sequence, const BnBNodeForWarmStart& parent,
                                               Data* instanceData,std::vector<double>& out_primals,std::vector<double>& out_slacks, std::vector<double>& out_duals){
    std::vector<Eigen::Vector3d> parent_turning_points=parent.turning_points();
    this->construct_initial_guess(current_sequence,parent.sequence,parent_turning_points,parent.duals,instanceData,out_primals,out_slacks,out_duals);
}

void WarmStartHandler::construct_initial_guess(const std::vector<int>& current_sequence, 
                             const std::vector<int>& parent_sequence, const std::vector<Eigen::Vector3d> & parent_turning_points, const Eigen::Ref<const Eigen::VectorX<double>>& parent_dual,
                             Data* instanceData,std::vector<double>& out_primals,std::vector<double>& out_slacks, std::vector<double>& out_duals){
    //the actual decision variables are [fi...xi], which are the travel distances and the turning point coordinates respectively
    //we use wi as shorthand for the displacement vector between x(i-1) and xi, wi=xi-x(i-1)
    //and vi as shorthand for the displacement vector from xi to ci (the center of neighborhood i), vi=ci-xi
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

    out_primals.resize(nvar);//TODO: might be able to get away with some kind of reserve shenanigans to avoid default initializing the output in construct_initial_guess
    out_slacks.resize(nslack);
    out_duals.resize(ndual);

    //figure out guesses for turning points (either copied from parent, or via a small SOCP)
    //for the new turning point, also computes the dual variable related to the constraint that the turning point stay in the neighborhood
    bool not_yet_inserted=true;
    bool not_in_parent=false;
    std::vector<int> index_in_parent(current_sequence.size());
    int cnid;
    int pnid;
    int pidx;
    Eigen::Vector3d depot={instanceData->getCoordx(current_sequence[0]),instanceData->getCoordy(current_sequence[0]),instanceData->getCoordz(current_sequence[0])};
    Eigen::Vector3d last_tp=depot;
    Eigen::Vector3d next_tp;
    for(size_t i=0;i<m;i++){//first entry in sequence is just the depot
        if (i<m-1&&not_yet_inserted){
            cnid=current_sequence.at(i+1);
            pnid=parent_sequence.at(i+1);
            pidx=i+1;
            not_in_parent=pnid!=cnid;
            //once the neighborhood in current_sequence differs from that in parent, we have found the spot where a new neighborhood was inserted
            not_yet_inserted=!not_in_parent;
        }else{
            cnid=current_sequence.at(i+1);
            pnid=parent_sequence.at(i);
            not_in_parent=pnid!=cnid;
            pidx=i;
        }
        if(not_in_parent){
            index_in_parent.at(i+1)=-1;
            //guess the turning point to be at the optimal location if all other turning points are held fixed
            if(i+2>=current_sequence.size()){
                next_tp=depot;
            }else{
                next_tp=parent_turning_points[i];
            }
            Eigen::Vector3d center={instanceData->getCoordx(cnid),instanceData->getCoordy(cnid),instanceData->getCoordz(cnid)};
            double radius=instanceData->getRadius(cnid);
            solve_insertion_problem(last_tp,next_tp,center,radius,out_primals.data()+xstart+ndim*i,out_duals.data()+i*per_turn);
        }
        else{
            index_in_parent.at(i+1)=pidx;
            //turning points in previously considered neighborhoods are guessed to not change
            last_tp=parent_turning_points.at(pidx-1);
            out_primals[xstart+ndim*i]=last_tp[0];
            out_primals[xstart+ndim*i+1]=last_tp[1];
            out_primals[xstart+ndim*i+2]=last_tp[2];
        }
    }

    //populate f guesses
    last_tp=depot;
    std::vector<Eigen::Vector3d> displacement_guesses;displacement_guesses.reserve(nf);
    for(size_t i=0;i<m;i++){
        next_tp=Eigen::Map<Eigen::Vector3d>(out_primals.data()+xstart+i*ndim);
        displacement_guesses.push_back(next_tp-last_tp);
        out_primals[i]=displacement_guesses[i].norm();
        last_tp=next_tp;
    }
    //last waypoint goes back to the depot
    displacement_guesses.push_back(depot-last_tp);
    out_primals[m]=displacement_guesses.at(m).norm();

    //slack guess
    //set the bounding magnitude of vi slack as s=[radii[i], ci-xguessi^T]^T
    for(size_t i=0;i<m;i++){
        int cnid=current_sequence[i+1];//first entry in sequence is just the depot
        out_slacks[i*per_turn]=instanceData->getRadius(cnid);

        out_slacks[i*per_turn+1]=instanceData->getCoordx(cnid)-out_primals[xstart+ndim*i];
        out_slacks[i*per_turn+2]=instanceData->getCoordy(cnid)-out_primals[xstart+ndim*i+1];
        out_slacks[i*per_turn+3]=instanceData->getCoordz(cnid)-out_primals[xstart+ndim*i+2];
    }

    //set the bounding magnitude of wi slack as [fi wguessi^T]
    size_t wmagstart=nx+m;
    for (size_t i=0;i<nf;i++){
        out_slacks[wmagstart+i*per_turn]=out_primals[i];
        out_slacks[wmagstart+i*per_turn+1]=displacement_guesses[i][0];
        out_slacks[wmagstart+i*per_turn+2]=displacement_guesses[i][1];
        out_slacks[wmagstart+i*per_turn+3]=displacement_guesses[i][2];
    }

    //dual guess
    //dual variables are paired with slack variables
    //for those constraints that persist from the parent we can just reuse the parent's dual solution
    //constraints are (in order that they appear in the constraint matrix A)
    //m SOC constraints bounding magnitude of ci-xi (ndim+1 rows per SOC constraint, so nx+m total rows)
    //nf SOC constraints relating magnitude of wi to fi (ndim+1 rows per SOC constraint, so nw+nf total rows)

    //ci-xi
    for(size_t i=0;i<m;i++){
        if (index_in_parent.at(i+1)>=0){//first entry in sequence is just the depot
            //reuse dual values for SOC constraints bounding the magnitude of ci-xi
            for(size_t j=0;j<per_turn;j++){
                out_duals[i*per_turn+j]=parent_dual[(index_in_parent.at(i+1)-1)*per_turn+j];
            }
        }else{
            //SOC constraints enforcing ||ci-xi||<=ri
            //at the solution the slacks will be [ri ci-xi] and the duals [lambda_i, -lambda_i*(ci-xi)/||ci-xi||]
            //where lamgda_i is the lagrange multiplier of the constraint ||ci-xu||<=ri
            //consider the optimization to pick the ith turning point holding the turn before and after it constant
            //this optimization has cost function (f_i+f_{i+1})/2 and its only constraint is ||ci-xi||**2<=ri**2
            //from the KKT conditions, specifically the 0 gradient of the Lagrangian, we obtain:
            //w_i/f_i/2-w_{i+1}/f_{i+1}/2+2*lambda_i*(x-c)=0
            Eigen::Map<Eigen::Vector3d> vguess(out_slacks.data()+wmagstart+i*per_turn+1);
            double lambda_i=out_duals[i*per_turn];//populated by solve_insertion_problem
            double d2center=vguess.norm();
            out_duals[i*per_turn+1]=-lambda_i*vguess[0]/d2center;
            out_duals[i*per_turn+2]=-lambda_i*vguess[1]/d2center;
            out_duals[i*per_turn+3]=-lambda_i*vguess[2]/d2center;
        }
    }

    //||wi||<=fi
    for(size_t i=0;i<nf;i++){
        //At the solution the slacks will be [fi wi.T] and the duals are [1 -wi.T/fi] so we use the wi and fi guess to construct the dual guess
        out_duals[wmagstart+i*per_turn]=1;
        
        out_duals[wmagstart+i*per_turn+1]=-displacement_guesses[i][0]/out_primals[i];
        out_duals[wmagstart+i*per_turn+2]=-displacement_guesses[i][1]/out_primals[i];
        out_duals[wmagstart+i*per_turn+3]=-displacement_guesses[i][2]/out_primals[i];
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
    insertion_problem_solver_ptr->update_b(new_b);

    insertion_problem_solver_ptr->solve();
    clarabel::DefaultSolution<double> solution=insertion_problem_solver_ptr->solution();

    for (size_t i=0;i<3;i++){
        out_turning_point[i]=solution.x[2+i];
    }
    *out_dual=solution.z[0];

}

BnBNodeForWarmStart::BnBNodeForWarmStart(const std::vector<int>& sequence_in,
                                        const Eigen::Map<Eigen::VectorX<double>>& primals,const Eigen::Map<Eigen::VectorX<double>>& duals):
                                        sequence(sequence_in),primals(primals),duals(duals){}

std::vector<Eigen::Vector3d> BnBNodeForWarmStart::turning_points() const{
    size_t ndim=3;
    size_t per_turn=ndim+1;
    size_t m=sequence.size()-1;
    size_t nf=m+1;

    size_t xstart=nf;

    std::vector<Eigen::Vector3d> tps;tps.reserve(m);
    for(size_t i=0;i<m;i++){
        tps.emplace_back(primals[xstart+i*ndim],primals[xstart+i*ndim+1],primals[xstart+i*ndim+2]);
    }
    return tps;
}