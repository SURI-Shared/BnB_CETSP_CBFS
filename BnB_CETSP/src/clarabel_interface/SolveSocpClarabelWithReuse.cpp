#include "SolveSocpClarabelWithReuse.h"

//constructors
SolveSocpClarabelWithReuse::SolveSocpClarabelWithReuse(Data * instance, vector<int>& in_sequence, bool reduced_first_correction):
    SolveSocpClarabel(instance,in_sequence,reduced_first_correction),
    info_struct(){
    solvers.insert({in_sequence.size(),solver_ptr});
}

SolveSocpClarabelWithReuse::SolveSocpClarabelWithReuse(Data * instance, int size_seq, bool reduced_first_correction):
    SolveSocpClarabel(instance,size_seq,reduced_first_correction),
    info_struct(){}

SolveSocpClarabelWithReuse::~SolveSocpClarabelWithReuse(){
    for(auto ptr: solvers){
        delete ptr.second;
    }
}
//methods
/*
void SolveSocpClarabel::createModel(vector<int>& sequence)

store a clarabel::DefaultSolver in this->solver and populate its data for the SOCP to visit neighbhoods in the order defined by sequence

Parameters: sequence : vector<int>&
                the order to visit neighborhoods in. The first element must be the depot, of radius 0. The path returns to the depot after visiting sequence[-1]
*/
void SolveSocpClarabelWithReuse::createModel(vector<int>& sequence){
    sizeProblem=sequence.size();
    if (solvers.count(sizeProblem)){
        solver_ptr=solvers.at(sizeProblem);
        this->sequence=sequence;
        update_b();
    }else{
        SolveSocpClarabel::createModel(sequence);
        solvers.insert({sizeProblem,solver_ptr});
        info_struct.solvers_made++;
    }
}
void SolveSocpClarabelWithReuse::solveSOCP(){
    SolveSocpClarabel::solveSOCP();
    accumulate_info();
}
void SolveSocpClarabelWithReuse::accumulate_info(){
    SolveSocpClarabel::accumulate_info(info_struct);
    info_struct.solvers_made=solvers.size();
}
void SolveSocpClarabelWithReuse::update_b(){
    size_t m=sequence.size()-1;
    size_t nf=m+1;
    size_t nx=SOCP_NDIM*m;
    size_t nvar=nf+nx;
    size_t xstart=nf;
    Eigen::VectorXd b(m*SOCP_NDIM+m+nf*SOCP_NDIM+nf);
    size_t row=0;

    /*||xi-ci||_2<=radii[i] thus becomes s==[radii[i]] in K_SOC
                                            [ci-xi   ]
    expressed as [0]xi+s==[radii[i]]=b
                 [I]      [ci      ]*/
    for(size_t i=0;i<m;i++){
        b[row]=objectData->getRadius( sequence[ i+1 ] );//lhs first entry is the radius
        row+=1;//row of zeros
        //x
        b[row]=objectData->getCoordx(sequence[i+1]);
        row+=1;
        //y
        b[row]=objectData->getCoordy(sequence[i+1]);
        row+=1;
        //z
        b[row]=objectData->getCoordz(sequence[i+1]);
        row+=1;
    }

    /*||x(i+1)j-xij||<=fi becomes s=[fi,x(i+1)-xi] in K_SOC
     expressed as [-1 0  0][fi]+s==0
                  [0  I -I][xi]
                           [x(i+1)]*/
    for(size_t i=0;i<nf;i++){
        b[row]=0;
        row+=1;
        if(i==0){
            //depot to first turn
            //||x1-depot||<=f0
            //s=[f0,x1-depot]
            //x
            b[row]=-objectData->getCoordx(sequence[0]);
            row+=1;
            //y
            b[row]=-objectData->getCoordy(sequence[0]);
            row+=1;
            //z
            b[row]=-objectData->getCoordz(sequence[0]);
            row+=1;

        }else if(i==m){
            //last turn to depot
            //||depot-xm||<=fi
            //s=[fi,depot-xm]
            //x
            b[row]=objectData->getCoordx(sequence[0]);
            row+=1;
            //y
            b[row]=objectData->getCoordy(sequence[0]);
            row+=1;
            //z
            b[row]=objectData->getCoordz(sequence[0]);
            row+=1;
        }else{
            //||x(i+1)j-xij||<=fi
            //s=[fi,x(i+1)-xi]
            for(size_t j=0;j<SOCP_NDIM;j++){
                b[row]=0;
                row+=1;
            }
        }
    }
    solver_ptr->update_b(b);
}
void SolveSocpClarabelWithReuse::clear_removable_constraints(){}
void SolveSocpClarabelWithReuse::clear_removable_constraints(int prev_pos, int curr_pos){}

/*provide a warm start to Clarabel
Parameters: primal_guess : double *
                guess for the primal variables. MUST be the correct length (# of columns in constraint matrix). Copied, not mutated.
            slack_guess : double *
                guess for the slack variables. MUST be the correct length (# of rows in constraint matrix). Copied, not mutated.
            dual_guess : double *
                guess for the dual variables. MUST be the correct length (# of rows in constraint matrix). Copied, not muted.
            mode : int
                how to use the initial guess. See Clarabel implementation for up to date options. If the value is unrecognized, it SHOULD default to just copying the guess and adjusting it slightly to be feasible.
            lambda : double
                for initializations that combine the guess with some other initialization, like Skajaa 2013, lambda is the weight to place on the guess and 1-lambda the weight to place on the other initialization
Mutates:    solution
            f_value
            violation
            m_num_solves
            *solver_ptr
*/
void SolveSocpClarabelWithReuse::solve_warm(double* primal_guess, double* slack_guess, double* dual_guess, int mode, double lambda){
    solver_ptr->solve_warm(primal_guess,slack_guess,dual_guess,mode,lambda);
    new (&solution) clarabel::DefaultSolution<double>(solver_ptr->solution());
    f_value=solution.obj_val;
    violation=solution.r_prim;
    m_num_solves+=1;
    accumulate_info();
}