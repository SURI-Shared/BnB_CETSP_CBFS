#include "SolveSocpSCSWithRecycling.h"

//constructors
SolveSocpSCSWithReuse::SolveSocpSCSWithReuse(Data * instance, vector<int>& in_sequence):
    SolveSocpSCS(instance,in_sequence),
    info_struct(){
    //zero out the info struct entries
    info_struct.m_num_solves=0;
    info_struct.solvers_made=0;
    info_struct.solve_time=0;
    info_struct.setup_time=0;
    info_struct.kktinit_time=0;
    info_struct.kkt_solve_time=0;
    info_struct.cone_proj_time=0;
    info_struct.accel_time=0;
    info_struct.iterations=0;
    solvers.insert({in_sequence.size(),workspace_ptr});
    size_t max_nvar=compute_nvar(instance->getSizeInst());
    size_t max_ncon=compute_ncon(instance->getSizeInst());
    solution.x=(double*)calloc(nvar,sizeof(double));
    solution.s=(double*)calloc(ncon,sizeof(double));
    solution.y=(double*)calloc(ncon,sizeof(double));
}

SolveSocpSCSWithReuse::SolveSocpSCSWithReuse(Data * instance, int size_seq):
    SolveSocpSCS(instance,size_seq),
    info_struct(){
    //zero out the info struct entries
    info_struct.m_num_solves=0;
    info_struct.solvers_made=0;
    info_struct.solve_time=0;
    info_struct.setup_time=0;
    info_struct.kktinit_time=0;
    info_struct.kkt_solve_time=0;
    info_struct.cone_proj_time=0;
    info_struct.accel_time=0;
    info_struct.iterations=0;
    size_t max_nvar=compute_nvar(instance->getSizeInst());
    size_t max_ncon=compute_ncon(instance->getSizeInst());
    solution.x=(double*)calloc(max_nvar,sizeof(double));
    solution.s=(double*)calloc(max_ncon,sizeof(double));
    solution.y=(double*)calloc(max_ncon,sizeof(double));
}

SolveSocpSCSWithReuse::SolveSocpSCSWithReuse(Data * instance, vector<int>& in_sequence, bool use_tridiagonal):
    SolveSocpSCS(instance,in_sequence,use_tridiagonal),
    info_struct(){
    //zero out the info struct entries
    info_struct.m_num_solves=0;
    info_struct.solvers_made=0;
    info_struct.solve_time=0;
    info_struct.setup_time=0;
    info_struct.kktinit_time=0;
    info_struct.kkt_solve_time=0;
    info_struct.cone_proj_time=0;
    info_struct.accel_time=0;
    info_struct.iterations=0;
    solvers.insert({in_sequence.size(),workspace_ptr});
    size_t max_nvar=compute_nvar(instance->getSizeInst());
    size_t max_ncon=compute_ncon(instance->getSizeInst());
    solution.x=(double*)calloc(nvar,sizeof(double));
    solution.s=(double*)calloc(ncon,sizeof(double));
    solution.y=(double*)calloc(ncon,sizeof(double));
}

SolveSocpSCSWithReuse::SolveSocpSCSWithReuse(Data * instance, int size_seq, bool use_tridiagonal):
    SolveSocpSCS(instance,size_seq, use_tridiagonal),
    info_struct(){
    //zero out the info struct entries
    info_struct.m_num_solves=0;
    info_struct.solvers_made=0;
    info_struct.solve_time=0;
    info_struct.setup_time=0;
    info_struct.kktinit_time=0;
    info_struct.kkt_solve_time=0;
    info_struct.cone_proj_time=0;
    info_struct.accel_time=0;
    info_struct.iterations=0;
    size_t max_nvar=compute_nvar(instance->getSizeInst());
    size_t max_ncon=compute_ncon(instance->getSizeInst());
    solution.x=(double*)calloc(max_nvar,sizeof(double));
    solution.s=(double*)calloc(max_ncon,sizeof(double));
    solution.y=(double*)calloc(max_ncon,sizeof(double));
}

SolveSocpSCSWithReuse::~SolveSocpSCSWithReuse(){
    for(auto ptr: solvers){
        scs_finish(ptr.second);
    }
    workspace_ptr=nullptr;
}
//methods
/*
void SolveSocpSCS::createModel(vector<int>& sequence)

store a ScsWork in this->workspace_ptr and populate its data for the SOCP to visit neighbhoods in the order defined by sequence

Parameters: sequence : vector<int>&
                the order to visit neighborhoods in. The first element must be the depot, of radius 0. The path returns to the depot after visiting sequence[-1]
*/
void SolveSocpSCSWithReuse::createModel(vector<int>& sequence){
    sizeProblem=sequence.size();
    if (solvers.count(sizeProblem)){
        workspace_ptr=solvers.at(sizeProblem);
        this->sequence=sequence;
        update_b();
    }else{
        SolveSocpSCS::createModel(sequence);
        solvers.insert({sizeProblem,workspace_ptr});
        info_struct.solvers_made++;
    }
}

void SolveSocpSCSWithReuse::clear_removable_constraints(){}//don't delete the old workspaces
void SolveSocpSCSWithReuse::clear_removable_constraints(int prev_pos, int curr_pos){
    clear_removable_constraints();
}

void SolveSocpSCSWithReuse::solveSOCP(){
    SolveSocpSCS::solveSOCP();
    accumulate_info();
}
void SolveSocpSCSWithReuse::accumulate_info(){
    SolveSocpSCS::accumulate_info(info_struct);
    info_struct.solvers_made=solvers.size();
}
void SolveSocpSCSWithReuse::update_b(){
    size_t m=sequence.size()-1;
    size_t nf=m+1;
    size_t nx=SOCP_NDIM*m;
    nvar=nf+nx;
    ncon=m*SOCP_NDIM+m+nf*SOCP_NDIM+nf;
    size_t xstart=nf;
    Eigen::VectorXd b(ncon);
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
    scs_update(workspace_ptr,b.data(),nullptr);
}

/*provide a warm start to SCS
Parameters: primal_guess : double *
                guess for the primal variables. MUST be the correct length (# of columns in constraint matrix). Copied, not mutated.
            slack_guess : double *
                guess for the slack variables. MUST be the correct length (# of rows in constraint matrix). Copied, not mutated.
            dual_guess : double *
                guess for the dual variables. MUST be the correct length (# of rows in constraint matrix). Copied, not muted.
Mutates:    solution
            f_value
            violation
            m_num_solves
            *workspace_ptr
*/
void SolveSocpSCSWithReuse::solve_warm(double* primal_guess, double* slack_guess, double* dual_guess){
    std::memcpy(solution.x,primal_guess,nvar*sizeof(double));
    std::memcpy(solution.s,slack_guess,ncon*sizeof(double));
    std::memcpy(solution.y,dual_guess,ncon*sizeof(double));
    scs_solve(workspace_ptr,&solution,&info,true);
    f_value=info.pobj;
    violation=info.res_pri;
    m_num_solves+=1;
    accumulate_info();
}