#include "SolveSocpSCS.h"

std::ostream& operator<<(std::ostream& os, const ScsInfo& status_obj){
    os << status_obj.status;
    return os;
}

std::ostream& operator<<(std::ostream& os, const SolveSocpSCSStatistics& info_struct){
    cout << "SOCPs solved: "<<info_struct.m_num_solves<<endl;
    cout << "Solvers made: "<<info_struct.solvers_made<<endl;
    cout << "SOCP Internal Iterations: "<<info_struct.iterations<<endl;
    cout << "SOCP Internal Solve Time: "<<info_struct.solve_time<<endl;
    cout << "SOCP Setup Time: "<<info_struct.setup_time<<endl;
    cout << "    SOCP KKT init Time: "<<info_struct.kktinit_time<<endl;
    cout << "SOCP Iteration Time: "<<info_struct.kkt_solve_time<<endl;
    return os;
}
typedef Eigen::Triplet<double> Triplet;

ScsMatrix* no_copy_eigen_SparseMatrix_to_ScsMatrix(Eigen::SparseMatrix<double>& in){
    // in.makeCompressed();
    ScsMatrix* out=new ScsMatrix();
    out->m=in.rows();
    out->n=in.cols();
    out->p=in.outerIndexPtr();
    out->i=in.innerIndexPtr();
    out->x=in.valuePtr();
    return out;
}

SolveSocpSCS::SolveSocpSCS(Data * instance, vector<int>& in_sequence):sequence(in_sequence),objectData(instance),sizeProblem(in_sequence.size()),f_value(0),solution(),workspace_ptr(),settings(){
    scs_set_default_settings(&settings);
    settings.verbose=false;
    //TODO make sure SCS tolerances are comparable to Clarabel
    createModel(sequence);
}

SolveSocpSCS::SolveSocpSCS(Data * instance, int size_seq):sequence(),objectData(instance),sizeProblem(size_seq),f_value(0),solution(),workspace_ptr(),settings(){
    scs_set_default_settings(&settings);
    settings.verbose=false;
}

/*
void SolveSocpSCS::createModel(vector<int>& sequence)

store a ScsWork in this->workspace_ptr and populate its data for the SOCP to visit neighbhoods in the order defined by sequence

Parameters: sequence : vector<int>&
                the order to visit neighborhoods in. The first element must be the depot, of radius 0. The path returns to the depot after visiting sequence[-1]
*/
void SolveSocpSCS::createModel(vector<int>& sequence){
    //decision vector is [f0 f2...fm x10 x11...x1ndim...xmndim]
    size_t m=sequence.size()-1;
    size_t nf=m+1;
    size_t nx=SOCP_NDIM*m;
    nvar=nf+nx;
    size_t xstart=nf;
    sizeProblem=nf;
    this->sequence=sequence;
    //Quadratic cost matrix is empty

    //Unit linear cost on f, 0 cost on xij
    scs_float* q=(scs_float*)calloc(nvar,sizeof(scs_float));
    for(size_t i=0;i<nf;i++){
        q[i]=1;
    }

    //Conic constraints must be expressed as Ax+s=b, s in K
    //first we have m SOC constraints, each with dimension ndim+1, to enforce ||xi-ci||_2<=radii[i]
    //each such constraint requires ndim non-zeros in A
    //then we have nf SOC constraints, each with dimension ndim+1, to enforce ||x(i+1)-xi||<=fi
    //each such constraint requires 1+2*ndim non-zeros in A
    //A is thus (m*ndim+m+nf*ndim+nf,nvar)
    //A will have m*ndim+nf*(1+2*ndim) nonzeros
    ncon=m*SOCP_NDIM+m+nf*SOCP_NDIM+nf;
    Eigen::SparseMatrix<double> A(ncon,nvar);
    scs_float* b=(scs_float*)malloc(sizeof(scs_float)*(ncon));
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(m*SOCP_NDIM+nf*(1+2*SOCP_NDIM));
    size_t row=0;

    /*||xi-ci||_2<=radii[i] thus becomes s==[radii[i]] in K_SOC
                                            [ci-xi   ]
    expressed as [0]xi+s==[radii[i]]=b
                 [I]      [ci      ]*/
    for(size_t i=0;i<m;i++){
        b[row]=objectData->getRadius( sequence[ i+1 ] );//lhs first entry is the radius
        row+=1;//row of zeros
        //x
        tripletList.push_back(Triplet(row,xstart+i*SOCP_NDIM,1));
        b[row]=objectData->getCoordx(sequence[i+1]);
        row+=1;
        //y
        tripletList.push_back(Triplet(row,xstart+i*SOCP_NDIM+1,1));
        b[row]=objectData->getCoordy(sequence[i+1]);
        row+=1;
        //z
        tripletList.push_back(Triplet(row,xstart+i*SOCP_NDIM+2,1));
        b[row]=objectData->getCoordz(sequence[i+1]);
        row+=1;
    }

    /*||x(i+1)j-xij||<=fi becomes s=[fi,x(i+1)-xi] in K_SOC
     expressed as [-1 0  0][fi]+s==0
                  [0  I -I][xi]
                           [x(i+1)]*/
    for(size_t i=0;i<nf;i++){
        tripletList.push_back(Triplet(row,i,-1));//fi
        b[row]=0;
        row+=1;
        if(i==0){
            //depot to first turn
            //||x1-depot||<=f0
            //s=[f0,x1-depot]
            //x
            tripletList.push_back(Triplet(row,xstart,-1));//x10
            b[row]=-objectData->getCoordx(sequence[0]);
            row+=1;
            //y
            tripletList.push_back(Triplet(row,xstart+1,-1));//x11
            b[row]=-objectData->getCoordy(sequence[0]);
            row+=1;
            //z
            tripletList.push_back(Triplet(row,xstart+2,-1));//x12
            b[row]=-objectData->getCoordz(sequence[0]);
            row+=1;

        }else if(i==m){
            //last turn to depot
            //||depot-xm||<=fi
            //s=[fi,depot-xm]
            //x
            tripletList.push_back(Triplet(row,xstart+(m-1)*SOCP_NDIM,1));//xm0
            b[row]=objectData->getCoordx(sequence[0]);
            row+=1;
            //y
            tripletList.push_back(Triplet(row,xstart+(m-1)*SOCP_NDIM+1,1));//xm1
            b[row]=objectData->getCoordy(sequence[0]);
            row+=1;
            //z
            tripletList.push_back(Triplet(row,xstart+(m-1)*SOCP_NDIM+2,1));//xm2
            b[row]=objectData->getCoordz(sequence[0]);
            row+=1;
        }else{
            //||x(i+1)j-xij||<=fi
            //s=[fi,x(i+1)-xi]
            for(size_t j=0;j<SOCP_NDIM;j++){
                tripletList.push_back(Triplet(row,xstart+(i-1)*SOCP_NDIM+j,1));//xij
                tripletList.push_back(Triplet(row,xstart+i*SOCP_NDIM+j,-1));//x(i+1)j
                b[row]=0;
                row+=1;
            }
        }
    }
    A.setFromTriplets(tripletList.cbegin(),tripletList.cend());

    ScsData data;
    data.m=A.rows();
    data.n=A.cols();
    data.A=no_copy_eigen_SparseMatrix_to_ScsMatrix(A);
    data.P=nullptr;
    data.b=b;
    data.c=q;

    //specify the m+nf by (ndim+1) SOCs
    ScsCone cones;
    cones.z=0;
    cones.l=0;
    cones.bsize=0;
    cones.bu=nullptr;
    cones.bl=nullptr;
    cones.qsize=m+nf;
    cones.q=(scs_int*)malloc((m+nf)*sizeof(scs_int));
    for(size_t i=0;i<m+nf;i++){
        cones.q[i]=SOCP_NDIM+1;
    }
    cones.ssize=0;
    cones.s=nullptr;
    cones.ep=0;
    cones.ed=0;
    cones.psize=0;
    cones.p=nullptr;
    workspace_ptr=scs_init(&data,&cones,&settings);
    free(b);
    free(q);
    free(cones.q);
    delete data.A;
}

void SolveSocpSCS::initialize_model(){}//nothing to do, provided to match API of SolveSocpCplex
void SolveSocpSCS::clear_removable_constraints(){
    scs_finish(workspace_ptr);
    sizeProblem=0;
}
void SolveSocpSCS::clear_removable_constraints(int prev_pos, int curr_pos){
    clear_removable_constraints();
}
void SolveSocpSCS::populate_removable_constraints(vector<int>& sequence){
    createModel(sequence);
}
void SolveSocpSCS::populate_removable_constraints(vector<int>& sequence, int prev_pos, int curr_pos){
    populate_removable_constraints(sequence);
}
void SolveSocpSCS::solveSOCP(vector<int>& sequence){//only call this if we used the constructor that doesn't call createModel!
    createModel(sequence);
    solveSOCP();
}

void SolveSocpSCS::solveSOCP(){
    scs_solve(workspace_ptr,&solution,&info,false);
    f_value=info.pobj;
    violation=info.res_pri;
    m_num_solves+=1;
}

void SolveSocpSCS::finishSOCP(){
    scs_finish(workspace_ptr);
    free(solution.x);
    free(solution.y);
    free(solution.s);
    sizeProblem=0;
};

void SolveSocpSCS::accumulate_info(SolveSocpSCSStatistics& info_struct){
    info_struct.m_num_solves++;
    info_struct.solve_time+=get_latest_solve_time();
    info_struct.setup_time+=get_latest_setup_time();
    info_struct.kktinit_time+=get_latest_kktinit_time();
    info_struct.kkt_solve_time+=get_latest_kkt_solve_time();
    info_struct.iterations+=get_latest_iterations();
}

double SolveSocpSCS::getF_value(){
    return f_value;
}

void SolveSocpSCS::printF_value(){
    cout << "Value of objective function: "<< f_value << endl;
}

double SolveSocpSCS::getSolutionX( int idx){
    if(idx==0){
        return objectData->getCoordx(sequence[0]);
    }
    return solution.x[sizeProblem+(idx-1)*SOCP_NDIM];
}
double SolveSocpSCS::getSolutionY( int idx){
    if(idx==0){
        return objectData->getCoordy(sequence[0]);
    }
    return solution.x[sizeProblem+(idx-1)*SOCP_NDIM+1];
}
double SolveSocpSCS::getSolutionZ( int idx){
    if(idx==0){
        return objectData->getCoordz(sequence[0]);
    }
    return solution.x[sizeProblem+(idx-1)*SOCP_NDIM+2];
}

void SolveSocpSCS::printSolution(vector<int>& sequence){
    cout<<"SolverStatus: "<<info.status<<endl;
    cout << "Solution: " << endl;
    for ( int j=0;j<sizeProblem;j++ ){
        cout << "( "<< getSolutionX(j) << ", " << getSolutionY(j) << ", " << getSolutionZ(j) << " )" << endl;
    }
    cout << "( "<< getSolutionX(0) << ", " << getSolutionY(0) << ", " << getSolutionZ(0) << " )" << endl;//print the depot coordinates again at the end
}