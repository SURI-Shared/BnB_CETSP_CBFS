#include "SolveSocpClarabel.h"

std::ostream& operator<<(std::ostream& os, const clarabel::SolverStatus& status_obj){
    switch(status_obj){
        case clarabel::SolverStatus::Unsolved: 
            os<<"Unsolved";
            break;
        case clarabel::SolverStatus::Solved:
            os<<"Solved";
            break;
        case clarabel::SolverStatus::PrimalInfeasible:
            os<<"PrimalInfeasible";
            break;
        case clarabel::SolverStatus::DualInfeasible:
            os<<"DualInfeasible";
            break;
        case clarabel::SolverStatus::AlmostSolved:
            os<<"AlmostSolved";
            break;
        case clarabel::SolverStatus::AlmostPrimalInfeasible:
            os<<"AlmostPrimalInfeasible";
            break;
        case clarabel::SolverStatus::AlmostDualInfeasible:
            os<<"AlmostDualInfeasible";
            break;
        case clarabel::SolverStatus::MaxIterations:
            os<<"MaxIterations";
            break;
        case clarabel::SolverStatus::MaxTime:
            os<<"MaxTime";
            break;
        case clarabel::SolverStatus::NumericalError:
            os<<"NumericalError";
            break;
        case clarabel::SolverStatus::ScalingError:
            os<<"ScalingError";
            break;
        case clarabel::SolverStatus::InsufficientProgress:
            os<<"InsufficientProgress";
            break;
        default:
            os<<"UnrecognizedStatus";
            break;
    }
    return os;
}
typedef Eigen::Triplet<double> Triplet;

SolveSocpClarabel::SolveSocpClarabel(Data * instance, vector<int>& in_sequence):sequence(in_sequence),objectData(instance),sizeProblem(in_sequence.size()),f_value(0),solution(),solver_ptr(){
    settings=clarabel::DefaultSettings<double>::default_settings();
    settings.verbose=false;
    settings.presolve_enable=false;
    createModel(sequence);
}

SolveSocpClarabel::SolveSocpClarabel(Data * instance, int size_seq):sequence(),objectData(instance),sizeProblem(size_seq),f_value(0),solution(),solver_ptr(){
    settings=clarabel::DefaultSettings<double>::default_settings();
    settings.verbose=false;
    settings.presolve_enable=false;
}

/*
void SolveSocpClarabel::createModel(vector<int>& sequence)

store a clarabel::DefaultSolver in this->solver and populate its data for the SOCP to visit neighbhoods in the order defined by sequence

Parameters: sequence : vector<int>&
                the order to visit neighborhoods in. The first element must be the depot, of radius 0. The path returns to the depot after visiting sequence[-1]
*/
void SolveSocpClarabel::createModel(vector<int>& sequence){
    //decision vector is [f0 f2...fm x10 x11...x1ndim...xmndim]
    size_t m=sequence.size()-1;
    size_t nf=m+1;
    size_t nx=SOCP_NDIM*m;
    size_t nvar=nf+nx;
    size_t xstart=nf;
    sizeProblem=nf;
    this->sequence=sequence;
    //Quadratic cost matrix is empty
    Eigen::SparseMatrix<double> P(nvar,nvar);

    //Unit linear cost on f, 0 cost on xij
    Eigen::VectorXd q(nvar);
    for(size_t i=0;i<nf;i++){
        q[i]=1;
    }
    for(size_t j=nf;j<nvar;j++){
        q[j]=0;
    }

    //Conic constraints must be expressed as Ax+s=b, s in K
    //first we have m SOC constraints, each with dimension ndim+1, to enforce ||xi-ci||_2<=radii[i]
    //each such constraint requires ndim non-zeros in A
    //then we have nf SOC constraints, each with dimension ndim+1, to enforce ||x(i+1)-xi||<=fi
    //each such constraint requires 1+2*ndim non-zeros in A
    //A is thus (m*ndim+m+nf*ndim+nf,nvar)
    //A will have m*ndim+nf*(1+2*ndim) nonzeros
    Eigen::SparseMatrix<double> A(m*SOCP_NDIM+m+nf*SOCP_NDIM+nf,nvar);
    Eigen::VectorXd b(m*SOCP_NDIM+m+nf*SOCP_NDIM+nf);
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

    //specify the m+nf by (ndim+1) SOCs
    vector<clarabel::SupportedConeT<double>> cones;cones.reserve(m+nf);
    for(size_t i=0;i<m+nf;i++){
        cones.push_back(clarabel::SecondOrderConeT<double>(SOCP_NDIM+1));
    }
    solver_ptr=new clarabel::DefaultSolver<double>(P,Eigen::Ref<Eigen::VectorXd>(q),A,Eigen::Ref<Eigen::VectorXd>(b),cones,settings);
}

void SolveSocpClarabel::initialize_model(){}//nothing to do, provided to match API of SolveSocpCplex
void SolveSocpClarabel::clear_removable_constraints(){
    delete solver_ptr;
    solver_ptr=NULL;
    sizeProblem=0;
}
void SolveSocpClarabel::clear_removable_constraints(int prev_pos, int curr_pos){
    clear_removable_constraints();
}
void SolveSocpClarabel::populate_removable_constraints(vector<int>& sequence){
    createModel(sequence);
}
void SolveSocpClarabel::populate_removable_constraints(vector<int>& sequence, int prev_pos, int curr_pos){
    populate_removable_constraints(sequence);
}
void SolveSocpClarabel::solveSOCP(vector<int>& sequence){//only call this if we used the constructor that doesn't call createModel!
    createModel(sequence);
    SolveSocpClarabel::solveSOCP();
}

void SolveSocpClarabel::solveSOCP(){
    solver_ptr->solve();
    new (&solution) clarabel::DefaultSolution<double>(solver_ptr->solution());
    f_value=solution.obj_val;
    violation=solution.r_prim;
    m_num_solves+=1;
}

void SolveSocpClarabel::finishSOCP(){};//nothing to do, provided to match API of SolveSocpCplex

double SolveSocpClarabel::getF_value(){
    return f_value;
}

void SolveSocpClarabel::printF_value(){
    cout << "Value of objective function: "<< f_value << endl;
}

double SolveSocpClarabel::getSolutionX( int idx){
    if(idx==0){
        return objectData->getCoordx(sequence[0]);
    }
    return solution.x[sizeProblem+(idx-1)*SOCP_NDIM];
}
double SolveSocpClarabel::getSolutionY( int idx){
    if(idx==0){
        return objectData->getCoordy(sequence[0]);
    }
    return solution.x[sizeProblem+(idx-1)*SOCP_NDIM+1];
}
double SolveSocpClarabel::getSolutionZ( int idx){
    if(idx==0){
        return objectData->getCoordz(sequence[0]);
    }
    return solution.x[sizeProblem+(idx-1)*SOCP_NDIM+2];
}

void SolveSocpClarabel::printSolution(vector<int>& sequence){
    cout<<"SolverStatus: "<<solution.status<<endl;
    cout << "Solution: " << endl;
    for ( int j=0;j<sizeProblem;j++ ){
        cout << "( "<< getSolutionX(j) << ", " << getSolutionY(j) << ", " << getSolutionZ(j) << " )" << endl;
    }
    cout << "( "<< getSolutionX(0) << ", " << getSolutionY(0) << ", " << getSolutionZ(0) << " )" << endl;//print the depot coordinates again at the end
}