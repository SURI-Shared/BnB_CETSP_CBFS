#include "clarabel_interface/SolveRedundantSocpClarabel.h"

typedef Eigen::Triplet<double> Triplet;

SolveRedundantSocpClarabel::SolveRedundantSocpClarabel(Data * instance, int size_seq):SolveSocpClarabel(instance,size_seq){}

SolveRedundantSocpClarabel::SolveRedundantSocpClarabel(Data * instance, int size_seq, bool reduced_first_correction):SolveSocpClarabel(instance,size_seq,reduced_first_correction){}

/*
void SolveRedundantSocpClarabel::createModel(vector<int>& sequence)

store a clarabel::DefaultSolver in this->solver and populate its data for the SOCP to visit neighbhoods in the order defined by sequence

Parameters: sequence : vector<int>&
                the order to visit neighborhoods in. The first element must be the depot, of radius 0. The path returns to the depot after visiting sequence[-1]
*/
void SolveRedundantSocpClarabel::createModel(vector<int>& sequence){
    //decision vector is [f0 f2...fm x11...x1ndim...xmndim w01...w0ndim...wmndim v11...v1ndim...vmndim]
    size_t m=sequence.size()-1;
    size_t nf=m+1;
    size_t nx=SOCP_NDIM*m;
    size_t nw=SOCP_NDIM*nf;
    size_t nv=SOCP_NDIM*m;
    size_t nvar=nf+nx+nw+nv;
    size_t xstart=nf;
    size_t wstart=xstart+nx;
    size_t vstart=wstart+nw;
    sizeProblem=nf;
    this->sequence=sequence;
    //Quadratic cost matrix is empty
    Eigen::SparseMatrix<double> P(nvar,nvar);

    //Unit linear cost on f, 0 cost on xij, wij, vij
    Eigen::VectorXd q(nvar);
    for(size_t i=0;i<nf;i++){
        q[i]=1;
    }
    for(size_t j=nf;j<nvar;j++){
        q[j]=0;
    }
    //Conic constraints must be expressed as Ax+s=b, s in K
    //
    //Affine equality constraints, with K=zero cone
    //w
    //first we have 1 affine equality constraint, with dimension ndim, to enforce w0-x1=-depot
    ///requiring 2*ndim non-zeros in A
    //then we have m-1 affine equality constraints, with dimension ndim, to enforce wi+xi-x(i+1)=0
    //requiring (m-1)*3*ndim non-zeros in A
    //then we have 1 affine equality constraint, with dimension ndim, to enforce wm+xm=depot
    //requiring 2*ndim non-zeros in A
    //v
    //m affine equality constraints, with dimension ndim, to enforce vi+xi=center[i]
    //requiring m*2*ndim non-zeros in A
    //
    //SOC constraints, with K=second order cone
    //v
    //first we have m SOC constraints, each with dimension ndim+1, to enforce ||vi||_2<=radii[i]
    //requiring m*ndim non-zeros in A
    //w
    //then we have nf SOC constraints, each with dimension ndim+1, to enforce ||wi||<=fi
    //requring nf*(1+ndim) non-zeros in A
    //A is thus (nf*ndim+m*ndim+m*(ndim+1)+nf*(ndim+1),nvar)
    //A will have 2*ndim+(m-1)*3*ndim+2*ndim+m*2*ndim+m*ndim+nf*(1+ndim) nonzeros
    //and there are, in order, nf+m zero cones of dimension ndim, followed by m+nf SOC cones of dimension 1+ndim
    vector<clarabel::SupportedConeT<double>> cones;cones.reserve(nf+m+m+nf);
    Eigen::SparseMatrix<double> A(nf*SOCP_NDIM+m*SOCP_NDIM+m*(SOCP_NDIM+1)+nf*(SOCP_NDIM+1),nvar);
    Eigen::VectorXd b(nf*SOCP_NDIM+m*SOCP_NDIM+m*(SOCP_NDIM+1)+nf*(SOCP_NDIM+1));
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(2*SOCP_NDIM+(m-1)*3*SOCP_NDIM+2*SOCP_NDIM+m*2*SOCP_NDIM+m*SOCP_NDIM+nf*(1+SOCP_NDIM));
    size_t row=0;

    //w0-x1==-depot
    //x
    tripletList.push_back(Triplet(row,wstart,1));//w00
    tripletList.push_back(Triplet(row,xstart,-1));//x10
    b[row]=-objectData->getCoordx(sequence[0]);
    row+=1;
    //y
    tripletList.push_back(Triplet(row,wstart+1,1));//w01
    tripletList.push_back(Triplet(row,xstart+1,-1));//x11
    b[row]=-objectData->getCoordy(sequence[0]);
    row+=1;
    //z
    tripletList.push_back(Triplet(row,wstart+2,1));//w02
    tripletList.push_back(Triplet(row,xstart+2,-1));//x12
    b[row]=-objectData->getCoordz(sequence[0]);
    row+=1;
    cones.push_back(clarabel::ZeroConeT<double>(SOCP_NDIM));

    //wi+xi-x(i+1)=0
    for(size_t i=1;i<m;i++){
        for(size_t j=0;j<SOCP_NDIM;j++){
            tripletList.push_back(Triplet(row,wstart+i*SOCP_NDIM+j,1));//wij
            tripletList.push_back(Triplet(row,(xstart-SOCP_NDIM)+i*SOCP_NDIM+j,1));//xij. Note that xstart points to x1, not x0, hence the -SOCP_NDIM.
            tripletList.push_back(Triplet(row,(xstart-SOCP_NDIM)+(i+1)*SOCP_NDIM+j,-1));//x(i+1)j
            b[row]=0;
            row+=1;
        }
        cones.push_back(clarabel::ZeroConeT<double>(SOCP_NDIM));
    }

    //wm+xm=depot
    //x
    tripletList.push_back(Triplet(row,wstart+m*SOCP_NDIM,1));//wm0
    tripletList.push_back(Triplet(row,xstart+(m-1)*SOCP_NDIM,1));//xm0. Note that xstart points to x1, not x0, hence the (m-1).
    b[row]=objectData->getCoordx(sequence[0]);
    row+=1;
    //y
    tripletList.push_back(Triplet(row,wstart+m*SOCP_NDIM+1,1));//wm1
    tripletList.push_back(Triplet(row,xstart+(m-1)*SOCP_NDIM+1,1));//xm1
    b[row]=objectData->getCoordy(sequence[0]);
    row+=1;
    //z
    tripletList.push_back(Triplet(row,wstart+m*SOCP_NDIM+2,1));//wm2
    tripletList.push_back(Triplet(row,xstart+(m-1)*SOCP_NDIM+2,1));//xm2
    b[row]=objectData->getCoordz(sequence[0]);
    row+=1;
    cones.push_back(clarabel::ZeroConeT<double>(SOCP_NDIM));

    //vi+xi=center[i]
    for(size_t i=0;i<m;i++){
        //x
        tripletList.push_back(Triplet(row,vstart+i*SOCP_NDIM,1));//vi0
        tripletList.push_back(Triplet(row,xstart+i*SOCP_NDIM,1));//xi0
        b[row]=objectData->getCoordx(sequence[i+1]);
        row+=1;

        //x
        tripletList.push_back(Triplet(row,vstart+i*SOCP_NDIM+1,1));//vi1
        tripletList.push_back(Triplet(row,xstart+i*SOCP_NDIM+1,1));//xi1
        b[row]=objectData->getCoordy(sequence[i+1]);
        row+=1;

        //x
        tripletList.push_back(Triplet(row,vstart+i*SOCP_NDIM+2,1));//vi2
        tripletList.push_back(Triplet(row,xstart+i*SOCP_NDIM+2,1));//xi2
        b[row]=objectData->getCoordz(sequence[i+1]);
        row+=1;

        cones.push_back(clarabel::ZeroConeT<double>(SOCP_NDIM));
    }

    /*||vi||_2<=radii[i] thus becomes s==[radii[i]] in K_SOC
                                            [vi   ]
    expressed as [ 0]vi+s==[radii[i]]=b
                 [-I]      [0      ]                    */
    for(size_t i=0;i<m;i++){
        b[row]=objectData->getRadius( sequence[ i+1 ] );//rhs first entry is the radius
        row+=1;//row of zeros
        //x
        tripletList.push_back(Triplet(row,vstart+i*SOCP_NDIM,-1));
        b[row]=0;
        row+=1;
        //y
        tripletList.push_back(Triplet(row,vstart+i*SOCP_NDIM+1,-1));
        b[row]=0;
        row+=1;
        //z
        tripletList.push_back(Triplet(row,vstart+i*SOCP_NDIM+2,-1));
        b[row]=0;
        row+=1;

        cones.push_back(clarabel::SecondOrderConeT<double>(SOCP_NDIM+1));
    }

    /*||wi||<=fi becomes s=[fi,wi] in K_SOC
     expressed as [-1  0][fi]+s==0
                  [0  -I][wi]               */
    for(size_t i=0;i<nf;i++){
        tripletList.push_back(Triplet(row,i,-1));//fi
        b[row]=0;
        row+=1;
        for(size_t j=0;j<SOCP_NDIM;j++){
            tripletList.push_back(Triplet(row,wstart+i*SOCP_NDIM+j,-1));//wij
            b[row]=0;
            row+=1;
        }
        cones.push_back(clarabel::SecondOrderConeT<double>(SOCP_NDIM+1));
    }
    A.setFromTriplets(tripletList.cbegin(),tripletList.cend());

    solver_ptr=new clarabel::DefaultSolver<double>(P,Eigen::Ref<Eigen::VectorXd>(q),A,Eigen::Ref<Eigen::VectorXd>(b),cones,settings);
}