#ifndef WARM_START_SOCP_H
#define WARM_START_SOCP_H

#include <vector>
#include <Clarabel>

#include"Data.h"

class BnBNodeForWarmStart{
    public:
        BnBNodeForWarmStart(BnBNodeForWarmStart*,const std::vector<int>&,const Eigen::Map<Eigen::VectorX<double>>& primals,const Eigen::Map<Eigen::VectorX<double>>& duals);

        BnBNodeForWarmStart* parent;//NULL if no parent
        const std::vector<int> sequence;
        const Eigen::VectorX<double> primals;
        const Eigen::VectorX<double> duals;

        std::vector<Eigen::Vector3d> turning_points() const;

};

class WarmStartHandler{
    public:
        WarmStartHandler();
        ~WarmStartHandler();
        void construct_initial_guess(const std::vector<int>& current_sequence, 
                             const std::vector<int>& parent_sequence, const std::vector<Eigen::Vector3d> & parent_turning_points,  const Eigen::Ref<const Eigen::VectorX<double>>& parent_dual,
                             Data* instanceData,std::vector<double>& out_primals,std::vector<double>& out_slacks, std::vector<double>& out_duals);
        void construct_initial_guess(const std::vector<int>& current_sequence,const BnBNodeForWarmStart& parent,Data* instanceData,std::vector<double>& out_primals,std::vector<double>& out_slacks, std::vector<double>& out_duals);
        void solve_insertion_problem(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&,double,double*,double*);
    protected:
        clarabel::DefaultSolver<double>* insertion_problem_solver_ptr;

};

#endif