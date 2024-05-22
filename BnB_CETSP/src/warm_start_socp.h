#ifndef WARM_START_SOCP_H
#define WARM_START_SOCP_H

#include <vector>
#include <Clarabel>

#include"Data.h"

class WarmStartHandler{
    public:
        WarmStartHandler();
        ~WarmStartHandler();
        void construct_initial_guess(const std::vector<int>& current_sequence, 
                             const std::vector<int>& parent_sequence, const Eigen::ArrayX3d & parent_turning_points, Eigen::Map<Eigen::VectorX<double>> parent_dual,
                             Data* instanceData,std::vector<double>& out_primals,std::vector<double>& out_slacks, std::vector<double>& out_duals);
        void solve_insertion_problem(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&,double,double*,double*);
    protected:
        clarabel::DefaultSolver<double>* insertion_problem_solver_ptr;

};

#endif