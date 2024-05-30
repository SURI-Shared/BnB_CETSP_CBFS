#ifndef WARM_START_SOCP_H
#define WARM_START_SOCP_H

#include <vector>
#include <Clarabel>

#include"Data.h"

class BnBNodeForWarmStart{
    public:
        BnBNodeForWarmStart(const std::vector<int>&,const Eigen::Map<Eigen::VectorX<double>>& primals,const Eigen::Map<Eigen::VectorX<double>>& duals,const Eigen::Map<Eigen::VectorX<double>>& slacks);

        const std::vector<int> sequence;
        const Eigen::VectorX<double> primals;
        const Eigen::VectorX<double> duals;
        const Eigen::VectorX<double> slacks;

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
        void construct_initial_guess(const std::vector<int>& current_sequence, const std::vector<int>& parent_sequence,
                                               const Eigen::Ref<const Eigen::VectorX<double>>& parent_primal, const Eigen::Ref<const Eigen::VectorX<double>>& parent_slack, const Eigen::Ref<const Eigen::VectorX<double>>& parent_dual,
                                               Data* instanceData,std::vector<double>& out_primals,std::vector<double>& out_slacks, std::vector<double>& out_duals);
        void solve_insertion_problem(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&,double, double* fa_out, double* fb_out, double* x_out, double* in_neighborhood_slack_out, double* fa_slack_out, double* fb_slack_out, double* in_neighborhood_dual_out, double* fa_dual_out,double* fb_dual_out);
        double get_init_time(){return init_time;}
        double get_construct_time(){return construct_time;}
        double get_solve_time(){return solve_time;}
        double get_total_time(){return init_time+construct_time;}
    protected:
        clarabel::DefaultSolver<double>* insertion_problem_solver_ptr;
        double init_time;
        double construct_time;
        double solve_time;

};

#endif