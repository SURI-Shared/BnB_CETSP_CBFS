#ifndef TSP_LB_H
#define TSP_LB_H
#define IL_STD 1

#include <stdio.h>
#include <ilcplex/ilocplex.h>
#include <pthread.h>
#define new concorde_new//workaround for use of new and class as argument names
#define class concorde_class
extern "C" {
#include <concorde.h>
}
#undef new
#undef class
#include "subgraph.h"

using namespace std;

class subgraph_tsp_solver
{
public:
    subgraph_tsp_solver() {}
    subgraph_tsp_solver(subgraph* graph) : m_subgraph(graph)
    {
        m_size_subgraph = m_subgraph->get_size_subgraph();
        m_sequence.resize(m_size_subgraph, 0);
        m_precedence_sequence = m_subgraph->m_precedence_sequence;
        //m_int_factor = m_subgraph->get_int_factor();
        //m_dist_addon = m_subgraph->get_dist_addon();
    }
    subgraph_tsp_solver(Data* d) : m_data(d)
    {
        m_size_subgraph = m_data->getSizeInst();
    }

    void solve_cplex_miqcp();
    void solve_cplex_weak();
    void solve_concorde();
    void solve_concorde(vector<vector<double>>& dist_mat);
    double get_solution() { return m_solution; }
    void set_int_factor(double int_factor) { m_int_factor = int_factor; }
    void set_dist_addon(int addon) { m_dist_addon = addon; }
    void set_time_lim(double time_lim) { m_time_lim = time_lim; }
    void get_tour_sequence(vector<int>& sequence);
    void set_init_solution(vector<int>& init_seq);

    vector<int> m_precedence_sequence;
    vector<int> m_init_sequence;
    vector<int> m_sequence;
    double m_solution = 0;
private:
    Data* m_data;
    int m_size_subgraph = 0;
    subgraph* m_subgraph;

    double m_int_factor = 1.0;
    int m_dist_addon = 0;
    double m_time_lim = 100000;
};

#endif