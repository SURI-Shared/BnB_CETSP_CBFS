#ifndef LOCAL_SEARCH_H
#define LOCAL_SEARCH_H

#include "SolveSocpCplex.h"
#include "subgraph.h"
#include "tsp_lb.h"
#include "util.h"

using namespace std;

class Local_search_cetsp
{
public:
    Local_search_cetsp(Data* d) : m_data(d)
    {}
    void initialize_sequence(vector<int>& seq, vector<double>& coord_x, vector<double>& coord_y, vector<double>& coord_z);
    void solve();
    double solve_new_solution();
    
    double m_init_solution = 0;
    vector<int> m_sequence;
private:
    struct coordinates 
    {
        double x;
        double y;
        double z;
    };
    void update_dist_mat(vector<double>& coord_x, vector<double>& coord_y, vector<double>& coord_z);
    double euclidianDistance(coordinates* p1, coordinates* p2);

    Data* m_data;
    int m_size_graph = 0;
    vector<vector<double>> m_dist_mat;
    
    int m_iter_lim = 100000;
    double m_time_lim = 360000.0;
    double m_incumbent_sol = 10000000.0;
    double m_equal_tol = 0.000001;
    vector<int> m_incumbent_sequence;
};

#endif