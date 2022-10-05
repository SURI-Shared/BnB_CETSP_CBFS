#ifndef SUBGRAPH_H
#define SUBGRAPH_H

#include <unordered_map>
#include <set>
#include "Data.h"

using namespace std;

class subgraph
{
public:
    subgraph(Data* data, vector<bool>& nodes, double int_factor);
    subgraph(Data* data);
    subgraph(vector<vector<double>>& dist_mat);

    double get_relaxed_dist(int i, int j);
    int get_size_subgraph();
    void set_int_factor(double i) { m_int_factor = i; }
    void set_dist_addon(int i) { m_dist_addon = i; }
    void set_atsp_to_tsp_addon(int i) { m_atsp_to_tsp_addon = i; }
    void set_special_to_index(int i);
    void set_special_depot_all_zero();
    double get_int_factor() { return m_int_factor; }
    int get_dist_addon() { return m_dist_addon; }
    void add_vertex_to_graph(vector<double> &dist_to_vertices);
    void set_precedence_constraint(int i, int j);
    void convert_atsp_to_tsp();

    vector<int> m_precedence_sequence;

private:
    Data* m_data;
    int m_size_subgraph = 0;
    vector<int> m_index_subgraph_to_relaxed_dist;
    vector<vector<double>> m_subgraph_dist_mat;
    unordered_map<int, int> m_tour_id_to_graph_id;

    double m_int_factor = 1.0;
    int m_dist_addon = 0;
    int m_atsp_to_tsp_addon = 1;
    int m_special_to_index = -1;
    bool m_special_depot = false;
    bool m_special_depot_all_zero = false;

    double m_large_dist_constant = 1000.0;

    bool m_use_custom_data_mat = false;
};

#endif