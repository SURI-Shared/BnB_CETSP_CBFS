#include "subgraph.h"

subgraph::subgraph(Data* data, vector<bool>& nodes, double int_factor)
{
    m_size_subgraph = 0;
    m_data = data;
    m_int_factor = int_factor;
    m_dist_addon = 0;
    int sizeInst = m_data->getSizeInst();
    
    for (int i = 0; i < sizeInst; i++)
    {
        if (nodes[i])
        {
            m_tour_id_to_graph_id[i] = m_size_subgraph;
            m_index_subgraph_to_relaxed_dist.push_back(i);
            m_size_subgraph++;
        }
    }

    m_subgraph_dist_mat.clear();
    m_subgraph_dist_mat.resize(m_size_subgraph);
    for (int i = 0; i < m_size_subgraph; i++)
    {
        m_subgraph_dist_mat[i].clear();
        m_subgraph_dist_mat[i].resize(m_size_subgraph, 0);
    }
    for (int i = 0; i < m_size_subgraph; i++)
    {
        for (int j = i+1; j < m_size_subgraph; j++)
        {
            m_subgraph_dist_mat[i][j] = m_data->get_relaxed_dist(m_index_subgraph_to_relaxed_dist[i], m_index_subgraph_to_relaxed_dist[j]);
            m_subgraph_dist_mat[j][i] =  m_subgraph_dist_mat[i][j];
        }
    }
    m_precedence_sequence.resize(1,0);
}

subgraph::subgraph(Data* data)
{
    m_data = data;
    m_size_subgraph = m_data->getSizeInst();

    m_subgraph_dist_mat.clear();
    m_subgraph_dist_mat.resize(m_size_subgraph);
    for (int i = 0; i < m_size_subgraph; i++)
    {
        m_subgraph_dist_mat[i].clear();
        m_subgraph_dist_mat[i].resize(m_size_subgraph, 0);

        m_tour_id_to_graph_id[i] = i;
        m_index_subgraph_to_relaxed_dist.push_back(i);
    }
    for (int i = 0; i < m_size_subgraph; i++)
    {
        for (int j = i+1; j < m_size_subgraph; j++)
        {
            m_subgraph_dist_mat[i][j] = m_data->get_relaxed_dist(m_index_subgraph_to_relaxed_dist[i], m_index_subgraph_to_relaxed_dist[j]);
            m_subgraph_dist_mat[j][i] =  m_subgraph_dist_mat[i][j];
        }
    }
    m_precedence_sequence.resize(1,0);
}

subgraph::subgraph(vector<vector<double>>& dist_mat)
{
    m_size_subgraph = dist_mat.size();

    m_subgraph_dist_mat.clear();
    m_subgraph_dist_mat.resize(m_size_subgraph);
    for (int i = 0; i < m_size_subgraph; i++)
    {
        m_subgraph_dist_mat[i].clear();
        m_subgraph_dist_mat[i].resize(m_size_subgraph, 0);
        m_tour_id_to_graph_id[i] = i;
        m_index_subgraph_to_relaxed_dist.push_back(i);
    }
    for (int i = 0; i < m_size_subgraph; i++)
    {
        for (int j = i+1; j < m_size_subgraph; j++)
        {
            m_subgraph_dist_mat[i][j] = dist_mat[i][j];
            m_subgraph_dist_mat[j][i] =  m_subgraph_dist_mat[i][j];
        }
    }
    m_precedence_sequence.resize(1,0);
}

double subgraph::get_relaxed_dist(int i, int j)
{
    if (m_special_depot)
    {
        // if (m_special_depot_all_zero && j == m_size_subgraph - 1)
        // {
        //     // No multiplier and no addon value
        //     return m_subgraph_dist_mat[i][j];
        // }
        // else if (j == m_size_subgraph - 1)
        // {
        //     if (i == 0 || m_index_subgraph_to_relaxed_dist[i] == m_special_to_index)
        //     {
        //         // No multiplier and no addon value
        //         return m_subgraph_dist_mat[i][j];
        //     }
        //     else
        //     {
        //         return m_subgraph_dist_mat[i][j] * m_int_factor + m_dist_addon;
        //     }
        // }
        // return m_subgraph_dist_mat[i][j] * m_int_factor + m_dist_addon;
    }
    else
    {
        return m_subgraph_dist_mat[i][j] * m_int_factor + m_dist_addon;
    }
}

int subgraph::get_size_subgraph()
{
    return m_size_subgraph;
}

// void subgraph::set_special_to_index(int i)
// {
//     m_special_to_index = i;

//     vector<double> temp_dists(m_size_subgraph + 1, m_large_dist_constant);

//     temp_dists[0] = 0;

//     temp_dists[m_index_subgraph_to_relaxed_dist[i]] = 0;
    
//     temp_dists[m_size_subgraph] = 0;

//     add_vertex_to_graph(temp_dists);
//     m_special_depot = true;
// }

// void subgraph::set_special_depot_all_zero()
// {
//     vector<double> temp_dists(m_size_subgraph + 1, 0);

//     add_vertex_to_graph(temp_dists);
//     m_special_depot = true;
//     m_special_depot_all_zero = true;
// }

void subgraph::add_vertex_to_graph(vector<double> &dist_to_vertices)
{
    m_subgraph_dist_mat.push_back(dist_to_vertices);
    for (int i = 0; i < m_size_subgraph; i++)
    {
        m_subgraph_dist_mat[i].push_back(dist_to_vertices[i]);
    }
    m_size_subgraph++;
}

void subgraph::set_precedence_constraint(int i, int j)
{
    m_subgraph_dist_mat[m_index_subgraph_to_relaxed_dist[i]][m_index_subgraph_to_relaxed_dist[j]] = m_large_dist_constant;
}
 
void subgraph::convert_atsp_to_tsp()
{
    int org_size_subgraph = m_size_subgraph;
    m_size_subgraph = m_size_subgraph * 2;
    vector<vector<double>> temp_dist_mat(m_size_subgraph);

    for (int i = 0; i < m_size_subgraph; i++)
    {
        temp_dist_mat[i].resize(m_size_subgraph);
        for (int j = 0; j < m_size_subgraph; j++)
        {
            if (i == j)
            {
                temp_dist_mat[i][j] = 0;
            }
            else if (i < org_size_subgraph && j < org_size_subgraph)
            {
                temp_dist_mat[i][j] = m_large_dist_constant;
            }
            else if (i < org_size_subgraph && j >= org_size_subgraph)
            {
                if (i == j - org_size_subgraph)
                    temp_dist_mat[i][j] = 0;
                else
                    temp_dist_mat[i][j] = m_subgraph_dist_mat[j - org_size_subgraph][i] + m_atsp_to_tsp_addon;
            }
            else if (i >= org_size_subgraph && j < org_size_subgraph)
            {
                if (i - org_size_subgraph == j)
                    temp_dist_mat[i][j] = 0;
                else
                    temp_dist_mat[i][j] = m_subgraph_dist_mat[i - org_size_subgraph][j] + m_atsp_to_tsp_addon;
            }
            else if (i >= org_size_subgraph && j >= org_size_subgraph)
            {
                temp_dist_mat[i][j] = m_large_dist_constant;
            }
        }
    }
    m_subgraph_dist_mat = temp_dist_mat;
}