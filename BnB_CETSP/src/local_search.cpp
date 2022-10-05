#include "local_search.h"

void Local_search_cetsp::solve()
{
    // initialize sequence
    if (m_sequence.empty())
    {
        // for (int i = 0; i < m_size_graph; i++)
        // {
        //     m_sequence.push_back(i);
        // }
        subgraph* temp_graph = new subgraph(m_data);
        subgraph_tsp_solver* temp_tsp_solve = new subgraph_tsp_solver(temp_graph);
        temp_tsp_solve->set_int_factor(10000);
        temp_tsp_solve->solve_concorde();
        temp_tsp_solve->get_tour_sequence(m_sequence);
    }

    int num_iter = 0;
    double init_time = cpuTime();
    double temp_sol = 0;
    double pre_sol = -1;

    vector<double> temp_x(m_size_graph, 0);
    vector<double> temp_y(m_size_graph, 0);
    vector<double> temp_z(m_size_graph, 0);

    while (num_iter < m_iter_lim && cpuTime() - init_time < m_time_lim)
    {
        for (int i = 0; i < m_sequence.size(); i++)
        {
            cout << m_sequence[i] << " ";
        }
        cout << endl;
        SolveSocpCplex *solve_socp = new SolveSocpCplex(m_data, m_sequence);

        solve_socp->solveSOCP(m_sequence);

        temp_sol = solve_socp->getF_value();
        solve_socp->finishSOCP();
        if (temp_sol < m_incumbent_sol - m_equal_tol)
        {
            m_incumbent_sol = temp_sol;
            m_incumbent_sequence = m_sequence;
        }
        else if (temp_sol <= pre_sol + m_equal_tol && temp_sol >= pre_sol - m_equal_tol)
        {
            break;
        }
        else
        {
            pre_sol = temp_sol;
        }

        for ( int i = 0; i < m_size_graph; i++ ){
            temp_x[i] = solve_socp->getSolutionX(i);
            temp_y[i] = solve_socp->getSolutionY(i);
            temp_z[i] = solve_socp->getSolutionZ(i);
        }

        update_dist_mat(temp_x, temp_y, temp_z);

        subgraph* temp_graph = new subgraph(m_dist_mat);
        subgraph_tsp_solver* temp_tsp_solve = new subgraph_tsp_solver(temp_graph);
        temp_tsp_solve->set_int_factor(10000);
        temp_tsp_solve->solve_concorde();
        temp_tsp_solve->get_tour_sequence(m_sequence);
        
        num_iter++;
        cout << "Iter: " << num_iter << ", Current_sol: " << temp_sol << ", Incumbent_sol: " << m_incumbent_sol << endl;
    }
}

double Local_search_cetsp::solve_new_solution()
{
    // subgraph* temp_graph = new subgraph(m_dist_mat);
    subgraph_tsp_solver* temp_tsp_solve = new subgraph_tsp_solver();
    temp_tsp_solve->set_int_factor(100000);
    vector<double> tempX(m_size_graph);
    vector<double> tempY(m_size_graph);
    vector<double> tempZ(m_size_graph);
    double temp_sol = 0;
    // temp_tsp_solve->set_init_solution(m_sequence);
    while (true)
    {
        temp_tsp_solve->solve_concorde(m_dist_mat);
        if (temp_tsp_solve->m_sequence.empty()) break;
        if (temp_tsp_solve->m_solution >= m_init_solution - 0.01) break;

        vector<int> concorde_sequence;
        temp_tsp_solve->get_tour_sequence(concorde_sequence);

        // m_sequence = concorde_sequence;
        vector<int> temp_sequence(m_size_graph);
        for (int i = 0; i < m_size_graph; i++)
        {
            temp_sequence[i] = m_sequence[concorde_sequence[i]];
        }
        SolveSocpCplex *solve_socp = new SolveSocpCplex(m_data, temp_sequence);
        solve_socp->solveSOCP(temp_sequence);
        temp_sol = solve_socp->getF_value();
        solve_socp->finishSOCP();

        for ( int i = 0; i < m_size_graph; i++ ){
            tempX[i] = solve_socp->getSolutionX(i);
            tempY[i] = solve_socp->getSolutionY(i);
            tempZ[i] = solve_socp->getSolutionZ(i);
        }
        update_dist_mat(tempX, tempY, tempZ);
        delete solve_socp;
        m_sequence = temp_sequence;
        m_init_solution = temp_sol;
    }

    delete temp_tsp_solve;
    return temp_sol;
}

void Local_search_cetsp::initialize_sequence(vector<int>& seq, vector<double>& coord_x, vector<double>& coord_y, vector<double>& coord_z)
{
    m_sequence = seq;
    m_size_graph = m_sequence.size();
    if (m_dist_mat.size() != m_size_graph)
    {
        m_dist_mat.resize(m_size_graph);
        for (int i = 0; i < m_size_graph; i++)
        {
            m_dist_mat[i].resize(m_size_graph, 0);
        }
    }
    update_dist_mat(coord_x, coord_y, coord_z);
}

void Local_search_cetsp::update_dist_mat(vector<double>& coord_x, vector<double>& coord_y, vector<double>& coord_z)
{
    coordinates p1, p2;
    //double temp_dist;
    for (int i = 0; i < m_size_graph; i++)
    {
        p1.x = coord_x[i]; p1.y = coord_y[i]; p1.z = coord_z[i];
        for (int j = i+1; j < m_size_graph; j++)
        {
            p2.x = coord_x[j]; p2.y = coord_y[j]; p2.z = coord_z[j];
            m_dist_mat[i][j] = euclidianDistance(&p1, &p2);
            m_dist_mat[j][i] = m_dist_mat[i][j];
        }
    }
}

double Local_search_cetsp::euclidianDistance(coordinates* p1, coordinates* p2)
{
   return sqrt( pow( p2->x - p1->x, 2) + pow( p2->y - p1->y, 2) + pow( p2->z - p1->z, 2) );
} 