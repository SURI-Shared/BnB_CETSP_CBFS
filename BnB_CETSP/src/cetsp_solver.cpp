#include "cetsp_solver.h"

CETSP_solver::CETSP_solver(Data* d, CETSP_Options& option)
{
    m_data = d;

    m_overlap_ratio = option.overlap_ratio;
    m_branching_rule = option.branching_rule;
    m_branching_strategy = option.branching_strategy;
    m_strong_branching_size = option.strong_branching_size;
    m_time_limit = option.time_limit;
    m_root_select_rule = option.root_select_rule;
    m_print_results_on = option.print_on;

    m_contour_option = option.contour_option;
    m_tie_brak_option = option.tie_break_option;
    m_uncov_cont_option = option.uncov_cont_option;
    m_uncov_cont_param = option.uncov_cont_param;
    m_measure_best_mode = option.measure_best_mode;

    if (m_branching_strategy != 4)
    {
        m_branching_strategy = 4;
        m_contour_option = 3;
    }

    m_size_inst = m_data->getSizeInst();
    m_best_known = m_data->getUbFromFile();
    m_best_ub = m_best_known;
    m_best_sol = DBL_MAX;
    m_bnbPtr = new BranchNBound(m_data);
}

bool CETSP_solver::init_root_node()
{
    m_init_total_bnb_time = monotonicClock();
    if( m_root_select_rule == 1 ) m_root_node->pts = m_bnbPtr->selectRoot();
    if( m_root_select_rule == 2 ) m_root_node->pts = m_bnbPtr->selectRoot2();
    if( m_root_select_rule == 3 ) m_root_node->pts = m_bnbPtr->selectRoot3();
    
    if (m_print_results_on)
    {
        cout << "Initial root: ";
        for ( int i = 0; i <m_root_node->pts.size() ; i++ ){
            cout << m_root_node->pts[ i ] << " ";
        }
        cout << endl;
    }
    
    SolveSocpCplex *solveSocpPtr = new SolveSocpCplex(m_data, m_root_node->pts);
    m_init_socp_time = monotonicClock();
    solveSocpPtr->solveSOCP( m_root_node->pts );
    m_total_socp_time += ( monotonicClock() - m_init_socp_time );
    m_count_SOCP_solved++;

    m_root_node->lb = solveSocpPtr->getF_value();
    m_root_lb = m_root_node->lb;
    m_root_node->s_lb = 1;
    m_root_node->id = 0; m_root_node->depth = 0;
    solveSocpPtr->printF_value();
    solveSocpPtr->finishSOCP();
    solveSocpPtr->printSolution(m_root_node->pts);
    m_sum_violation += solveSocpPtr->violation;

    bool feasibilityTest = false;
    vector<vector<double>> solutionXYZ;
    vector<double> tempX;
    vector<double> tempY;
    vector<double> tempZ;
    for ( int i = 0; i < m_root_node->pts.size(); i++ )
    {
        tempX.push_back(solveSocpPtr-> getSolutionX(i));
        tempY.push_back(solveSocpPtr-> getSolutionY(i));
        tempZ.push_back(solveSocpPtr-> getSolutionZ(i));
    }
    feasibilityTest = m_bnbPtr->check_feasibility_Q( m_root_node, tempX, tempY, tempZ );

    if (m_print_results_on)
    {
        cout << endl;
        cout << "Not covered in root: ";
        for ( int i = 0; i < m_bnbPtr->notCoveredBalls.size(); i++ ){
            cout << m_bnbPtr->notCoveredBalls[ i ] << " ";
        }
        cout << endl;
    }

    //check not covered clients
    for ( int j = 0; j < m_bnbPtr->notCoveredBalls.size(); j++ ){
        m_root_node->notCovered.push_back( m_bnbPtr->notCoveredBalls[ j ] );
    }
    m_root_node->insert_spans.push_back(make_pair(0, m_root_node->pts.size()));

    delete solveSocpPtr;
    if (feasibilityTest)
    {
        if (m_print_results_on) cout << "FEASIBLE ROOT" << endl;
        m_best_lb = m_root_lb;
        m_best_sol = m_root_lb;

        m_solution_sequence = m_root_node->pts;
        m_solution_coords.clear();
        m_solution_coords.push_back(tempX); m_solution_coords.push_back(tempY); m_solution_coords.push_back(tempZ);
        
        m_computation_time = monotonicClock() - m_init_total_bnb_time;
        m_gap_root = ((m_best_ub - m_root_lb)/ m_best_ub)*100;

        delete m_root_node;
    }
    else
    {
        if (m_print_results_on) cout << "INFEASIBLE ROOT" << endl;
    }
    if (m_print_results_on) cout << endl;
    
    m_delete_root_node = true;
    return feasibilityTest;
}

int CETSP_solver::solve()
{
    if (m_bnbPtr == nullptr)
    {
        m_bnbPtr = new BranchNBound(m_data);
    }
    if (m_root_node == nullptr)
    {
        bool is_feasible = init_root_node();
        if (is_feasible) return 0;
    }
    else if (m_root_node->notCovered.empty())
    {
        m_best_ub = m_root_node->lb;
        return 0;
    }
    
    m_num_solves++;
    // m_size_inst = m_root_node->notCovered.size();
    
    int precision = 4;
    int printHeader = 0;
    int optimalFound = 0;
    // int print_counter = floor( m_size_inst/2 );
    // list <node*> open;
    list <node*>::iterator itOpen;

    long int itCount = m_root_node->id;
    int insert_id = 0;
    vector< double > tempX;
    vector< double > tempY;
    vector< double > tempZ;

    m_root_seq = convert_seq_to_string(m_root_node->pts);
    m_sequence_checked[m_root_seq] = false;
    m_sequence_lb[m_root_seq] = m_root_node->lb;
//     SolveSocpCplex* solveSocpPtr = new SolveSocpCplex(m_data, m_root_node->pts);
//     solveSocpPtr->solveSOCP( m_root_node->pts );
//     solveSocpPtr->finishSOCP();
//     for ( int i = 0; i < m_root_node->pts.size(); i++ )
//     {
//       tempX.push_back( solveSocpPtr-> getSolutionX( i ) );
//       tempY.push_back( solveSocpPtr-> getSolutionY( i ) );
//       tempZ.push_back( solveSocpPtr-> getSolutionZ( i ) );
//    }
//    bool feasibilityTest = m_bnbPtr->check_feasibility_Q( m_root_node, tempX, tempY, tempZ );
//    delete solveSocpPtr;
//    if (feasibilityTest)
//    {
//        m_best_ub = m_root_node->lb;
//        return 0;
//    }
//    m_root_node->notCovered.clear();
//    for ( int j = 0; j < m_bnbPtr->notCoveredBalls.size(); j++ )
//    {
//       m_root_node->notCovered.push_back(m_bnbPtr->notCoveredBalls[j]);
//    }

    // open.push_back(m_root_node);
    m_bnbPtr->add_lb(m_root_node->lb);
    m_best_lb = m_root_node->lb;

    CBFS* cbfs = new CBFS(m_data, m_contour_option, m_tie_brak_option);
    if (m_contour_option == 2)
    {
        switch (m_uncov_cont_option)
        {
        case 1:
            cbfs->setContourBinSize(m_uncov_cont_param);
            break;
        case 2:
            cbfs->setContourNumCont(m_uncov_cont_param);
            break;
        }
    }
    cbfs->setMeasureBestMode(m_measure_best_mode);
    cbfs->addNode(m_root_node);
    cbfs->resetCurContour();
    
    // int strongBranchingSize = 0;
    // if ( m_branching_rule == 1 )
    //     strongBranchingSize = 1;
    // else if ( m_branching_rule == 2 )
    //     strongBranchingSize = m_strong_branching_size;
    int strongBranchingSize = 1;
    
    int sum = 0;
    vector<sbAuxStruct> vectorOfChildren;
    vector<vector<string>> all_children_seqs;
    vectorOfChildren.resize(1);
    all_children_seqs.resize(1);
    vector<int>::iterator stBrchit;
    m_init_total_bnb_time = monotonicClock();
    // begin iteration
    // while (!open.empty() && monotonicClock() - m_init_total_bnb_time <= m_time_limit)
    while (cbfs->m_num_unexplrd_nodes != 0 && monotonicClock() - m_init_total_bnb_time <= m_time_limit)
    {
        node* current;
        if (m_branching_strategy == 4)
        {
            current = cbfs->getNextNode();
            // open.remove(current);
            m_bnbPtr->remove_lb(current->lb);
        }

        if (!current->notCovered.empty() && current->lb < m_best_ub)
        {
            stBrchit = current->notCovered.begin();
            //control strong branching size
            // if (m_branching_rule == 2)
            // {
            //     strongBranchingSize = m_strong_branching_size;
            //     strongBranchingSize -= current->pts.size() - 3;
            //     strongBranchingSize = max(1, strongBranchingSize);
            // }
            double initialTimeSB = monotonicClock();
            // if (strongBranchingSize > 1 && m_print_results_on)
            // {
            //     cout << "Doing Strong Branching: " << strongBranchingSize << endl;
            // }
            string parent_seq = convert_seq_to_string(current->pts);
            for (int t = 0; t < strongBranchingSize && t < current->notCovered.size(); t++)
            {
                insert_id = (*stBrchit);
                stBrchit++;
                sbAuxStruct tempSbStr;
                tempSbStr.index = insert_id;
                tempSbStr.sum = 0;
                SolveSocpCplex *solveSocpPtr2 = new SolveSocpCplex(m_data, current->pts.size()+1);
                solveSocpPtr2->initialize_model();
                vector<string> child_seqs;

                for (int span_id = 0; span_id < current->insert_spans.size(); span_id++)
                {
                    int start_pos = current->insert_spans[span_id].first;
                    int end_pos = current->insert_spans[span_id].second;

                    int prev_insert_pos = start_pos + 1;
                    int curr_insert_pos = start_pos + 1;
                    for (int pos = start_pos; pos < end_pos; pos++)
                    {
                        node * child = new node;
                        curr_insert_pos = pos + 1;
                        itCount++;
                        child->add_uncovered_lb = current->add_uncovered_lb;
                        child->s_lb = 0;
                        child->depth = current->depth + 1;
                        child->id = itCount;
                        child->pts = current->pts;
                        child->pts.insert(child->pts.begin() + (pos+1), insert_id);

                        m_init_socp_time = monotonicClock();
                        if (solveSocpPtr2->m_num_solves == 0)
                        {
                            solveSocpPtr2->populate_removable_constraints(child->pts);
                        }
                        else
                        {
                            solveSocpPtr2->clear_removable_constraints(prev_insert_pos, curr_insert_pos);
                            solveSocpPtr2->populate_removable_constraints(child->pts, prev_insert_pos, curr_insert_pos);
                        }
                        solveSocpPtr2->solveSOCP();
                        // solveSocpPtr2->clear_removable_constraints();
                        prev_insert_pos = curr_insert_pos;
                        m_total_socp_time += ( monotonicClock() - m_init_socp_time );
                        m_count_SOCP_solved++;
                        child->lb = solveSocpPtr2->getF_value();
                        // child->lb = max(child->lb, current->lb);
                        
                        string child_seq = convert_seq_to_string(child->pts);
                        child_seqs.push_back(child_seq);
                        m_sequence_checked[child_seq] = false;
                        // if (m_sequence_lb.count(child_seq) == 0 || m_sequence_lb[child_seq] < child->lb)
                            m_sequence_lb[child_seq] = child->lb;
                        // else
                        //     child->lb = m_sequence_lb[child_seq];

                        if (child->lb >= m_best_ub)
                        {
                            tempSbStr.sum += m_best_ub;
                            delete child;
                            continue;
                        }
                        else
                            tempSbStr.sum += child->lb;

                        m_sum_violation += solveSocpPtr2->violation;
                        int temp_size_curr_seq = child->pts.size();
                        tempX.resize(temp_size_curr_seq);
                        tempY.resize(temp_size_curr_seq);
                        tempZ.resize(temp_size_curr_seq);
                        for ( int i2 = 0; i2 < child->pts.size(); i2++ ){
                            tempX[i2] = solveSocpPtr2-> getSolutionX(i2);
                            tempY[i2] = solveSocpPtr2-> getSolutionY(i2);
                            tempZ[i2] = solveSocpPtr2-> getSolutionZ(i2);
                        }
                        //check feasibility
                        bool feasibilityTest;
                        feasibilityTest = m_bnbPtr->check_feasibility_Q( current, child, pos, insert_id, tempX, tempY, tempZ );
                        child->feasible = (feasibilityTest) ? 1 : 0;
                        if (child->feasible == 1)
                        {
                            if (child->lb < m_best_sol)
                            {
                                m_best_sol = child->lb;
                                m_solution_sequence = child->pts;
                                m_solution_coords.clear();
                                m_solution_coords.push_back(tempX); m_solution_coords.push_back(tempY); m_solution_coords.push_back(tempZ);
                                m_iter_to_incumbent = itCount;
                                if (m_best_sol < m_best_known)
                                {
                                    m_best_ub = m_best_sol;
                                    ++m_num_fea_sol_found;                           
                                }
                            }
                            delete child;
                        }
                        else
                        {
                            child->notCovered = m_bnbPtr->notCoveredBalls;
                            //============================================================
                            child->insert_spans.push_back(make_pair(start_pos, end_pos+1));
                            //============================================================
                            tempSbStr.candidates.push_back( child );
                            child->uncovered_node_min_dist_to_edges.clear();
                        }
                    }
                }

                m_bnbPtr->m_uncov_nodes_dists_to_edges.clear();
                all_children_seqs[0] = child_seqs;
                vectorOfChildren[0] = tempSbStr;
                delete solveSocpPtr2;
            }
            m_sb_computation_time +=  monotonicClock() - initialTimeSB;

            int pos = 0;

            m_sequence_children[parent_seq] = all_children_seqs[pos];
            // do branching for the selected vertex
            for ( int s = 0; s < vectorOfChildren[ pos ].candidates.size(); s++ )
            {
                if( vectorOfChildren[ pos ].candidates[ s ]->lb < m_best_ub )
                {
                    // open.push_back( vectorOfChildren[ pos ].candidates[ s ] );
                    m_bnbPtr->add_lb(vectorOfChildren[pos].candidates[s]->lb);
                    cbfs->addNode(vectorOfChildren[pos].candidates[s]);
                }
                else
                {
                    delete vectorOfChildren[ pos ].candidates[ s ];
                }
            }
            // update global lower bound
            m_best_lb = m_bnbPtr->get_current_best_lb();
            // vectorOfChildren.clear();
            // all_children_seqs.clear();
        }
        m_iter_count++;
        cbfs->delNode(current);
        delete current;
    }
    m_computation_time = monotonicClock() - m_init_total_bnb_time;
    if (m_computation_time > m_time_limit) optimalFound = 1;
    
    // update all lbs based on the recoreded lb info
    update_all_lb();
    // clean up
    m_bnbPtr->clear_lb();
    cbfs->clean_up();
    delete cbfs;
    
    return optimalFound;
}

int CETSP_solver::solve_keep_node()
{
    if (m_bnbPtr == nullptr)
    {
        m_bnbPtr = new BranchNBound(m_data);
    }
    if (m_root_node == nullptr)
    {
        bool is_feasible = init_root_node();
        if (is_feasible) return 0;
    }
    else if (m_root_node->notCovered.empty())
    {
        m_best_ub = m_root_node->lb;
        return 0;
    }
    
    m_num_solves++;    
    int precision = 4;
    int printHeader = 0;
    int optimalFound = 0;
    list <node*>::iterator itOpen;

    long int itCount = m_root_node->id;
    int insert_id = 0;
    vector< double > tempX;
    vector< double > tempY;
    vector< double > tempZ;
    // m_root_seq = convert_seq_to_string(m_root_node->pts);
    // m_sequence_checked[m_root_seq] = false;
    // m_sequence_nodes[m_root_seq] = m_root_node;

    m_bnbPtr->add_lb(m_root_node->lb);
    m_best_lb = m_root_node->lb;

    CBFS* cbfs = new CBFS(m_data, m_contour_option, m_tie_brak_option);
    if (m_contour_option == 2)
    {
        switch (m_uncov_cont_option)
        {
        case 1:
            cbfs->setContourBinSize(m_uncov_cont_param);
            break;
        case 2:
            cbfs->setContourNumCont(m_uncov_cont_param);
            break;
        }
    }
    cbfs->setMeasureBestMode(m_measure_best_mode);
    cbfs->addNode(m_root_node);
    cbfs->resetCurContour();
    
    int strongBranchingSize = 1;
    
    int sum = 0;
    vector<sbAuxStruct> vectorOfChildren;
    // vector<vector<string>> all_children_seqs;
    vectorOfChildren.resize(1);
    // all_children_seqs.resize(1);
    vector<int>::iterator stBrchit;
    // vector<node*> temp_store_all_nodes;
    m_init_total_bnb_time = monotonicClock();
    // begin iteration
    while (cbfs->m_num_unexplrd_nodes != 0 && monotonicClock() - m_init_total_bnb_time <= m_time_limit)
    {
        node* current;
        if (m_branching_strategy == 4)
        {
            current = cbfs->getNextNode();
            // open.remove(current);
            m_bnbPtr->remove_lb(current->lb);
        }

        if (!current->notCovered.empty() && current->lb < m_best_ub)
        {
            stBrchit = current->notCovered.begin();
            double initialTimeSB = monotonicClock();
            // string parent_seq = convert_seq_to_string(current->pts);
            for (int t = 0; t < strongBranchingSize && t < current->notCovered.size(); t++)
            {
                insert_id = (*stBrchit);
                stBrchit++;
                sbAuxStruct tempSbStr;
                tempSbStr.index = insert_id;
                tempSbStr.sum = 0;
                SolveSocpCplex *solveSocpPtr2 = new SolveSocpCplex(m_data, current->pts.size()+1);
                solveSocpPtr2->initialize_model();
                vector<string> child_seqs;

                for (int span_id = 0; span_id < current->insert_spans.size(); span_id++)
                {
                    int start_pos = current->insert_spans[span_id].first;
                    int end_pos = current->insert_spans[span_id].second;

                    int prev_insert_pos = start_pos + 1;
                    int curr_insert_pos = start_pos + 1;
                    for (int pos = start_pos; pos < end_pos; pos++)
                    {
                        node * child = new node;
                        curr_insert_pos = pos + 1;
                        itCount++;
                        child->add_uncovered_lb = current->add_uncovered_lb;
                        child->s_lb = 0;
                        child->depth = current->depth + 1;
                        child->id = itCount;
                        child->pts = current->pts;
                        child->pts.insert(child->pts.begin() + (pos+1), insert_id);

                        m_init_socp_time = monotonicClock();
                        if (solveSocpPtr2->m_num_solves == 0)
                        {
                            solveSocpPtr2->populate_removable_constraints(child->pts);
                        }
                        else
                        {
                            solveSocpPtr2->clear_removable_constraints(prev_insert_pos, curr_insert_pos);
                            solveSocpPtr2->populate_removable_constraints(child->pts, prev_insert_pos, curr_insert_pos);
                        }
                        solveSocpPtr2->solveSOCP();
                        prev_insert_pos = curr_insert_pos;
                        m_total_socp_time += ( monotonicClock() - m_init_socp_time );
                        m_count_SOCP_solved++;
                        child->lb = solveSocpPtr2->getF_value();
                        // LAH records
                        // string child_seq = convert_seq_to_string(child->pts);
                        // child_seqs.push_back(child_seq);
                        // m_sequence_checked[child_seq] = false;

                        if (child->lb >= m_best_ub - dbl_compare_constant)
                        {
                            tempSbStr.sum += m_best_ub;
                            // LAH records
                            // temp_store_all_nodes.push_back(nullptr);

                            delete child;
                            continue;
                        }
                        else
                        {
                            tempSbStr.sum += child->lb;
                        }

                        m_sum_violation += solveSocpPtr2->violation;
                        int temp_size_curr_seq = child->pts.size();
                        tempX.resize(temp_size_curr_seq);
                        tempY.resize(temp_size_curr_seq);
                        tempZ.resize(temp_size_curr_seq);
                        for ( int i2 = 0; i2 < child->pts.size(); i2++ ){
                            tempX[i2] = solveSocpPtr2-> getSolutionX(i2);
                            tempY[i2] = solveSocpPtr2-> getSolutionY(i2);
                            tempZ[i2] = solveSocpPtr2-> getSolutionZ(i2);
                        }
                        //check feasibility
                        bool feasibilityTest;
                        feasibilityTest = m_bnbPtr->check_feasibility_Q( current, child, pos, insert_id, tempX, tempY, tempZ );
                        child->feasible = (feasibilityTest) ? 1 : 0;
                        // LAH records
                        // temp_store_all_nodes.push_back(child);

                        if (child->feasible == 1)
                        {
                            if (m_use_fsi)
                            {
                                double alter_lb = m_bnbPtr->find_new_solution_from_curr_sequence(tempX, tempY, tempZ, child);
                                if (alter_lb > 0)
                                    child->lb = alter_lb;
                            }
                            if (child->lb < m_best_sol - dbl_compare_constant)
                            {
                                m_best_sol = child->lb;
                                m_solution_sequence = child->pts;
                                m_solution_coords.clear();
                                m_solution_coords.push_back(tempX); m_solution_coords.push_back(tempY); m_solution_coords.push_back(tempZ);
                                m_iter_to_incumbent = itCount;
                                if (m_best_sol < m_best_known - dbl_compare_constant)
                                {
                                    m_best_ub = m_best_sol;
                                    ++m_num_fea_sol_found;                              
                                }
                            }
                            delete child;
                        }
                        else
                        {
                            child->notCovered = m_bnbPtr->notCoveredBalls;
                            //============================================================
                            child->insert_spans.push_back(make_pair(start_pos, end_pos+1));
                            //============================================================
                            tempSbStr.candidates.push_back( child );
                            child->uncovered_node_min_dist_to_edges.clear();
                        }
                    }
                }

                m_bnbPtr->m_uncov_nodes_dists_to_edges.clear();
                // all_children_seqs[0] = child_seqs;
                vectorOfChildren[0] = tempSbStr;
                delete solveSocpPtr2;
            }
            m_sb_computation_time +=  monotonicClock() - initialTimeSB;

            int pos = 0;
            // m_sequence_children[parent_seq] = all_children_seqs[pos];
            // do branching for the selected vertex
            for ( int s = 0; s < vectorOfChildren[ pos ].candidates.size(); s++ )
            {
                if( vectorOfChildren[ pos ].candidates[ s ]->lb < m_best_ub )
                {
                    // open.push_back( vectorOfChildren[ pos ].candidates[ s ] );
                    m_bnbPtr->add_lb(vectorOfChildren[pos].candidates[s]->lb);
                    cbfs->addNode(vectorOfChildren[pos].candidates[s]);
                }
                else
                {
                    delete vectorOfChildren[ pos ].candidates[ s ];
                }
            }
            // update global lower bound
            m_best_lb = m_bnbPtr->get_current_best_lb();
        }
        m_iter_count++;
        cbfs->delNode(current);
        delete current;
    }
    m_computation_time = monotonicClock() - m_init_total_bnb_time;
    if (m_computation_time > m_time_limit) optimalFound = 1;
    
    if (optimalFound == 1)
    {
        cbfs->get_all_unexplored_nodes(m_unexplored_nodes);
    }
    // if (optimalFound == 0)
    // {
    //     for (auto item : temp_store_all_nodes)
    //     {
    //         if (item != nullptr) delete item;
    //     }
    // }
    // else
    // {
    //     for (auto item : temp_store_all_nodes)
    //     {
    //         if (item != nullptr)
    //         {
    //             string temp_sequence = 
    //         }
    //     }
    // }
    
    // update all lbs based on the recoreded lb info
    // update_all_lb();
    // clean up
    m_bnbPtr->clear_lb();
    cbfs->mContours.clear();
    // cbfs->clean_up();
    delete cbfs;
    
    return optimalFound;
}

int CETSP_solver::solve_keep_history()
{
    if (m_bnbPtr == nullptr)
    {
        m_bnbPtr = new BranchNBound(m_data);
    }
    if (m_root_node == nullptr)
    {
        bool is_feasible = init_root_node();
        if (is_feasible) return 0;
    }
    else if (m_root_node->notCovered.empty())
    {
        m_best_ub = m_root_node->lb;
        return 0;
    }
    
    m_num_solves++;    
    int precision = 4;
    int printHeader = 0;
    int optimalFound = 0;
    list <node*>::iterator itOpen;

    long int itCount = m_root_node->id;
    int insert_id = 0;
    vector< double > tempX;
    vector< double > tempY;
    vector< double > tempZ;
    m_root_seq = convert_seq_to_string(m_root_node->pts);
    m_sequence_checked[m_root_seq] = false;

    m_bnbPtr->add_lb(m_root_node->lb);
    m_best_lb = m_root_node->lb;

    CBFS* cbfs = new CBFS(m_data, m_contour_option, m_tie_brak_option);
    if (m_contour_option == 2)
    {
        switch (m_uncov_cont_option)
        {
        case 1:
            cbfs->setContourBinSize(m_uncov_cont_param);
            break;
        case 2:
            cbfs->setContourNumCont(m_uncov_cont_param);
            break;
        }
    }
    cbfs->setMeasureBestMode(m_measure_best_mode);
    cbfs->addNode(m_root_node);
    cbfs->resetCurContour();
    
    int strongBranchingSize = 1;
    
    int sum = 0;
    vector<sbAuxStruct> vectorOfChildren;
    vector<vector<string>> all_children_seqs;
    vectorOfChildren.resize(1);
    all_children_seqs.resize(1);
    vector<int>::iterator stBrchit;
    vector<node*> temp_store_all_nodes;
    temp_store_all_nodes.push_back(m_root_node);
    m_init_total_bnb_time = monotonicClock();
    // begin iteration
    while (cbfs->m_num_unexplrd_nodes != 0 && monotonicClock() - m_init_total_bnb_time <= m_time_limit)
    {
        node* current;
        if (m_branching_strategy == 4)
        {
            current = cbfs->getNextNode();
            m_bnbPtr->remove_lb(current->lb);
        }

        if (!current->notCovered.empty() && current->lb < m_best_ub)
        {
            stBrchit = current->notCovered.begin();
            double initialTimeSB = monotonicClock();
            string parent_seq = convert_seq_to_string(current->pts);
            for (int t = 0; t < strongBranchingSize && t < current->notCovered.size(); t++)
            {
                insert_id = (*stBrchit);
                stBrchit++;
                sbAuxStruct tempSbStr;
                tempSbStr.index = insert_id;
                tempSbStr.sum = 0;
                SolveSocpCplex *solveSocpPtr2 = new SolveSocpCplex(m_data, current->pts.size()+1);
                solveSocpPtr2->initialize_model();
                vector<string> child_seqs;

                for (int span_id = 0; span_id < current->insert_spans.size(); span_id++)
                {
                    int start_pos = current->insert_spans[span_id].first;
                    int end_pos = current->insert_spans[span_id].second;

                    int prev_insert_pos = start_pos + 1;
                    int curr_insert_pos = start_pos + 1;
                    for (int pos = start_pos; pos < end_pos; pos++)
                    {
                        node * child = new node;
                        curr_insert_pos = pos + 1;
                        itCount++;
                        child->add_uncovered_lb = current->add_uncovered_lb;
                        child->s_lb = 0;
                        child->depth = current->depth + 1;
                        child->id = itCount;
                        child->pts = current->pts;
                        child->pts.insert(child->pts.begin() + (pos+1), insert_id);

                        m_init_socp_time = monotonicClock();
                        if (solveSocpPtr2->m_num_solves == 0)
                        {
                            solveSocpPtr2->populate_removable_constraints(child->pts);
                        }
                        else
                        {
                            solveSocpPtr2->clear_removable_constraints(prev_insert_pos, curr_insert_pos);
                            solveSocpPtr2->populate_removable_constraints(child->pts, prev_insert_pos, curr_insert_pos);
                        }
                        solveSocpPtr2->solveSOCP();
                        prev_insert_pos = curr_insert_pos;
                        m_total_socp_time += ( monotonicClock() - m_init_socp_time );
                        m_count_SOCP_solved++;
                        child->lb = solveSocpPtr2->getF_value();
                        // LAH records
                        string child_seq = convert_seq_to_string(child->pts);
                        child_seqs.push_back(child_seq);
                        m_sequence_checked[child_seq] = false;
                        temp_store_all_nodes.push_back(child);

                        if (child->lb >= m_best_ub - dbl_compare_constant)
                        {
                            tempSbStr.sum += m_best_ub;
                            // delete child;
                            continue;
                        }
                        else
                        {
                            tempSbStr.sum += child->lb;
                        }

                        m_sum_violation += solveSocpPtr2->violation;
                        int temp_size_curr_seq = child->pts.size();
                        tempX.resize(temp_size_curr_seq);
                        tempY.resize(temp_size_curr_seq);
                        tempZ.resize(temp_size_curr_seq);
                        for ( int i2 = 0; i2 < child->pts.size(); i2++ ){
                            tempX[i2] = solveSocpPtr2-> getSolutionX(i2);
                            tempY[i2] = solveSocpPtr2-> getSolutionY(i2);
                            tempZ[i2] = solveSocpPtr2-> getSolutionZ(i2);
                        }
                        //check feasibility
                        bool feasibilityTest;
                        feasibilityTest = m_bnbPtr->check_feasibility_Q( current, child, pos, insert_id, tempX, tempY, tempZ );
                        child->feasible = (feasibilityTest) ? 1 : 0;

                        if (child->feasible == 1)
                        {
                            if (m_use_fsi)
                            {
                                double alter_lb = m_bnbPtr->find_new_solution_from_curr_sequence(tempX, tempY, tempZ, child);
                                if (alter_lb > 0)
                                    child->lb = alter_lb;
                            }
                            if (child->lb < m_best_sol - dbl_compare_constant)
                            {
                                m_best_sol = child->lb;
                                m_solution_sequence = child->pts;
                                m_solution_coords.clear();
                                m_solution_coords.push_back(tempX); m_solution_coords.push_back(tempY); m_solution_coords.push_back(tempZ);
                                m_iter_to_incumbent = itCount;
                                if (m_best_sol < m_best_known - dbl_compare_constant)
                                {
                                    m_best_ub = m_best_sol;
                                    ++m_num_fea_sol_found;                              
                                }
                            }
                            // delete child;
                        }
                        else
                        {
                            child->notCovered = m_bnbPtr->notCoveredBalls;
                            //============================================================
                            child->insert_spans.push_back(make_pair(start_pos, end_pos+1));
                            //============================================================
                            tempSbStr.candidates.push_back( child );
                            child->uncovered_node_min_dist_to_edges.clear();
                        }
                    }
                }

                m_bnbPtr->m_uncov_nodes_dists_to_edges.clear();
                all_children_seqs[0] = child_seqs;
                vectorOfChildren[0] = tempSbStr;
                delete solveSocpPtr2;
            }
            m_sb_computation_time +=  monotonicClock() - initialTimeSB;

            int pos = 0;
            m_sequence_children[parent_seq] = all_children_seqs[pos];
            // do branching for the selected vertex
            for ( int s = 0; s < vectorOfChildren[ pos ].candidates.size(); s++ )
            {
                if( vectorOfChildren[ pos ].candidates[ s ]->lb < m_best_ub )
                {
                    m_bnbPtr->add_lb(vectorOfChildren[pos].candidates[s]->lb);
                    cbfs->addNode(vectorOfChildren[pos].candidates[s]);
                }
                // else
                // {
                //     delete vectorOfChildren[ pos ].candidates[ s ];
                // }
            }
            // update global lower bound
            m_best_lb = m_bnbPtr->get_current_best_lb();
        }
        m_iter_count++;
        cbfs->delNode(current);
        // delete current;
    }
    m_computation_time = monotonicClock() - m_init_total_bnb_time;
    if (m_computation_time > m_time_limit) optimalFound = 1;
    
    if (optimalFound == 0)
    {
        for (auto item : temp_store_all_nodes)
        {
            if (item != nullptr) delete item;
        }
    }
    else
    {
        for (auto item : temp_store_all_nodes)
        {
            if (item != nullptr)
            {
                string temp_sequence_str = convert_seq_to_string(item->pts);
                // m_sequence_lb[temp_sequence_str] = item->lb;
                m_sequence_nodes[temp_sequence_str] = lah_node_record(item->lb,
                    item->depth, item->feasible, item->notCovered, 
                    item->not_in_sequence_covered_node_locations);
                delete item;
            }
        }
        // update all lbs based on the recoreded lb info
        // update_all_lb();
        update_all_lb_by_node_info();
    }
    
    // clean up
    m_bnbPtr->clear_lb();
    // cbfs->mContours.clear();
    delete cbfs;
    
    return optimalFound;
}

CETSP_solver::~CETSP_solver()
{
    if (m_bnbPtr != nullptr)
        delete m_bnbPtr;
}

void CETSP_solver::set_intersect_tol(double tol)
{
    if (m_bnbPtr != nullptr)
        m_bnbPtr->m_intersect_tol = tol;
}

void CETSP_solver::set_root_node(node* root)
{
    copy_node(root);
    m_root_lb = m_root_node->lb;
}

void CETSP_solver::copy_node(node* org_node)
{
    m_root_node = new node;
    m_root_node->pts = org_node->pts;
    // m_root_node->solXYZ = org_node->solXYZ;
    m_root_node->notCovered = org_node->notCovered;
    m_root_node->lb = org_node->lb;
    m_root_node->s_lb = org_node->s_lb;
    m_root_node->feasible = org_node->feasible;
    m_root_node->contour = org_node->contour;
    m_root_node->depth = org_node->depth;
    m_root_node->id = org_node->id;
    m_root_node->uncov_est = org_node->uncov_est;
    m_root_node->add_uncovered_lb = org_node->add_uncovered_lb;
    m_root_node->insert_spans = org_node->insert_spans;
    m_root_node->uncovered_node_min_dist_to_edges = org_node->uncovered_node_min_dist_to_edges;
    m_root_node->not_in_sequence_covered_node_locations = org_node->not_in_sequence_covered_node_locations;
}

void CETSP_solver::clean_up()
{
    m_sequence_checked.clear();
    m_sequence_children.clear();
    m_best_ub = m_best_known;
    m_best_sol = DBL_MAX;
    m_count_SOCP_solved = 0;
    // Time related parameters
    m_init_total_bnb_time = 0;
    m_computation_time = 0;
    m_init_socp_time = 0;
    m_sb_computation_time = 0;
    m_total_socp_time = 0;
    // Iteration related paramters
    m_iter_count = 0;
    m_iter_to_incumbent = 0;
    m_num_lnodes = 0;
    m_num_nodes = 0;
    m_sum_violation = 0;
    m_num_fea_sol_found = 0;
}

void CETSP_solver::update_all_lb()
{
    update_all_lb_helper(m_root_seq);
}

double CETSP_solver::update_all_lb_helper(string current_seq)
{
    double min_temp_lb = DBL_MAX;
    if (m_sequence_children.count(current_seq) == 0 || m_sequence_children[current_seq].empty())
    {
        min_temp_lb = m_sequence_lb[current_seq];
        m_sequence_checked[current_seq] = true;
        return min_temp_lb;
    }
    for (auto child_seq : m_sequence_children[current_seq])
    {
        double temp_lb = DBL_MAX;
        if (m_sequence_children.count(child_seq) == 0 || m_sequence_children[child_seq].empty())
        {
            temp_lb = m_sequence_lb[child_seq];
            m_sequence_checked[child_seq] = true;
        }
        else if (m_sequence_checked.count(child_seq) == 0 || !m_sequence_checked[child_seq])
        {
            temp_lb = update_all_lb_helper(child_seq);
        }
        else
        {
            // Should nenver go into this condition if working correctly
            temp_lb = m_sequence_lb[child_seq];
            m_sequence_checked[child_seq] = true;
        }
        min_temp_lb = min(min_temp_lb, temp_lb);
    }

    // If for some reason min_temp_lb is not updated, just set lb to 0
    // min_temp_lb = (min_temp_lb == DBL_MAX) ? 0 : min_temp_lb;

    m_sequence_checked[current_seq] = true;
    m_sequence_lb[current_seq] = min_temp_lb;
    // if (m_sequence_nodes.count(current_seq) > 0)
    // {
    //     m_sequence_nodes[current_seq].lb = min_temp_lb;
    // }
    return min_temp_lb;
}

void CETSP_solver::update_all_lb_by_node_info()
{
    update_all_lb_by_node_info_helper(m_root_seq);
}

double CETSP_solver::update_all_lb_by_node_info_helper(string current_seq)
{
    double min_temp_lb = DBL_MAX;
    if (m_sequence_children.count(current_seq) == 0 || m_sequence_children[current_seq].empty())
    {
        min_temp_lb = m_sequence_nodes[current_seq].lb;
        m_sequence_checked[current_seq] = true;
        return min_temp_lb;
    }
    for (auto child_seq : m_sequence_children[current_seq])
    {
        double temp_lb = DBL_MAX;
        if (m_sequence_children.count(child_seq) == 0 || m_sequence_children[child_seq].empty())
        {
            temp_lb = m_sequence_nodes[child_seq].lb;
            m_sequence_checked[child_seq] = true;
        }
        else if (m_sequence_checked.count(child_seq) == 0 || !m_sequence_checked[child_seq])
        {
            temp_lb = update_all_lb_by_node_info_helper(child_seq);
        }
        else
        {
            // Should nenver go into this condition if working correctly
            temp_lb = m_sequence_nodes[child_seq].lb;
            m_sequence_checked[child_seq] = true;
        }
        min_temp_lb = min(min_temp_lb, temp_lb);
    }

    m_sequence_checked[current_seq] = true;
    m_sequence_nodes[current_seq].lb = min_temp_lb;

    return min_temp_lb;
}

bool CETSP_solver::reduce_uncovered_at_root()
{
    // double ratio_notcovered = static_cast<double>(m_root_node->notCovered.size()) / 
    //                 static_cast<double>(m_size_inst);
    // int num_notcovered_to_keep = m_root_node->notCovered.size() * (1 - ratio_notcovered) / 2;
    int num_notcovered_to_keep = min(m_root_node->pts.size(), m_root_node->notCovered.size());
    if (num_notcovered_to_keep <= 2)
        return false;
    m_root_node->notCovered.resize(num_notcovered_to_keep);
    return true;
}

void CETSP_solver::reset_root_node()
{
    delete m_root_node;
}