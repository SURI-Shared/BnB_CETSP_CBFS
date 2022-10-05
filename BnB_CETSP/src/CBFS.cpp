//======= WENDA CHANGE =======
#include"CBFS.h"

void CBFS::addNode(node* nodeData)
{
    nodeData->contour = calContour(nodeData);
    double best;
    switch (m_measure_best_mode)
    {
    case 1:
        best = nodeData->lb;
        break;
    case 2:
        best = nodeData->lb + nodeData->uncov_est;
        break;
    case 3:
        best = floor(nodeData->lb);
        break;
    default:
        best = nodeData->lb;
    }
    mContours[nodeData->contour].insert({best, nodeData});

    // if (mLBCount.count(nodeData->lb) == 0)
    //     mLBCount[nodeData->lb] = 1;
    // else
    //     mLBCount[nodeData->lb]++;
    // mNumNodes++;
    m_num_unexplrd_nodes++;
    m_max_num_unexplrd_nodes = max(m_num_unexplrd_nodes, m_max_num_unexplrd_nodes);
    //cout << "Node " << nodeData->id << " is inserted. Depth is " << nodeData->depth << endl;
}

void CBFS::delNode(node* nodeData)
{
    double best;
    switch (m_measure_best_mode)
    {
    case 1:
        best = nodeData->lb;
        break;
    case 2:
        best = nodeData->lb + nodeData->uncov_est;
        break;
    case 3:
        best = floor(nodeData->lb);
        break;
    default:
        best = nodeData->lb;
    }
    auto range = mContours[nodeData->contour].equal_range(best);
    // Remove node with exact same id as nodeData
    for (auto i = range.first; i != range.second; i++)
    {
        if (i->second->id == nodeData->id)
        {
            mContours[nodeData->contour].erase(i);
            break;
        }
    }
    //cout << "Node " << nodeData->id << " is deleted. Depth is " << nodeData->depth << endl;
    // Remove contour if it becomes empty
    if (mContours[nodeData->contour].empty())
    {
        if (mCurContour == mContours.find(nodeData->contour))
        {
            if (mCurContour == mContours.begin())
                mCurContour = mContours.end();
            --mCurContour;
            m_keep_curr_contour = false;
        }
        mContours.erase(nodeData->contour);
    }
    m_num_unexplrd_nodes--;
}

int CBFS::calContour(node* nodeData)
{
    int contour;
    int cur_num_uncov;
    switch (mCbfsMode)
    {
        // depth contour
        case 1:
            contour = nodeData->depth;
            break;

        case 2:
            contour = m_num_cont_uncov - nodeData->notCovered.size() / m_bin_size;
            break;
        // case 3:
        //     cur_num_uncov = nodeData->notCovered.size();
        //     if (cur_num_uncov > 20)
        //     {
        //         contour = (m_size_inst - cur_num_uncov) / 5;
        //     }
        //     else
        //     {
        //         contour = (m_size_inst - 20) / 5 + (21 - cur_num_uncov);
        //     }
        //     break;
        case 3:
            contour = 0;
            break;
        default:
            contour = 0;
            break;
    }
    return contour;
}
 
node* CBFS::getNextNode()
{
    node* next;
    vector<node*> candidates;
    // for now we just use cycling

    if (!m_keep_curr_contour)
    {
        mCurContour++;
    }
    else
    {
        m_keep_curr_contour = false;
    }
        
    if (mCurContour == mContours.end())
    {
        // m_prev_cycle_root_depth = m_current_cycle_root_depth;
        mCurContour = mContours.begin();
        m_from_top_level = true;
        m_num_cycle++;
    }
    else
    {
        m_from_top_level = false;
    }
        
    double search = mCurContour->second.begin()->first;
    auto range = mCurContour->second.equal_range(search);
    switch (mTiebreakMode)
    {
        // FIFO
        case 1:
            next = range.first->second;
            break;
        // LIFO
        case 2:
            next = (--range.second)->second;
            break;
        // Arbitrary
        case 3:
            for (auto iter = range.first; iter != range.second; iter++)
                candidates.push_back(iter->second);
            next = candidates[rand() % candidates.size()];
            break;
        // Should be FIFO as well?
        default:
            next = mCurContour->second.begin()->second; 
    }
    //cout << "Node " << next->id << " is retrived." << endl;
    // if (m_from_top_level)
    // {
    //     m_current_cycle_root_depth = next->depth;
    // }
    return next;
}

// int CBFS::getNumNodes(double lb)
// {
//     int count = 0;
//     for (auto const& item : mLBCount)
//     {
//         if (item.first >= lb)
//             break;
//         count += item.second; 
//     }
//     return count;
// }

void CBFS::setContourBinSize(int bin_size)
{
    m_bin_size = bin_size;
    m_num_cont_uncov = m_size_inst / m_bin_size;
}

void CBFS::setContourNumCont(int num_cont)
{
    int temp_bin_size = num_cont;
    
    if (temp_bin_size < 1)
    {
        temp_bin_size = m_size_inst;
    }
    m_num_cont_uncov = temp_bin_size;
    m_bin_size = m_size_inst / m_num_cont_uncov;
}

void CBFS::clean_up()
{
    for (auto it = mContours.begin(); it != mContours.end(); it++)
    {
        for (auto item = it->second.begin(); item != it->second.end(); item++)
        {
            delete item->second;
        }
        // (it->second).clear();
    }
    // mContours.clear();
    // mLBCount.clear();
}

void CBFS::clean_up(double curr_ub)
{
    for (auto it = mContours.begin(); it != mContours.end();)
    {
        for (auto item = it->second.lower_bound(curr_ub); item != it->second.end();)
        {
            delete item->second;
            item = it->second.erase(item);
            m_num_unexplrd_nodes--;
        }
        if (it->second.empty())
        {
            if (mCurContour == it)
            {
                if (mCurContour == mContours.begin())
                    mCurContour = mContours.end();
                --mCurContour;
                m_keep_curr_contour = false;
            }
            it = mContours.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

void CBFS::free_percentage_unexplored_nodes(double percentage)
{
    percentage = min(1.0, max(0.0, percentage));
    for (auto it = mContours.begin(); it != mContours.end();)
    {
        int num_to_remove = static_cast<int>(it->second.size() * percentage);
        while (num_to_remove > 0)
        {
            auto item = it->second.end();
            item--;
            delete item->second;
            it->second.erase(item);
            m_num_unexplrd_nodes--;
            num_to_remove--;
        }
        if (it->second.empty())
        {
            if (mCurContour == it)
            {
                if (mCurContour == mContours.begin())
                    mCurContour = mContours.end();
                --mCurContour;
                m_keep_curr_contour = false;
            }
            it = mContours.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

void CBFS::get_all_unexplored_nodes(vector<node*>& all_nodes)
{
    all_nodes.clear();
    for (auto it = mContours.begin(); it != mContours.end(); it++)
    {
        for (auto item = it->second.begin(); item != it->second.end(); item++)
        {
            all_nodes.push_back(item->second);
        }
    }
}
//======= END WENDA CHANGE =======