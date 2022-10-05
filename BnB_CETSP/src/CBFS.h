//======= WENDA CHANGE =======
#ifndef CBFS_H
#define CBFS_H

#include <vector>
#include <map>
#include <algorithm>

#include "Data.h"
#include "structs.h"

using namespace std;

typedef multimap<double, node*> CbfsHeap;
typedef map<int, CbfsHeap> ContourMap;

class CBFS
{
public:
    CBFS(Data* d, int c,int t) : m_data(d), mCbfsMode(c), mTiebreakMode(t), mNumNodes(0) 
    {
        m_size_inst = m_data->getSizeInst();
    }

    void addNode(node* nodeData);
    void delNode(node* nodeData);
    node* getNextNode();
    int calContour(node* nodeData);
    void resetCurContour() { mCurContour = mContours.begin(); }
    void setContourBinSize(int bin_size);
    void setContourNumCont(int num_cont);
    void setMeasureBestMode(int m) { m_measure_best_mode = m; }
    // int getNumNodes() { return mNumNodes; }
    // int getNumNodes(double lb);
    void clean_up();
    void clean_up(double curr_ub);
    void free_percentage_unexplored_nodes(double percentage);
    void get_all_unexplored_nodes(vector<node*>& all_nodes);

    ContourMap mContours;
    ContourMap::iterator mCurContour;
    bool m_from_top_level = false;
    int m_num_cycle = 0;
    unsigned long int m_num_unexplrd_nodes = 0;
    bool m_keep_curr_contour = false;
    unsigned long int m_max_num_unexplrd_nodes = 0;
    // int m_current_cycle_root_depth = 0;
    // int m_prev_cycle_root_depth = 0;

private:
    Data* m_data;
    map<double, int> mLBCount;

    int mNumNodes;
    int mCbfsMode;
    int mTiebreakMode;
    int m_size_inst;

    // uncover contour parameter
    int m_num_cont_uncov = 20;      // actual num of contour should be m_num_cont_uncov + 1
    int m_bin_size = 1;

    int m_measure_best_mode = 1;    // if 1: best lower bound; 2: best estimate
};

#endif
//======= END WENDA CHANGE =======