// For compiling with Mersenne twister random number generator include MTWISTER
// compiling directive

#ifndef UTIL_H
#define UTIL_H

#include <time.h>
#include <unistd.h>
#include <sys/times.h>
#include <sys/timeb.h>
#include <sys/resource.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include "structs.h"

using namespace std;

void randomize();

void setSeed(const unsigned int seed);

/* generates an integer i in {0,...,maxValue-1} */
unsigned int intRandom(const unsigned int maxValue);

double doubleRandom(const double maxValue);

double wallClock();

double cpuTime();

string convert_seq_to_string(vector<int>& sequence);

void convert_string_to_seq(string seq_string, vector<int>& seq, int num_in_seq);

double euclidianNorm( point * p1, point * p2 );

double Norm_2( vector< double > aVector );

vector< double > make_vector( point * p1, point * p2 );

double dot_product( vector< double > vector1, vector< double > vector2 );

vector< double > cross_product( vector< double > vector1, vector< double > vector2 );

vector< double > difference( vector< double > vector1, vector< double > vector2 );

vector< double > sum_vector( vector< double > vector1, vector< double > vector2 );

vector< double > scalar_product( double lambda, vector< double > vector1 );

void generate_rev_strings(vector<int> org_seq, vector<string>& all_strings, vector<pair<int,int>>& rev_candidates, vector<pair<int,int>>& rev_candidates_endids);

void generate_rev_strings_helper(vector<int> curr_seq, int candidate_pos, int num_reversed, vector<string>& all_strings, vector<pair<int,int>>& rev_candidates, vector<pair<int,int>>& rev_candidates_endids);

long get_curr_memory_consumption();

const double dbl_compare_constant = 0.00001;

#endif /* ifndef UTIL_H */
