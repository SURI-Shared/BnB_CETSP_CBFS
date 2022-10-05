
#include "util.h"


#ifdef MTWISTER
/* Code from Mersenne twister random number generator */
/* http://www.math.keio.ac.jp/~matumoto/emt.html      */

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s) {
	mt[0]= s & 0xffffffffUL;
	for (mti=1; mti<N; mti++) {
		mt[mti] =
			(1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].                        */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		mt[mti] &= 0xffffffffUL;
		/* for >32 bit machines */
	}
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void) {
	unsigned long y;
	static unsigned long mag01[2]={0x0UL, MATRIX_A};
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= N) { /* generate N words at one time */
		int kk;

		if (mti == N+1)   /* if init_genrand() has not been called, */
			init_genrand(5489UL); /* a default initial seed is used */

		for (kk=0;kk<N-M;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (;kk<N-1;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}

	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void) {
	return genrand_int32()*(1.0/4294967295.0);
	/* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void) {
	return genrand_int32()*(1.0/4294967296.0);
	/* divided by 2^32 */
}

#endif

void randomize() {
    unsigned int seed;
	seed = (unsigned int)time((time_t *)NULL);

#ifdef MTWISTER
	printf("Using mersenne twister random number generator.\n");
	init_genrand(seed);
#else
	printf("Using default random number generator.\n");
	srand(seed);
#endif
	printf("Random seed is %d\n", seed);
}

void setSeed(const unsigned int seed) {
#ifdef MTWISTER
	init_genrand(seed);
#else
	srand(seed);
#endif
}

inline unsigned int intRandom(const unsigned int maxValue) {
#ifdef MTWISTER
	return ((double)genrand_real2())*((double)maxValue);
#else
	static unsigned int res;
	static unsigned int rgen;
	static double factor;
	rgen = rand();
	factor = ((double)rgen/(double)INT_MAX);

	res = (unsigned int) ( maxValue * factor );
	if (res==maxValue) res--;
	return res;
#endif
}

inline double doubleRandom(const double maxValue) {
#ifdef MTWISTER
	return maxValue * genrand_real1();
#else
	return (((double)rand())/((double)INT_MAX))*maxValue;
#endif
}

double wallClock() {
	struct timeb tp;
	double mili;

	ftime(&tp);
	mili = (double)( (tp.time)+((double)tp.millitm)/1000);

	return mili;
}

double cpuTime() {
	static struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	return ((double)usage.ru_utime.tv_sec)+(((double)usage.ru_utime.tv_usec)/((double)1000000));
}

string convert_seq_to_string(vector<int>& sequence)
{
	string res;
	for (int i = 0; i < sequence.size(); i++)
	{
		res += to_string(sequence[i]);
		if (i != sequence.size()-1)
			res += "-";
	}
	return res;
}

void convert_string_to_seq(string seq_string, vector<int>& seq, int num_in_seq)
{
	seq.resize(num_in_seq);
	int prev_pos = 0;
	int curr_pos = 0;
	for (int i = 0; i < num_in_seq - 1; i++)
	{
		curr_pos = seq_string.find('-', prev_pos);
		seq[i] = stoi(seq_string.substr(prev_pos, curr_pos-1));
		prev_pos = curr_pos+1;
	}
	seq[num_in_seq - 1] = stoi(seq_string.substr(prev_pos));
}

double euclidianNorm( point * p1, point * p2 )
{
   return sqrt( pow( p2->x - p1->x, 2) + pow( p2->y - p1->y, 2) + pow( p2->z - p1->z, 2) );
}

double Norm_2( vector< double > aVector )
{
   return sqrt( pow( aVector[ 0 ], 2) + pow( aVector[ 1 ], 2) + pow( aVector[ 2 ], 2) );
}

vector< double > make_vector( point * p1, point * p2 )
{
   vector < double > aVector;

   aVector.push_back( p2->x - p1->x );
   aVector.push_back( p2->y - p1->y );
   aVector.push_back( p2->z - p1->z );

   return aVector;
}

double dot_product( vector< double > vector1, vector< double > vector2 )
{
   return ( vector1[ 0 ]*vector2[ 0 ] + vector1[ 1 ]*vector2[ 1 ] + vector1[ 2 ]*vector2[ 2 ] );
}	

vector< double > cross_product( vector< double > vector1, vector< double > vector2 )
{
   vector< double > crossProduct;
   crossProduct.resize( 3 );

   crossProduct[ 0 ] = vector2[ 1 ]*vector1[ 2 ] - vector2[ 2 ]*vector1[ 1 ];
   crossProduct[ 1 ] = vector2[ 2 ]*vector1[ 0 ] - vector2[ 0 ]*vector1[ 2 ];
   crossProduct[ 2 ] = vector2[ 0 ]*vector1[ 1 ] - vector2[ 1 ]*vector1[ 0 ];

   return crossProduct;
}

vector< double > difference( vector< double > vector1, vector< double > vector2 )
{
   vector< double > dif;
   dif.resize( 3 );

   dif[ 0 ] = vector2[ 0 ] - vector1[ 0 ];
   dif[ 1 ] = vector2[ 1 ] - vector1[ 1 ];
   dif[ 2 ] = vector2[ 2 ] - vector1[ 2 ];

   return dif;
}

vector< double > sum_vector( vector< double > vector1, vector< double > vector2 )
{
   vector< double > sum;
   sum.resize( 3 );

   sum[ 0 ] = vector2[ 0 ] + vector1[ 0 ];
   sum[ 1 ] = vector2[ 1 ] + vector1[ 1 ];
   sum[ 2 ] = vector2[ 2 ] + vector1[ 2 ];

   return sum;
}

vector< double > scalar_product( double lambda, vector< double > vector1 )
{
   vector< double > prod;
   prod.resize( 3 );

   prod[ 0 ] = lambda*vector1[ 0 ];
   prod[ 1 ] = lambda*vector1[ 1 ];
   prod[ 2 ] = lambda*vector1[ 2 ];

   return prod;
}

void generate_rev_strings(vector<int> org_seq, vector<string>& all_strings, 
						  vector<pair<int,int>>& rev_candidates, vector<pair<int,int>>& rev_candidates_endids)
{
	generate_rev_strings_helper(org_seq, 0, 0, all_strings, rev_candidates, rev_candidates_endids);
}

void generate_rev_strings_helper(vector<int> curr_seq, int candidate_pos, int num_reversed,
								 vector<string>& all_strings, vector<pair<int,int>>& rev_candidates, 
								 vector<pair<int,int>>& rev_candidates_endids)
{
	if (candidate_pos >= rev_candidates.size()) 
	{
		if (num_reversed > 0)
			all_strings.push_back(convert_seq_to_string(curr_seq));
		return;
	}
	else
	{
		generate_rev_strings_helper(curr_seq, candidate_pos+1, num_reversed, all_strings, rev_candidates, rev_candidates_endids);

		int first_ind = rev_candidates[candidate_pos].first;
		int second_ind = rev_candidates[candidate_pos].second;
		int first_end_id = rev_candidates_endids[candidate_pos].first;
		int second_end_id = rev_candidates_endids[candidate_pos].second;
		
		if (first_end_id == curr_seq[first_ind] && second_end_id == curr_seq[second_ind])
		{
			double temp = curr_seq[first_ind + 1];
			curr_seq[first_ind + 1] = curr_seq[first_ind + 2];
			curr_seq[first_ind + 2] = temp;
			num_reversed++;
			generate_rev_strings_helper(curr_seq, candidate_pos+1, num_reversed, all_strings, rev_candidates, rev_candidates_endids);
		}
	}
}

long get_curr_memory_consumption()
{
	rusage sys_monitor_struct;
	getrusage(RUSAGE_SELF, &sys_monitor_struct);
	return sys_monitor_struct.ru_maxrss;
}