//tsp_lb.h expects a <concorde.h> but the latest Concorde release doesn't actually provide one
#ifndef CC_CONCORDE_MAIN
#define CC_CONCORDE_MAIN
#include "machdefs.h"
#include "util.h"
#include "tsp.h"

#define CC_JUST_SUBTOUR (1)
#define CC_JUST_BLOSSOM (2)
#define CC_JUST_SUBTOUR_AND_BLOSSOM (3)
#define CC_JUST_FAST_CUTS (4)

static int norm = CC_EUCLIDEAN;
static char *datfname        = (char *) NULL;
static char *edgegenfname    = (char *) NULL;
static char *problname       = (char *) NULL;
static char *probfname       = (char *) NULL;
static char *edgefname       = (char *) NULL;
static char *fullfname       = (char *) NULL;
static char *tourfname       = (char *) NULL;
static char *masterfname     = (char *) NULL;
static char *poolfname       = (char *) NULL;
#ifdef CCtsp_USE_DOMINO_CUTS
static char *dominopoolfname = (char *) NULL;
#endif
static char *restartfname    = (char *) NULL;
static char *xfname          = (char *) NULL;
static char *outfname        = (char *) NULL;
static char *filecutname     = (char *) NULL;
static int seed                 = 0;
static int nnodes_want          = 0;
static int binary_in            = 0;
static int tsplib_in            = 1;
static int gridsize             = 0;
static int just_cuts            = 0;
static int dontcutroot          = 0;
static int usetighten           = 0;
static int usedominos           = 0;
static int maxchunksize         = 16;
static int multiple_chunker     = 0;
static int valid_edges          = 0;
static int dfs_branching        = 0;
static int bfs_branching        = 1;
static int simple_branching     = 0;
static int usebranchcliques     = 1;  
static int tentative_branch_num = 0;
static int complete_price       = 0;
static int want_rcnearest       = 0;
static int output_tour_as_edges = 0;
static int run_silently         = 1;
static int be_nethost           = 0;
static int unlink_files         = 0;
static double initial_ub = CCtsp_LP_MAXDOUBLE;
static unsigned short hostport = CCtsp_HOST_PORT;
static char *grunthostname = (char *) NULL;
static char *cutbossname   = (char *) NULL;
static char *dombossname   = (char *) NULL;
static int eliminate_edges = -1;   /* Set to 1 to force elim, 0 to not elim */
static int eliminate_sparse = 0;   /* Set to 1 to elim from full edge list  */
                                   /* if the full edge list is valid        */
static int longedge_branching = 1; /* Set to 0 to turn off           */
static int save_proof = 0;         /* Set to 1 to save the proof     */
static int standalone_branch = 0;  /* Set to 1 to do a manual branch */


static void
    adjust_upbound (double *bound, int ncount, CCdatagroup *dat, int silent),
    usage (char *f);

static int
    handle_just_cuts (CCtsp_lp *lp, int the_cuts, CCrandstate *rstate,
       int silent),
    run_hk (int ncount, CCdatagroup *dat, int *hk_tour, int silent),
    build_edges (int *p_ecount, int **p_elist, int **p_elen,
        int ncount, int *ptour, CCdatagroup *dat, char *in_edgefname,
        char *in_edgegenfname, int in_just_cuts, int silent,
        CCrandstate *rstate),
    build_fulledges (int *p_excount, int **p_exlist, int **p_exlen,
        int ncount, int *ptour, char *in_fullfname),
    parseargs (int ac, char **av),
    find_tour (int ncount, CCdatagroup *dat, int *perm, double *ub,
            int trials, int silent, CCrandstate *rstate),
    getedges (CCdatagroup *dat, CCedgegengroup *plan, int ncount, int *ecount,
            int **elist, int **elen, int silent, CCrandstate *rstate),
    dump_rc (CCtsp_lp *lp, int count, char *pname, int usesparse);
int solve (int ncount, int *out_tour, double *optval, CCdatagroup *datp, int silent, char* problname, int* in_tour);
int set_out_tour(int ncount, int *perm, int *tour, int *out_tour);
#endif