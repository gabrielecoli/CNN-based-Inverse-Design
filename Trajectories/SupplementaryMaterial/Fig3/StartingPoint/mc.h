#define NVT
#define TARGET_ACC_VOL 0.20
#define TARGET_ACC_DIS 0.45
#define TOL 0.05

#define epsilon 1.0
#define kappa 10.0
#define rc 3.5
#define Nadjust 100
#define Nbins 100

#define rc_sann 3.5
#define m_MAX 24
#define MAX_NEIGH 100

#define SymmetryIndex 5
#define Target_QC 1.00
#define SQR(x) ((x)*(x))

/************** STRUCTS ******************/
typedef struct {
  double x;
  double y;
  double z;
} vector;

typedef struct {
  double x;
  double y;
  double z;
  double V;
} box;

typedef struct {
  double x;
  double y;
  double z;
  double r;
  int id;
} neighbor;

typedef struct {
  double Re;
  double Im;
} complex;

/*************** FUNCTIONS ***************/
void init_param(box *b, double *DeltaPos, double *DeltaL, double *Dr);
void init_rand_conf(vector *par, box b);
void init_sc_conf(vector *par, box b);
void init_hex_conf(vector *par, box b);
void read_init_conf(vector *par, box *b);
void print_snapshot(vector *par, box b, unsigned long long int n, char ** argv);
void volume_move(vector *par, box *b, double DeltaL, unsigned long long int *Nattempts, unsigned long long int *Naccepted);
void displacement_move(vector *par, box b, double DeltaPos, unsigned long long int *Nattempts, unsigned long long int *Naccepted);
void rescale(vector *par, box *b, vector scaling, double newV);
double energy_particle(vector newpos, int i, vector *par, box b, vector scaling, int k);
int partial_overlap_check(vector newpos, int i, vector *par, box b, vector scaling);
int overlap_check(vector newpos, int i, vector *par, box b, vector scaling);
vector pbc(vector dist, box b, vector scaling);
double ran();
// sampling
void gr(double *g, double Dr, vector *par, box b, int *N);
void normalize_gr(double *g, double Dr, int Ng, double rho, char ** argv);
double energy_system(vector *par, box b, vector scaling);
// vectors
vector setVector(double x, double y, double z);
vector sumVectors(vector a, vector b);
vector distanceVector(vector a, vector b);
double dotProduct(vector a, vector b);
// BOPS
int compare(const void *nb1, const void *nb2);
void neighbors_sann(vector *par, box b, int *Nb, double *radius, neighbor **neighbors, double **chi, int i);
void bops(vector *par, int *Nb, double **chi, neighbor **neighbors, int i);

/*************** GLOBAL VARIABLES ***************/
int Npart, seed;
double sigma1;
double betaP=1.0, InitDensity, beta, T;
unsigned long long int Niterations, Nequilibrations;
int Nprint;
double lambda;
int InitConf;
