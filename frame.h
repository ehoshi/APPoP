#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <armadillo>


#define MAXAUX  2
#define MAXDIM  100  // number of cols in input
#define MAXVERT MAXDIM+1  // number of rows in input


// debug output file for GlobalPrintDebug
extern std::ofstream GLOBAL_debuglog;

//******WRK stat is equivalence of point_status
/*
   WORKER_status
   undef = worker waiting to start working: simulation has not begun yet.
   pending  = running simulation, but results are NOT meaningful
   active = running simulation, and results are meaningful
   done = finished evaluation
   finished = simulation has finished running (normally). This should not
              happen often, as simulation should be infintely long.
   error    = something went wrong for the processor.
   stall    = some error, but not fatal (e.g. stopped)
*/
enum WORKER_status { undef, pending, active,
                     done, finished, error, stall, infected };

/*
   MASTER_command = simplex action
   start      = command the worker to start the simulation from scratch.
   reinitate  = command the worker to restart the simulation from scratch.
                This is most likely caused by error in the simulation.
   question   = command the worker to return the status.
                the response will be different depending on the WRK response.
   nothing    = do nothing. This is here to prevent dead/live lock...and to keep 
                the code simple (somewhat).
   GMTV       = The only all-caps command. Stands for Give Me The Value.
                ...asks for workers to send the g_VALUE and g_UNC.
   contagious = used for book keeping. Tells zombified process to set its stat
                to zombie and do no further action.
   stop       = command the worker to quit working. This only happens if:
                a fatal error occurs or normal termination happens.
*/
enum MASTER_comand {start, reinitiate, question, nothing, GMTV, contagious ,stop};

/*
   actions that simplex can recommend:
   initiate (begin a set of simulations)
   swait (come back later, hopefully with more data)
   renamed because of a name conflict in sys/wait.h
   restart (end a set of simulations and begin again)
   terminate (end a set of simulations)
   collapse - this is in here, because it is somewhat of a special case
   quit (shut down fully because we're done)
   crash - unknown condition
   none - used in the beginning to prevent unecessary printing
*/
enum simplex_action { initiate, simplex_wait, restart, shrink, quit, crash, nothingyet};

/*
* I think it is self-explanatory, but this is point type
  verte = current vertex point
  auxil = current auxillary point
* This is used for Master process loop to determine what action to
  take on some of the branches.
* Basically, this is used for looping over process number, rather
  than vertex number.
* Also, the abbriviations are awkward so that there are no overlap
  for variable names.
*/
enum Point_kind {verte, auxil};

/*
 *Basically tells wether the vertex has simulation running
 or not
 */
enum Point_stat { running, not_running };
/*
   point-gives descrption of a single point
   coord(i) = ith coordinate of point
   dim = number of dimensions of space
   value = function value at point
   error = uncertainty (standard error) of fxn at that point
   pointID = id of a job-NOT processor id in parallel
   ProcID = processor id in parallel. AKA rank
   VertID = vertex number. Should not go above number of vertex.
            Does NOT count for aux ponts
   Ptype = point type. Either 0 = normal vertex, 1 = aux vertex.
   status = status of that worker (point)
   ComEx = Command to be Excuted. Place to store Mcomm.
*/
struct point
{
   arma::vec coord;
   arma::uword dim;
   double value;
   double error;
   int PointID;
   int ProcID;
   int VertID;
   int  orderID;
   bool overlap;
   Point_stat Pstat;
   Point_kind Pkind;
   WORKER_status status;
   MASTER_comand ComEx;
};

/*
   the state of various pieces of the simplex can be:
   stale (the associated contents may not be trusted)
   current (the associated contents match the current simplex coordinates)
*/
enum simplex_status { stale, current };

/*
   the second auxiliary vertex can hold points of two types:
   none (the auxiliary vertex may not be trusted)
   extension (the auxiliary vertex holds an extension)
   contraction (the auxiliary vertex holds a contraction)
*/
enum aux_type { none, extension, contraction };

/*
   data structure for a simplex:
   vertices = number of vertices in simplex (dimensionality + 1)
   vertex(k) = kth vertex of simplex
   aux(j) = jth auxiliary point associated with simplex (reflections, etc)
   base(i) = coordinates of base point in ith dimension
   lowest = vertex with smallest function value
   highest = vertex with largest function value
   sechigh = vertex with second highest function value
   age = number of generations of simplex
   sort_status = whether simplex sorting status and base are current
   reflect_status = whether reflection is current
   aux_holds = which type of provisional vertex is stored in aux[1]
*/
struct simplex
{
   int vertices;
   point vertex[MAXVERT];
   point aux[MAXAUX];
   double base[MAXDIM];
   int lowest;
   int highest;
   int sechigh;
   int age;
   enum simplex_status sort_status;
   enum simplex_status base_status;
   enum simplex_status reflect_status;
   enum aux_type aux_holds;
};

/*
   Apparently complex has special meaning in C++. To prevent confusion, the 
   name has been changed to Notsimplex instead of complex
   data structure for Notsimplex = complex:
   vertex(k) = kth vertex of complex
   auxR(j) = jth auxiliary points associated with reflection of complex
   aux2(j) = jth auxiliary points associated with secondary move of complex
             i.e. contraction and extension
   n_vert = number of vertices in complex
   n_dim = number of dimension = n_cols
   fSize = size of the fraction. Usually n_vert/n_dim
   baseLow(i) = coordinates of base for fLow
   baseHigh(i) = coordinates of base for fHigh
   baseSecHigh(i) = coordinates of base for fSecHigh
   baseCent(i) = coordinates of base for centroid f \neq fHigh
   y = objective function values of the respective fractions
   u = noise/error of the objective function values of the respective fractions
   age = number of generations of simplex
   sort_status = whether simplex sorting status and base are current
   reflect_status = whether reflection is current
   aux_holds = which type of provisional vertex is stored in aux2
*/
struct Notsimplex
{
   std::vector <point> vertex;
   std::vector <point> auxR;
   std::vector <point> aux2;

   arma::uword n_vert;
   arma::uword n_dim;
   arma::uword fSize;
   arma::vec baseLow;
   arma::vec baseHigh;
   arma::vec baseSecHigh;
   arma::vec baseCent;

   double yfCent;
   double yfLow;
   double yfSecHigh;
   double yfHigh;
   double yfR;
   double yf2;

   double ufCent;
   double ufLow;
   double ufSecHigh;
   double ufHigh;
   double ufR;
   double uf2;

   int age;
   enum simplex_status sort_status;
   enum simplex_status base_status;
   enum simplex_status reflect_status;
   enum aux_type aux_holds;
};

/*
   text verbosity levels:
   optimization (report optimization result)
   vertex (report vertex position and value)
   simplexStatus (report status of simplex)
   vertexStatus (report status of a vertex)
   progress (report simplex action)
   recommend (recommended simplex action)
   debug (debug statement)
*/
enum verbosity { optimization = 1, vertex = 2, simplexStatus = 3,
                 vertexStatus = 4, progress = 5,
                 recommend = 7, debug = 9 };

/*
   NoStart = not started
   EqStart = equilbration started/working
   ProdStart = production started
   ProdAvail = data is ready for analysis
   ProdFinish = production finished, eventhough it is not supposed to
*/
enum WRK_situation {NoStart, EqStart, ProdStart, ProdAvail, ProdFin};


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     LIST OF FUNCTIONS USED SOMEWHERE IN THE PROGRAM
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// master and worker
void MASTER(int, int, arma::mat, double);
void WORKER(int);

// calculation and other functions
int gcalc(double &, double &, WORKER_status &, std::ofstream &, unsigned long int, double, double);
int gcalc2(double &, double &, WORKER_status &, std::ofstream &, unsigned long int, double, double);
int gr_calc(std::string, std::string, double &, double&);
bool has_nan(const arma::mat &);
bool is_file_useful(std::ostringstream &);
void GlobalPrintDebug(std::string, const int);
void PrintDebug(const std::string &, std::ofstream &, const int);
void SetVec(arma::vec &, const int, const Notsimplex);
bool ChkErr(const Notsimplex blob, bool );
void calcFrac(Notsimplex &);

// simplex related
void reflect(Notsimplex *);
void extend(Notsimplex *);
void contract(Notsimplex *);
void collapse(Notsimplex *);
point slide(point, point, double);
int order1(simplex *);
int order2(Notsimplex *, bool);
void getbase(Notsimplex *);
void report(std::string, verbosity);
void dart(point *);
void printpoint(std::string, point);
enum simplex_action optimize(Notsimplex *, std::vector<int> &, bool &, double);
int compare(point, point, double);
int compare2(double, double, double, double, double);
void swapPoint(Notsimplex &, int);
void resetRound(Notsimplex &);

// fork-exec and status reltated functions (mainly used by MASTER)
int StatChk(const pid_t, WORKER_status &, const int);
void killJob(const pid_t);
void PrintWsit(const WRK_situation, std::ofstream &);
enum MASTER_comand MASwork(const WORKER_status, const int, const Notsimplex);
void updateMWstat(const WORKER_status, const int, Notsimplex &);
void updateVALUEs(const int, const double, const double, Notsimplex &);
MASTER_comand GetMcomm(const int, const Notsimplex);
void printsimplex(const Notsimplex &, std::ostream &);
bool DetectChangeWstat(const Notsimplex, const WORKER_status, int);
WORKER_status getPrevWstat(const Notsimplex, int);

//Debug/error checking functions
void CheckProcID(const Notsimplex);

// Only function other than worker that calls vfork&exec
int DoScript(std::string);


/*
 * Local Variables:
 * mode: c++
 * c-file-style: "linux"
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 *
 * vim: set autoindent expandtab shiftwidth=3 softtabstop=3 tabstop=3:
 */
