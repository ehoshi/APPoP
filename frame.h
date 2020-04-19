#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <mpi.h>
#include <armadillo>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#include <errno.h>
#include <sys/time.h>
#include <sstream>
#include <cstdio>
#include <ctime>
#include <limits>
#include <algorithm>
#include <sys/stat.h>


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
   birthday = time that point was initiated
   as_of = time at which the value was evaluated
*/
struct point
{
   double coord[MAXDIM];
   int dim;
   double value;
   double error;
   int PointID;
   int ProcID;
   int VertID;
   int orderID;
   bool overlap;
   Point_kind Pkind;
   WORKER_status status;
   MASTER_comand ComEx;
   timeval birthday;
   timeval as_of;
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
   point vertex[MAXVERT];
   point aux[MAXAUX];
   double base[MAXDIM];

   int vertices;
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


#if 0
std::ostream& operator<<(std::ostream& os, const point& p)
{
   for (int i = 0; i < p.dim; ++i) {
      os << p.coord[i] << std::endl;
   }

   return os;
}

std::ostream& operator<<(std::ostream& os, const simplex& s)
{
   using std::endl;

   for (int i = 0; i < s.vertices; ++i) {
      os << "The coord of vertex" << i << ": \n";
      for (int j = 0; j < s.vertex[i].dim; j++) {
         os << s.vertex[i].coord[j] << "\t";
      }
      os << endl;
   }

   os << "number of vertices: " << s.vertices << endl;
   os << "vert with lowest value: " << s.lowest << endl;
   os << "vert with highest value: " << s.highest << endl;

   return os;
}
#endif


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     LIST OF FUNCTIONS USED SOMEWHERE IN THE PROGRAM
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// calculation and other functions
int gcalc(double &, double &, WORKER_status &, std::ofstream &, unsigned long int, double, double);
int gr_calc(std::string, std::string, double &, double&);
bool has_nan(const arma::mat &);
bool is_file_useful(std::ostringstream &);
void GlobalPrintDebug(std::string, const int);
void PrintDebug(const std::string &, std::ofstream &, const int);
void SetVec(arma::vec &, const int, const int, const simplex);
bool ChkErr(const simplex blob, bool );

// simplex related
void reflect(simplex *);
void extend(simplex *);
void contract(simplex *);
void collapse(simplex *);
point slide(point, point, double);
int order1(simplex *);
int order2(simplex *, bool);
void getbase(simplex *);
void report(std::string, verbosity);
void dart(point *);
void printpoint(std::string, point);
enum simplex_action optimize(simplex *, point **, bool &);
int compare(point, point);
void swapPoints(simplex &, int, int);
void resetRound(simplex &);

// fork-exec and status reltated functions (mainly used by MASTER)
int StatChk(const pid_t, WORKER_status &, const int);
void killJob(const pid_t);
void PrintWsit(const WRK_situation, std::ofstream &);
enum MASTER_comand MASwork(const WORKER_status, const int, const simplex);
void updateMWstat(const WORKER_status, const int, simplex &);
void updateVALUEs(const int, const double, const double, simplex &);
MASTER_comand GetMcomm(const int, const simplex);
void printsimplex(const simplex &, std::ostream &);
bool DetectChangeWstat(const simplex, const WORKER_status, int);
WORKER_status getPrevWstat(const simplex, int);

// Only function other than worker that calls vfork&exec
int DoScript(std::string);
