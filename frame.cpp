#include "frame.h"

#include <vector>


void MASTER(int, int, arma::mat);
void WORKER(int);


std::ofstream GLOBAL_debuglog;


///////////////////////Signal Handler for waitpid///////////////////////
static void
sigStuff(int signo)
{
   if (SIGCHLD == signo) {
      signal(SIGCHLD, sigStuff);
   }
//   std::cerr << "Signal catch" << std::endl;
}
////////////////////////////////////////////////////////////////////////

int
main(int argc, char *argv[])
{
   int Proc;  // total number of processes
   int Rank;  // rank of the process
   arma::mat inP;  //XXX temorary complete input file

   // MPI initializations
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &Proc);
   MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

   // check commandline argument
   if (argc != 2) {
      if (0 == Rank) {
         std::cerr << "MASTER ERROR: numnber of argument" << std::endl;
         std::cerr << "Usage: mpirun -options ./frame inputfile" << std::endl;
      }

      MPI_Finalize();
      return 1;
   }

   /*
     All of the processors will read a same file. However, the file is small enoguh
     that it should not cause huge proglems. Although all the processors will read 
     the input, only the master will use it to prepare the full input file.
     tabunn korede iito omoimasu...
   */

   if ( !inP.load(argv[1]) ) {
      if (0 == Rank) {
         std::cerr << "FATAL ERROR: cannot load input file " << std::endl;
      }

      MPI_Finalize();
      return 2;
   }

   if (static_cast<arma::uword>(Proc) != inP.n_rows + 3) {
      if (0 == Rank) {
        std::cerr << "FATAL ERROR: There are not enough number of parallel processor.\n "
                  << "Minimum number of processor required is: ndim + 1 + 1 + 2\n" 
                  << "ndim = number of dimensions, or number of parameters" << std::endl;
      }

      MPI_Finalize();
      return 3;
   }
   // XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
   // PLEASE PLEASE PLEASE REMEMBER that this is modifying the signal handler.
   // some part in MPI relies on SIGCHLD signal handler, so there is a possibility
   // that this handler below is breaking the MPI at some point or place.
   // To fix thsi issue, we need to come up with a plan to NOT use waitpid() on
   // the parent of the worker process.
   if (signal(SIGCHLD, sigStuff) == SIG_ERR) {
      std::cerr << "CRAP, I DON'T EVEN KNOW WHAT TO DO ANYMORE" << std::endl;
   }

   if (0 == Rank) {
      std::cout << "This is test print statment before master starts working" << std::endl;
      MASTER(Rank, Proc, inP);
      std::cout << "Print after MASTER() function, success! in main" << std::endl;
   }
   else {
      std::cout << "This is test print statment before worker starts working from proc " << Rank << std::endl;
      WORKER(Rank);
      std::cout << "Print after WORKER() function, success! in main from proc" << Rank << std::endl;
   }

   MPI_Finalize();
   return 0;
}


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
