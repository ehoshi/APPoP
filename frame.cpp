#include "frame.h"

#include <csignal>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <armadillo>
#include "cxxopts/cxxopts.hpp"
#include <mpi.h>


std::ofstream GLOBAL_debuglog;


/* \brief Called by atexit() handler before normal termination
 */
static void
cleanup(void)
{
   MPI_Finalize();
}

///////////////////////Signal Handler for waitpid///////////////////////
extern "C" void
sigStuff(int signo)
{
   if (SIGCHLD == signo) {
      std::signal(SIGCHLD, sigStuff);
   }
}
////////////////////////////////////////////////////////////////////////

int
main(int argc, char *argv[])
{
   int Proc;  // total number of processes
   int Rank;  // rank of the process
   arma::mat inP;  //XXX temorary complete input file

   // MPI initializations
   // MPI can modify arguments so this step must occur before general
   // argument processing.
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &Proc);
   MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
   std::atexit(cleanup);

   // check command-line arguments
   cxxopts::Options options("frame", "Automated Parallel Parameterization of Parameters");
   options.add_options("", {
         {"h,help", "Print usage"},
         {"m,mult", "Control size of the error bars",
          cxxopts::value<double>()->default_value("0.0")},
      });
   options.custom_help("[OPTION...] inputfile");
   cxxopts::ParseResult parsed = options.parse(argc, argv);
   if (parsed.count("help")) {
      if (0 == Rank) {
         std::cout << options.help() << std::endl;
      }

      // MPI_Finalize called by atexit handler
      return 0;
   }
   if (argc != 2) {
      if (0 == Rank) {
         std::cerr << "MASTER ERROR: number of arguments" << std::endl;
         std::cerr << options.help() << std::endl;
      }

      // MPI_Finalize called by atexit handler
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

      // MPI_Finalize called by atexit handler
      return 2;
   }

   // as an option with a default value, mult should always be available
   double mult = parsed["mult"].as<double>();
   if (mult < 0) {
      if (0 == Rank) {
         std::cerr << "mult must be >= 0" << std::endl;
      }

      // MPI_Finalize called by atexit handler
      return 10;
   }

   arma::uword tempFsize = inP.n_rows/(inP.n_cols+1);
   if(0 == Rank){
      std::cout << "nProc = " << Proc << std::endl;
      std::cout << "setting mult to " << mult << std::endl;
      std::cout << "fSize = " << tempFsize << std::endl;
   }
   //temporary calculate the fSize for Proc error checking

   //n_proc required is: inP.n_rows for each vertex
   //2*tempFsize for auxR and aux2
   //+1 for the master process
   if (static_cast<arma::uword>(Proc) < inP.n_rows + tempFsize*2 + 1) {
      if (0 == Rank) {
        std::cerr << "FATAL ERROR: There are wrong number of parallel processor.\n "
                  << "The number of processor required is: n_rows + 2*fSize + 1\n" 
                  << "n_rows = number of vertex, fSize = fraction size" << std::endl;
      }

      // MPI_Finalize called by atexit handler
      return 3;
   }

   // XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
   // PLEASE PLEASE PLEASE REMEMBER that this is modifying the signal handler.
   // some part in MPI relies on SIGCHLD signal handler, so there is a possibility
   // that this handler below is breaking the MPI at some point or place.
   // To fix thsi issue, we need to come up with a plan to NOT use waitpid() on
   // the parent of the worker process.
   if (std::signal(SIGCHLD, sigStuff) == SIG_ERR) {
      std::cerr << "CRAP, I DON'T EVEN KNOW WHAT TO DO ANYMORE" << std::endl;
   }

   if (0 == Rank) {
      std::cout << "This is test print statment before master starts working" << std::endl;
      MASTER(Rank, Proc, inP, mult);
      std::cout << "Print after MASTER() function, success! in main" << std::endl;
   }
   else {
      std::cout << "This is test print statment before worker starts working from proc " << Rank << std::endl;
      WORKER(Rank);
      std::cout << "Print after WORKER() function, success! in main from proc" << Rank << std::endl;
   }

   // MPI_Finalize called by atexit handler
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
