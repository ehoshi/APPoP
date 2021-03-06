//#define DEBUG1
//#define DEBUG2
#include "frame.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <signal.h>
#include <sstream>
#include <string>
#include <sys/errno.h>
#include <sys/stat.h>
#include <unistd.h>

#include <armadillo>
#include <mpi.h>


/*! \brief Launch a script in the background and return its pid.
 *
 * The script must be located in the current directory and its path
 * will be formed by appending ".sh" to the end of \a name. The script
 * will be started with `argv[] = {name}`.
 *
 * \param name the name of the script to launch.
 * \return the process id of the launched script or -1 on error.
 */
static pid_t
LaunchScript(std::string name)
{
   //TODO remove when switching to spdlog, possibly by adding a custom
   // logging format that automatically adds the rank id to all
   // messages
   int rank = -1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   pid_t pid = ::vfork();
   if (!pid) {  // child process
      std::cout << "\nFRAME: This is child process on process " << rank << ".\n"
                << "About to start " << name << std::endl;

      std::string path = "./"+name+".sh";
      int chkEXC1 = ::execl(path.c_str(), name.c_str(), reinterpret_cast<char*>(0));
      if (-1 == chkEXC1) {
         std::perror("execl");
         std::cout << "FRAME ERROR: failed at execl (child) by proc " << rank << std::endl;
         ::_exit(1);
      }
      assert(false);  // should never be reached
   }

   return pid;
}

void
WORKER(int id)
{
   int dimension = 0;
   int simcount = 0;  // number of simulation "starts"
   bool done_prog = false;  // is optimization done?
   double g_VALUE = 0.0;  // obj fxn value
   double g_UNC   = 0.0;  // error associated with obj fxn
   unsigned long int temp_count = 0;

   double prev_g = 0;  // previous value of g
   double prev_u = 0;  // previous value of unc-used to min. printing
   // boolean to print particular message once per simulation
   bool Printed1 = false;

   // situations of a worker
   WRK_situation Wsit = NoStart;  // situation a worker is currently in
   // this is used for printing purposes. If the worker situation changed
   // between the loop, then it will print out what the situation is.
   WRK_situation PrevWsit = NoStart;

   // variables to keep track of MPI messages
   MPI_Status StatMPI;
   int CountMPI = 0;

   // other MPI variables
   MASTER_comand Wcomm;  // this stores the received command
   WORKER_status Wstat = undef;  // this stores the current worker status
   pid_t parent; //  parent pid
   pid_t pid = 0; //  child pid

   // declare the stringstream used for reading/making/checking the files
   std::ostringstream sysstr6;  // process directory

   // string for creating directory and moving the working directory (only applies to workers)
   sysstr6 << "Process" << id;

   // make new directory
   int chkSys6 = mkdir(sysstr6.str().c_str(), 0777);
   if (-1 == chkSys6) {
      //TODO this does not take care of a case where there is already a directory exisiting
      std::perror("mkdir");
      std::cout << "FRAME ERROR: failed at making directory by proc"<< id <<std::endl;
      Wstat = error;
      //send error status to MASTER
   }
   else {
      // change the working directory
      int chkCHDIR = ::chdir(sysstr6.str().c_str());
      if (-1 == chkCHDIR) {
         std::perror("chdir");
         std::cout << "FRAME ERROR: cannot change directory by processor" << id <<'\n'
                   << "--" << std::strerror(errno) << std::endl;
         Wstat = error;
         // send error status to MASTER
      }
      else {
         int chkcpy = std::system("cp ../runfiles/* .");
         if (-1 == chkcpy) {
            std::perror("runfiles");
            std::cout << "FRAME ERROR: cannot copy templates by processor" << id << std::endl;
            Wstat = error;
            //send error status to MASTER
         }
      }
   }

   /*
     TODO GRRRRR...!!! more error chekcing that is missing...
     MASTER does not recognize above error messages above. Therefore, edit MASTER fxn
     so that it accepts these error messages
   */

   // open the log files
   std::ostringstream outFstr;
   outFstr << "Resultid_" << id << "-" << simcount << ".log";

   GLOBAL_debuglog.open( outFstr.str().c_str() );  // open the master log file
   std::ofstream& outlog = GLOBAL_debuglog;

   // log files for keeping track of all of the g values & unc.
   std::ofstream glog;
   std::ostringstream outGstr;
   outGstr << "gLog_" << id << "_" << simcount << ".log";
   glog.open(outGstr.str().c_str());

   // coord files for generalized version. only print out coordinates, and external 
   // script will make it ready for simulation

   // recieve the dimension of the problem from MASTER
   MPI_Recv(&dimension, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &StatMPI);  // [A1]
   MPI_Get_count(&StatMPI, MPI_INT, &CountMPI);
   {
      std::ostringstream SomeString;
      SomeString << "*MPI_Status check:\n"
                 << "Recieved " << CountMPI << " count from Proc "
                 << StatMPI.MPI_SOURCE << "\n tag = " << StatMPI.MPI_TAG
                 << std::endl << "\nThe dimension fo this problem = "
                 << dimension << std::endl;
      PrintDebug(SomeString.str(), outlog, 1);
   }

   outlog << "\n========================================\n"
          << "Setup complete. Starting optimization..."
          << "\n========================================\n" << std::endl;

   while (!done_prog) {
      // Wsit is not initialized between different simulation
      // runs. This only mess up the print-time, but FIX.  at end of
      // file. i.e. when the file is closed.
      PrevWsit = Wsit;  // this is used for printing at the end of the while loop.

      PrintDebug("About to recv Mcomm", outlog, 9);
      MPI_Recv(&Wcomm, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &StatMPI);  // LINE TAG A
      MPI_Get_count(&StatMPI, MPI_INT, &CountMPI);
      {
         std::ostringstream SomeString;
         SomeString << "The Wcomm after MPI_Recv is " << Wcomm << "\n"
                    << "*MPI_Status check:\n"
                    << "Recieved " << CountMPI << " count from Proc "
                    << StatMPI.MPI_SOURCE << "\n tag = " << StatMPI.MPI_TAG;
         PrintDebug(SomeString.str(), outlog, 9);
      }

      if (start == Wcomm || reinitiate == Wcomm) {
         outlog << "initiate or reinitiate called" << std::endl;
         Wstat = pending;

         if (reinitiate == Wcomm) {
            outlog << "Reinitate called on proc " << id << std::endl;
            outlog << "Copying neccessary files..." << std::endl;

            if (EqStart == Wsit) {
               int chkSys0 = DoScript("./cpFiles_eq.sh");

               if (-1 == chkSys0) {
                  outlog << "ERROR: failed to copy dyn files by proc "<< id <<std::endl;
                  Wstat = error;
               }
            }

            if (ProdStart == Wsit || ProdAvail == Wsit) {
               int chkSys0 = DoScript("./cpFiles_prod.sh");
               if (-1 == chkSys0) {
                  outlog << "ERROR: failed to copy dyn files by proc "<< id <<std::endl;
                  Wstat = error;
               }
            }

            // This resets the Wsit, just in case
            Wsit = NoStart;

            outlog << "killing old job and reinitating new process..." <<  std::endl;
            errno = 0;  // what is this??????????????
            outlog << "pid that is terminated = " << pid << std::endl;

            killJob(pid);
            outlog << "success on terminating old simulation (?)" << std::endl;

            // close the output file from previous run
            if (outlog.is_open()) {
               const time_t Ctime = std::time(0);
               outlog << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
               outlog << "This is end of out put on process " << id
                      << "-" << simcount -1 << std::endl;
               outlog << Ctime << '\t' << std::asctime( std::localtime(&Ctime) ) << '\t' << Ctime << std::endl;
               outlog << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

               GLOBAL_debuglog.close();  // open the new log file
               glog.close();
               // open new log files
               std::ostringstream outFstr;
               outFstr << "Resultid_" << id << "-" << simcount << ".log";
               outlog.open( outFstr.str().c_str() );  // open the new log file

               std::ostringstream outGstr;
               outGstr << "gLog_" << id << "_" << simcount << ".log";
               glog.open(outGstr.str().c_str());
            }
         }  // end Wcomm == reinitate

         // this resets the g value and uncertainty from previous run
         g_VALUE = -1.0;
         g_UNC   = -1.0;

         // declare the vector for coordinate and recieve the coordinate.
         assert(dimension > 0);
         arma::vec coordinate(static_cast<arma::uword>(dimension));

         MPI_Recv(coordinate.memptr(), dimension, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &StatMPI);
         MPI_Get_count(&StatMPI, MPI_DOUBLE, &CountMPI);
         {
            std::ostringstream SString;
            SString << "*MPI_Status check:\n"
                    << "Recieved " << CountMPI << " count from Proc "
                    << StatMPI.MPI_SOURCE << "\n tag = " << StatMPI.MPI_TAG
                    << std::endl;
            PrintDebug(SString.str(), outlog, 9);
         }

         // setup sed command. putting declaration here prevents 'doubleing up' of command
         // print all coordinates to
         // name: coord.out

         std::ofstream coordOut;
         coordOut.open("coord.out");
         std::ostringstream sysstr1;
         coordOut << coordinate(0) << '\n'
                  << coordinate(1) << '\n'
                  << coordinate(2) << std::endl;
         coordOut.close();
/*
        std::ostringstream sysstr1;
        sysstr1  << "sed -e s/SIG_PLACE/" << coordinate(0) << "/ "
                 << "-e s/EPS_PLACE/" << coordinate(1) << "/ "
                 << "-e s/qH_PLACE/" << coordinate(2) << "/ "
                 << "-e s/qO_PLACE/" << (-2.0*coordinate(2)) << "/ "
                 << "temp_models.dat > models.dat";
*/
         //XXX test message:
         outlog << "+++++++++++++test message++++++++++++++++\n" << std::endl;
         outlog << "This is the inP for process" << id << '\n' << coordinate << std::endl;
         outlog << "+++++++++++++test message++++++++++++++++" << '\n' << std::endl;

//         PrintDebug("about to delete old models.dat" , outlog, 9);
         // addition of extra lines in case this statement has to run multiple times.
//         char filename[20] = "models.dat";  // XXX at the end, make this user specified name
//         int status = remove(filename);

         //TODO in here, write shell script that reads coord.out and put original sysstr1 
         //     (sed command) in this shell script
         // chkSys1-Prepare unique executable for each processor
         // sed command to make script unique for each process
         int chkSys1 = DoScript( "./PrepSIM.sh" );
         if (-1 == chkSys1) {
            outlog << "ERROR: chkSys1-failed at sed to preping PrepSIM by proc "<< id <<std::endl;
            Wstat = error;
         }
         else {
            PrintDebug("successful on PrepSIM ", outlog, 9);
            /*
              record the bg job id for later use to the log file
              nux -> this should make the first line of the
              log file to be the job id of the simulation (not the current fg job).
            */

            // fork the process. Since fork causes a mess using MPI, vfork is used
            parent = ::getpid();  // parent pid
            pid = LaunchScript("eqi");

            if (pid < 0){
               std::perror("fork");
               outlog << "ERROR: failed to fork by process " << id << std::endl;
               Wstat = error;
            }
            /*
              parent process
              The parent constantly chekcs for child's exit status every time the master questions.
              If exit status is 1, that means execl in child process failed.
              If there is any abnormal termination in the program (NOT the exec) the program 
              exits with error status of the external program
              If for some reason the external process is stopped, kill the external program, 
              and restart the external process on child process (TODO).
            */
            else {
               outlog    << "Printed from parallel process id " << id << '\n'
                         << "The pID of the child process(eqi) is "  << pid <<'\n'
                         << "The pID of the parent process (eqi) is " << parent << '\n' << std::endl;
               ::sleep(1);  // this gives time for child process to start execl

               int chkOut = StatChk(pid, Wstat, id);
               if (0 == chkOut) {  // set Wsit to equilibrium start
                  Wsit = EqStart;
               }
            }  // end parent process
        }  // end else chkSys1
         Printed1 = false;
         simcount++;
      }  // end Wcomm==start
      else if (nothing == Wcomm) {
         // Nothing, really.
      }
      else if (question == Wcomm && infected == Wstat) {
         Wstat = infected;
         // nothing, but be infected
      }
      else if (question == Wcomm && Wstat != infected) {
         {
            std::ostringstream someString;
            someString << "Processing Wsit = " << Wsit << " in Wcomm = question" << std::endl;
            PrintDebug(someString.str(), outlog, 9);
         }

         if (NoStart == Wsit) {  // Nothing has started
            // this is normal for Aux for a while
            Wstat = undef;
         }

         if (EqStart == Wsit) {
            int chkOut = StatChk(pid, Wstat, id);
            {
               std::ostringstream SomeString;
               SomeString << "TEST MESSAGE EQSTART: chkOut = " << chkOut;
               PrintDebug(SomeString.str(), outlog, 9);
               // chkOut = 0 : equilibration step is still running
               // chkOut = 1 : starting production
            }

            // begin prod
            if (1 == chkOut) {
               /*
                 At this point, Wstat will be 'finished' since that is the default output from StatChk.
                 However, Wstat should be changed in the upcoming FXNs,which will preveint this process
                 from sending the wrong status to MASTER.
               */
               Wsit = ProdStart;
               outlog << "-- starting production... " << std::endl;

               parent = ::getpid();  // parent pid
               pid = LaunchScript("prod");

               if (pid < 0) {
                  std::perror("fork");
                  outlog << "ERROR: failed to fork by process " << id << std::endl;
                  Wstat = error;
               }
               else {  // parent process 
                  std::ostringstream Pstring;
                  Pstring    << "Printed from parallel process id " << id << '\n'
                             << "The pID of the child process(prod) is "  << pid <<'\n'
                             << "The pID of the parent process (prod) is " << parent;
                  PrintDebug(Pstring.str(), outlog, 1);
                  ::sleep(1);  // this gives time for child process to start execl

                  int chkOut = StatChk(pid, Wstat, id);
                  {
                     // hopefully, Wstat is set to 'pending'.
                     const time_t Ctime = std::time(0);
                     std::ostringstream SomeString;
                     SomeString << Ctime << '\t' << std::asctime( std::localtime(&Ctime) ) << std::endl
                                << "TEST MESSAGE EQSTART-BeginProd: chkOut = " << chkOut 
                                << " and Wstat = " << Wstat;
                     PrintDebug(SomeString.str(), outlog, 5);
                  }
               }  // end parent process
            }
         }  // end Wsit == eqstart

         if (ProdStart == Wsit) {
            int chkOut = StatChk(pid, Wstat, id);
            {
               std::ostringstream SomeString;
               SomeString << "TEST MESSAGE PRODSTART: chkOut = " << chkOut << std::endl;
               PrintDebug(SomeString.str(), outlog, 9);
            }

            if (0 == chkOut) {
               int B = DoScript("./PrepObjF.sh");
               if (-1 == B) {
                  // Wstat should be pending at this point
                  PrintDebug("This MESSAGE might show up for testing", outlog, 4);
//                  Wstat = error;
               }
               else {
                  PrintDebug("DEBUG BEFOER gcalc-in PRODstart", outlog, 9);
                  int gcalc_Chk = gcalc(g_VALUE, g_UNC, Wstat, glog, temp_count, -1, -1);
                  {
                     std::ostringstream Sstring;
                     Sstring << "DEBUG AFTER gcalc-in PRODstart gcalc_Chk = " << gcalc_Chk;
                     PrintDebug(Sstring.str(), outlog, 9);
                  }
                  if (0 == gcalc_Chk) {
                     // gcalc sets Wstat = active
                     // normal execution
                     Wsit = ProdAvail;
                     std::ostringstream Sstring;
                     Sstring << "PRINTED FROM ProdStart. objFxn = "
                             << g_VALUE << " +/- " << g_UNC;
                     PrintDebug(Sstring.str(), outlog, 8);
                  }
               }  // END else system B
            }
         }  // end ProdStart

         if (ProdAvail == Wsit) {
            int chkOut = StatChk(pid, Wstat, id);
            // By default output, Wstat will be whatever StatChk spitsout.
            // But this Wstat will be modified upon sccuess
            // by gcalc before it is sent via MPI.
            {
               std::ostringstream Sstring;
               Sstring << "PRINTED FROM Wsit = ProdAvail: chkOut = " << chkOut << std::endl;
               PrintDebug(Sstring.str(), outlog, 9);
            }

            if (0 == chkOut) {
               // this means no error
               PrintDebug("PRINTED FROM INSIDE chkOut = 0", outlog, 9);

//               int B = system("./props_Linux < props.in > props.out 2>&1");
               int B = DoScript("./PrepObjF.sh");
               if (-1 == B){
                  // NOTE: XXX this does not take care of the case where props_Linux don't exist
                  outlog << "ERROR in props." << "\n Since this is happeneing after Wsit = ProdAvail, \n"
                         << "something has happened that that deleted or corrupted the file" << std::endl;
                  Wstat = error;
               }
               else {
                  prev_g = g_VALUE;  // previous value of g
                  prev_u = g_UNC;  // previous value of unc-used to min. printing

                  PrintDebug("PRINT BEFORE gcalc-in PRODavail", outlog, 9); 
                  int gcalc_Chk = gcalc(g_VALUE, g_UNC, Wstat , glog, temp_count, prev_g, prev_u);
                  {
                     std::ostringstream Sstring;
                     Sstring << "Print after gcalc:"  << " The value of gcalcChk = " 
                             << gcalc_Chk << std::endl;
                     PrintDebug(Sstring.str(), outlog, 9);
                  }
                  if (0 == gcalc_Chk) {
                     // this will set Wstat to active
                     // normal execution
                     Wstat = active;  // failsafe???
                     std::ostringstream Sstring;
                     Sstring << "PRINTED FROM ProdAvail-AFTER gcalc. objFxn = "
                             << g_VALUE << " +/- " << g_UNC;
                     PrintDebug(Sstring.str(), outlog, 8);
                  }
               }  // END else system B
            }  // end chkOut == 0
         }  // end prod avail
      }  // end Wcomm == question
      else if (GMTV == Wcomm) {
         Wstat = active;
         if (!Printed1) {
            std::ostringstream Pstring;
            Pstring << "The results are ready\n Wstat is now active ";
            const time_t Ctime = std::time(0);
            Pstring << Ctime << '\t' << std::asctime( std::localtime(&Ctime) ) << std::endl;
            PrintDebug(Pstring.str(), outlog, 5);
            Printed1 = true;
         }

         // send g_VALUE and g_UNC
         PrintDebug("sending g_VALUE back to MASTER-tag2", outlog, 9);
         MPI_Send(&g_VALUE, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
         PrintDebug("succes on sending g_VALUE", outlog, 9);

         PrintDebug("sending g_UNC back to MASTER-tag3", outlog, 9);
         MPI_Send(&g_UNC, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
         PrintDebug("succes on sending g_UNC", outlog, 9);

         if (std::abs(prev_g - g_VALUE) >= prev_g * 0.01 || std::abs(prev_u - g_UNC) >= prev_u * 0.01) {
            const time_t Ctime = std::time(0);
            std::ostringstream SomeString;
            SomeString << "The value of objFxn has changed to: " << g_VALUE << " +- "
                       << g_UNC << " " << Ctime << '\t' << std::asctime( std::localtime(&Ctime) ) ;
            PrintDebug(SomeString.str(), outlog, 5);
         }
      }  // end Wcomm == GMTV
      else if (stop == Wcomm) {
         Wstat = done;
         outlog << "This is processor " << id
                << ".\nI have received message to quit" << std::endl;
         // kill background process
         outlog << "pid that is terminated = " << pid << std::endl;
         killJob(pid);

         done_prog = true;
      }
      else if (contagious == Wcomm) {
         Wstat = infected;
      }
      else {
         std::cout << "FRAME: This is processor " << id
                   << ". Something went wrong."
                   << "\nI have received message to FORCEFULLY quit. " << std::endl;

         done_prog = true;
         ::kill(pid, SIGKILL);

         Wstat = error;
         // something went wrong... forcing exit
      }  // end else Wcomm == stop

      // Reporting time

      {
         const time_t Ctime = std::time(0);
         std::ostringstream SomeString;
         SomeString << "sending status=" << Wstat
                    << " back to MASTER ";
         SomeString << Ctime<<  '\t' << std::asctime( std::localtime(&Ctime) ) << std::endl;
         PrintDebug(SomeString.str(), outlog, 9);
      }

      MPI_Send(&Wstat, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      PrintDebug("Done sending status", outlog, 9);

      if (PrevWsit != Wsit) {
         const time_t Ctime = std::time(0);
         outlog << "=======================================\n"
                << Ctime << '\t' << std::asctime( std::localtime(&Ctime) ) << std::endl
                << "Wsit has changed. ";
         PrintWsit(Wsit, outlog);
      }

      temp_count++;
   }  // end while

   //test message to the final output file

   outlog << "Proc " << id << " successfully quit out of while loop" << std::endl;

   // close the file if file is open
   if (outlog.is_open()) {
      const time_t Ctime = std::time(0);
      outlog << "+++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      outlog << "This is the final out put on process " << id << std::endl;
      outlog << "The while loop was performed " << temp_count << " times" << std::endl;
      outlog << Ctime << '\t' << std::asctime( std::localtime(&Ctime) ) << std::endl;
      outlog << "+++++++++++++++++++++++++++++++++++++++++++" << std::endl;

      outlog.close();
   }

   if (glog.is_open()) {
      glog.close();
   }
}  // end worker


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
