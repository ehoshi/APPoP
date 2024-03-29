//#define TEST1
//#define DEBUG1
//#define AUXDEBUG
#include "frame.h"

#include <cassert>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <thread>

#include <armadillo>
#include <mpi.h>


verbosity outputLevel;


void
MASTER(int id, int Proc, arma::mat inP, double mult)
{
   // rows and columns will be cast to int later
   assert(inP.n_rows < static_cast<arma::uword>(std::numeric_limits<int>::max()));
   assert(inP.n_cols < static_cast<arma::uword>(std::numeric_limits<int>::max()));

   bool done_prog = false;
   double g_VALUE = 0.0;
   double g_UNC = 0.0;
   MASTER_comand Mcomm = start;  // command to send to processors
   std::ofstream SIMPlog;  // only prints g_VALUE and g_UNC
   bool IsItInit;
   bool DoIprint1 = false;  // flg for when simplex takes step
   bool DoIprint2 = false;  // flg for when simplex takes step

   // variables used to keep track of MPI messages
   int CountMPI = 0;  // number of stuff recieved
   MPI_Status StatMPI;  // status of recieve, which includes tags and stuff
   WORKER_status Wtemp;  // this is used for place holder for receive Wstat

   // declare simplex-related variables
//   point *target = 0;  // target siimplex point of action
   std::vector<int> target;  // vector of index/procID for point of action
   Notsimplex blob;  // where complex info is
   enum simplex_action task = nothingyet;
   outputLevel = optimization;
//   outputLevel = debug;
   blob.sort_status = stale;
   blob.reflect_status = stale;
   blob.aux_holds = none;
   blob.age = 0;

   //XXX initialize Notsimplex

   //set n_dim
   blob.n_dim = inP.n_cols;

   // get fSize
   blob.fSize = (inP.n_rows/(inP.n_cols+1));

   //set the size of vertex vector = n_vert
   blob.vertex.resize(inP.n_rows);
 
  //Set the size of aux vectors
   blob.auxR.resize(blob.fSize);
   blob.aux2.resize(blob.fSize);
 
   //Set the size of mid-point (centroid) coord
   blob.baseLow.zeros(inP.n_cols);
   blob.baseSecHigh.zeros(inP.n_cols);
   blob.baseHigh.zeros(inP.n_cols);
   blob.baseCent.zeros(inP.n_cols);
   
   //Set the size of coord vectors
   for(arma::uword i = 0; i < inP.n_rows; i++){
      blob.vertex[i].coord.zeros(inP.n_cols);
   }

   for(arma::uword i = 0; i < blob.fSize; i++){
      blob.auxR[i].coord.zeros(inP.n_cols);
   }

   for(arma::uword i = 0; i < blob.fSize; i++){
      blob.aux2[i].coord.zeros(inP.n_cols);
   }

   blob.n_vert = inP.n_rows;
   std::cout << "FRAME: blob.n_vert =  " << blob.n_vert << std::endl;

   for (arma::uword ivert = 0; ivert < blob.n_vert; ++ivert) {
      blob.vertex.at(ivert).dim = inP.n_cols;
      blob.vertex.at(ivert).value = -1.0;
      blob.vertex.at(ivert).error = -1.0;
      blob.vertex.at(ivert).status = undef;
      blob.vertex.at(ivert).PointID = static_cast<int>(ivert);
      blob.vertex.at(ivert).ProcID  = static_cast<int>(ivert + 1);
      blob.vertex.at(ivert).VertID = static_cast<int>(ivert);
      blob.vertex.at(ivert).orderID = static_cast<int>(ivert);
      blob.vertex.at(ivert).overlap = true;
      blob.vertex.at(ivert).Pkind = verte;
      blob.vertex.at(ivert).ComEx = nothing;
   }

   for (arma::uword ivert = 0; ivert < blob.fSize; ++ivert) {
      // initialize the coord for aux. B/c it will not have
      // any coord at the beginning of the program and will
      // spit/print out garbege
      blob.auxR.at(ivert).dim = inP.n_cols;
      blob.auxR.at(ivert).value = -1.0;
      blob.auxR.at(ivert).error = -1.0;
      blob.auxR.at(ivert).status = undef;
      blob.auxR.at(ivert).PointID = static_cast<int>(ivert);
      blob.auxR.at(ivert).ProcID = static_cast<int>(blob.n_vert + ivert + 1);
      blob.auxR.at(ivert).VertID = 0;
      blob.auxR.at(ivert).orderID = -1;
      blob.auxR.at(ivert).overlap = true;
      blob.auxR.at(ivert).Pkind = auxil;
      blob.auxR.at(ivert).Pstat = not_running;
      blob.auxR.at(ivert).ComEx = nothing;
   }

   for (arma::uword ivert = 0; ivert < blob.fSize; ++ivert) {
      // initialize the coord for aux. B/c it will not have
      // any coord at the beginning of the program and will
      // spit/print out garbege
      blob.aux2.at(ivert).dim = inP.n_cols;
      blob.aux2.at(ivert).value = -1.0;
      blob.aux2.at(ivert).error = -1.0;
      blob.aux2.at(ivert).status = undef;
      blob.aux2.at(ivert).PointID = static_cast<int>(ivert);
      blob.aux2.at(ivert).ProcID = static_cast<int>(blob.n_vert + blob.fSize + ivert + 1);
      blob.aux2.at(ivert).VertID = 0;
      blob.aux2.at(ivert).orderID = -1;
      blob.aux2.at(ivert).overlap = true;
      blob.aux2.at(ivert).Pkind = auxil;
      blob.aux2.at(ivert).Pstat = not_running;
      blob.aux2.at(ivert).ComEx = nothing;
   }

   GLOBAL_debuglog.open("MASTERlog.out");  // open the master log file
   std::ofstream& MASlog = GLOBAL_debuglog;
   SIMPlog.open("SIMPlog.out");  // open the g, g_UNC log file

   MASlog << "======Starting initialization...======" << std::endl;
   IsItInit = true;

   {
      std::ostringstream Mstring;
      Mstring << "This is the inP for MASTER process:\n" << inP << std::endl;
      PrintDebug(Mstring.str(), MASlog, 2);
   }

   // Send initial start signal, size of the dimension, and
   // actual parameters to each worker processors
   // note: dimension is only sent once to all processors including aux.
   int temp_dim = static_cast<int>(inP.n_cols);
   for (int i = 0; i < Proc-1; i++) {
      std::ostringstream Mstring;
      Mstring    << "MASTER = proc" << id <<  ": dimension "
                 << " to proc "<< i+1;
      PrintDebug(Mstring.str(), MASlog, 9);
      MPI_Send(&temp_dim, 1, MPI_INT, i+1, 0, MPI_COMM_WORLD);
   }

   Mcomm = start;

   // send initialization messages
   arma::vec tempvec(inP.n_cols);
   for (arma::uword i = 0; i < blob.n_vert; i++) {
      {
         std::ostringstream Mstring;
         Mstring    << "MASTER = proc" << id <<  ": Sending MASTER_command = "
                    << Mcomm << " to proc "<< blob.vertex[i].ProcID;
         PrintDebug(Mstring.str(), MASlog, 9);
      }

      // send start message
      PrintDebug("Now sending Mcomm = start", MASlog, 9);
      MPI_Send(&Mcomm, 1, MPI_INT, blob.vertex.at(i).ProcID, 1, MPI_COMM_WORLD);
      PrintDebug("Done sending above", MASlog, 9);

      tempvec = inP.row(i).t();

      PrintDebug("Now sending coordinates", MASlog, 9);
      MPI_Send(tempvec.memptr(), temp_dim, MPI_DOUBLE, blob.vertex[i].ProcID, 2, MPI_COMM_WORLD);//[A2]
      PrintDebug("Done sending above", MASlog, 9);

      // convert armadillo row vec to an array for simplex part of the code
      for(arma::uword j = 0; j < inP.n_cols; j++) {
         blob.vertex.at(i).coord(j) = tempvec(j);
      }
   }

   // This prints out the input matrix
   {
      std::ostringstream Mstring;
      Mstring << "\n###############CHECKING STATEMENTS FOR INPUTS##################" << std::endl;
      Mstring << "Here are the print out of the blob.vertex[i].coord[j]" << std::endl;
      for (arma::uword i = 0; i < blob.n_vert; i++) {
         for (arma::uword j = 0; j < inP.n_cols; j++) {
            Mstring << blob.vertex.at(i).coord(j) << "\t";
         }
         Mstring << "\n";
      }
      PrintDebug(Mstring.str(), MASlog, 5);
   }
   PrintDebug("\n", MASlog, 5);
   // This prints out the process number, ProcID
   {
      std::ostringstream Mstring;
      Mstring << "\n###############CHECKING STATEMENTS FOR ProcID##################" << std::endl;
      Mstring << "Here are the print out of the blob.vertex[i].ProcID" << std::endl;
      for (arma::uword i = 0; i < blob.n_vert; i++) {
         Mstring << "vertex " << i << ": ProcID = " << blob.vertex.at(i).ProcID << std::endl;
      }
      for (arma::uword i = 0; i < blob.fSize; i++) {
         Mstring << "auxR "<< i << ": ProcID = " <<blob.auxR.at(i).ProcID << std::endl;
      }
      for (arma::uword i = 0; i < blob.fSize; i++) {
         Mstring << "aux2 "<< i << ": ProcID = " <<blob.aux2.at(i).ProcID << std::endl;
      }
      PrintDebug(Mstring.str(), MASlog, 5);
   }

   PrintDebug("################################################################", MASlog, 5);

   MASlog    << "\n===============initiallization status below==================\n"
             << "MASTER = proc " << id << std::endl;

   // Receive the status from initialization
   for (arma::uword i = 0; i < blob.n_vert; i++) {
      std::ostringstream Mstring1;
      std::ostringstream Mstring2;

      Mstring1 << "Receiving from proc " << i+1 << std::endl;
      PrintDebug(Mstring1.str(), MASlog, 9);

      MPI_Recv(&Wtemp, 1, MPI_INT, static_cast<int>(i+1), MPI_ANY_TAG, MPI_COMM_WORLD, &StatMPI);
      MPI_Get_count(&StatMPI, MPI_INT, &CountMPI);

      Mstring2  << "Recieved Wstat = " << Wtemp
                << " from proc "<< i+1 << "\n"
                << "which should equal source\n"
                << "Source = " << StatMPI.MPI_SOURCE << '\n'
                << "tag = " << StatMPI.MPI_TAG;
      PrintDebug(Mstring2.str(), MASlog, 9);

      blob.vertex[i].status = Wtemp;
   }

   // check and make sure that all processor did NOT return error
   // if all processors return error at initial setup, exit the program
   if (ChkErr(blob, IsItInit)) {
      Mcomm = stop;
      std::cout << "FRAME MASTER ERROR: error returned by ALL processors at setup." << std::endl;
      PrintDebug( "MASTER ERROR: error returned by ALL processors at setup.", MASlog, 1);
      for (int i = 0; i < Proc-1; i++) {
         std::ostringstream Mstring;
         Mstring    << "MASTER: Sending MASTER_command stop-" << Mcomm
                    << " to proc "<< i+1 << std::endl;
         PrintDebug(Mstring.str(), MASlog, 9);

         MPI_Send(&Mcomm, 1, MPI_INT, i+1, 3, MPI_COMM_WORLD);
         done_prog = true;
      }
   }
   else {
      MASlog << "\nSuccess on initialization. Moving on..." << std::endl;
      IsItInit = false;
   }
////////////////////////// end initialization //////////////////////////

   while (!done_prog) {
      // ask question to each processor
      // and update status in UpdateWstat
      // This only loops over Process (Rank) number
      for(int Pnum = 1; Pnum < Proc; Pnum++) {  //loop over Processor ID
         WORKER_status PREV_wstat = getPrevWstat(blob, Pnum);  // previous status: before recv. Check one by one
         std::ostringstream Mstring1;
         std::ostringstream Mstring2;
         std::ostringstream MstringR;
         std::ostringstream MstringR2;
         std::ostringstream MstringR3;

         // send 'question'
         Mstring1 << "Asking question to proc " << Pnum << "-tag4";
         PrintDebug(Mstring1.str(), MASlog, 9);

         Mcomm = question;
         MPI_Send(&Mcomm, 1, MPI_INT, Pnum, 4, MPI_COMM_WORLD);

         Mstring2 << "Finish sending question to proc " << Pnum;
         PrintDebug(Mstring2.str(), MASlog, 9);

         // recieve Wtemp/Wstat
         MstringR << "--Recieving response from proc " << Pnum;
         PrintDebug(MstringR.str(), MASlog, 9);

         MPI_Recv(&Wtemp, 1, MPI_INT, Pnum, MPI_ANY_TAG, MPI_COMM_WORLD, &StatMPI);
         MPI_Get_count(&StatMPI, MPI_INT, &CountMPI);

         MstringR2 << "Recieved Wtemp = "<< Wtemp << " from proc " << Pnum;
         PrintDebug(MstringR2.str(), MASlog, 9);

         MstringR3 << "*MPI_Status check:\n"
                   << "Recieved " << CountMPI << " count from Proc "
                   << StatMPI.MPI_SOURCE << "\n tag = " << StatMPI.MPI_TAG;
         PrintDebug(MstringR3.str(), MASlog, 9);

         updateMWstat(Wtemp, Pnum, blob);
/*
         {//XXX test message
            std::ostringstream Mstring;
            Mstring << "Current Proc number = " << Pnum << '\n'
                       << "Wtemp = " << Wtemp << std::endl;
            PrintDebug(Mstring.str(), MASlog, 1);
         }
*/
         if (DetectChangeWstat(blob, PREV_wstat, Pnum)) {
            std::ostringstream SCWstr;
            const time_t Ctime = std::time(0);
            SCWstr << ' ' << __LINE__ << " Detected change in Wstat "
                   << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            PrintDebug(SCWstr.str(), MASlog, 1);
            printsimplex(blob, MASlog);
         }
      }  // end of loop, Pnum

      // check Wlist.
      // if all processor returns error, quit with error message
      // this will be done each and everytime this loop runs
      PrintDebug("Checking for 'all process returning error'", MASlog, 9);
      if (ChkErr(blob, IsItInit)) {
         Mcomm = stop;
         std::cout << "FRAME MASTER ERROR: error returned by ALL processors." << std::endl;
         for (int i = 0; i < Proc-1; i++) {
            std::cout << "FRAME MASTER: Sending MASTER_command = " << Mcomm
                      << " to proc "<< i << std::endl;
            MPI_Send(&Mcomm, 1, MPI_INT, i+1, 16, MPI_COMM_WORLD);
         }

         done_prog = true;
         break;
      }
   
      //XXX test message
      {
         std::ostringstream Mstring;
               Mstring << "The current Mcomm set to each process are:" << std::endl;
         for(arma::uword i = 0; i < blob.n_vert; i++){
               Mstring << "vertex" << i << '\t' << "Proc" << blob.vertex.at(i).ProcID<< '\t'
                       <<  "Mcomm = " << blob.vertex.at(i).ComEx << std::endl;
         }
         for(arma::uword i = 0; i < blob.fSize; i++){
               Mstring << "auxR_" << i << '\t' << "Proc" << blob.auxR.at(i).ProcID<< '\t'
                       <<  "Mcomm = " << blob.auxR.at(i).ComEx << std::endl;
         }
         for(arma::uword i = 0; i < blob.fSize; i++){
               Mstring << "aux2_" << i << '\t' << "Proc" << blob.aux2.at(i).ProcID<< '\t'
                       <<  "Mcomm = " << blob.aux2.at(i).ComEx << std::endl;
         }
         PrintDebug(Mstring.str(), MASlog, 9);
      }


      PrintDebug("So far, not 'all error'. \n starting optimize", MASlog, 8);

////////////////////////////////////////////////////////////////////////

      std::cerr << __FILE__ << ' ' << __LINE__ << " before CheckProcID" << std::endl;
      CheckProcID(blob);

      task = optimize(&blob, target, DoIprint1, mult);

      std::cerr << __FILE__ << ' ' << __LINE__ << " before CheckProcID" << std::endl;
      CheckProcID(blob);

      switch (task) {
      case initiate:
      case restart:
      {
         /*NOTE: the check to see whether the processor is busy is done by
           the function, optimize.
           begin a set of simulations at target->coord[],
           and set the status of the point to pending
           This is the first time the aux processes are initiated
           also, send the new parameters to the starting processor
           and it SHOULD NOT touch the currently running simulations
           uses both ProcID and Pstat to determine wether the vertex needs
           to initiate or reinitiate
           loops over all processors
         */
            for(arma::uword i = 0; i < blob.n_vert; i++){
               for(arma::uword j = 0; j < target.size(); j++){
                  if(blob.vertex.at(i).ProcID == target.at(j) && blob.vertex.at(i).Pstat == not_running){
                     blob.vertex.at(i).ComEx = start;
                     std::ostringstream Mstring;
                     Mstring << "Processor" << blob.vertex.at(i).ProcID << " needs to initiate.";
                     PrintDebug(Mstring.str(), MASlog, 5);
                  }
                  if(blob.vertex.at(i).ProcID == target.at(j) && blob.vertex.at(i).Pstat == running){
                     blob.vertex.at(i).ComEx = reinitiate;
                     std::ostringstream Mstring;
                     Mstring << "Processor" << blob.vertex.at(i).ProcID << " needs to reinitiate.";
                     PrintDebug(Mstring.str(), MASlog, 5);
                  }
               }
            }
            for(arma::uword i = 0; i < blob.fSize; i++){
               for(arma::uword j = 0; j < target.size(); j++){
                  if(blob.auxR.at(i).ProcID == target.at(j) && blob.auxR.at(i).Pstat == not_running){
                     blob.auxR.at(i).ComEx = start;
                     std::ostringstream Mstring;
                     Mstring << "Processor" << blob.auxR.at(i).ProcID << " needs to initiate.";
                     PrintDebug(Mstring.str(), MASlog, 5);
                  }
                  if(blob.auxR.at(i).ProcID == target.at(j) && blob.auxR.at(i).Pstat == running){
                     blob.auxR.at(i).ComEx = reinitiate;
                     std::ostringstream Mstring;
                     Mstring << "Processor" << blob.auxR.at(i).ProcID << " needs to reinitiate.";
                     PrintDebug(Mstring.str(), MASlog, 5);
                  }
                  if(blob.aux2.at(i).ProcID == target.at(j) && blob.aux2.at(i).Pstat == not_running){
                     blob.aux2.at(i).ComEx = start;
                     std::ostringstream Mstring;
                     Mstring << "Processor" << blob.aux2.at(i).ProcID << " needs to initiate.";
                     PrintDebug(Mstring.str(), MASlog, 5);
                  }
                  if(blob.aux2.at(i).ProcID == target.at(j) && blob.aux2.at(i).Pstat == running){
                     blob.aux2.at(i).ComEx = reinitiate;
                     std::ostringstream Mstring;
                     Mstring << "Processor" << blob.aux2.at(i).ProcID << " needs to reinitiate.";
                     PrintDebug(Mstring.str(), MASlog, 5);
                  }
               }
            }
         break;
      }
      case simplex_wait:
      {
         // no need to do anything
         //     (σ･ω･)σ
/*         for(arma::uword i = 0; i < blob.n_vert; i++){
            blob.vertex.at(i).ComEx = nothing;
//            std::ostringstream Mstring;
//          Mstring << "Processor" << blob->vertex(i).ProcID << " needs to initiate.";
//          PrintDebug(Mstring.str(), MASlog, 5);
         }
         for(arma::uword i = 0; i < blob.fSize; i++){
            blob.auxR.at(i).ComEx = nothing;
            blob.aux2.at(i).ComEx = nothing;
         }
*/
         break;
      }
      case shrink:
      {
         PrintDebug("Blob needs to shrink", MASlog, 1);
         //loop over all vertices to set to 'reinitiate' initially, then go back and set
         //the command to nothing for bottom fraction
         for(arma::uword i = 0; i < blob.n_vert; i++){
            blob.vertex.at(i).ComEx = reinitiate;
         }
         // Since the vertices should be in order, vertex 0 to fSize is the lowest frac
         for (arma::uword i = 0; i < blob.fSize; i++) {
            blob.vertex.at(i).ComEx = nothing;
         }  // end for
         break;
      }  // end case collapse
      case quit:
      case crash:
      {
         done_prog = true;
         // set all ComEx to stop
         for (arma::uword i = 0; i < blob.n_vert; i++) {
            blob.vertex.at(i).ComEx= stop;
         }
         for (arma::uword i = 0; i < blob.fSize; i++) {
            blob.auxR.at(i).ComEx= stop;
         }
         for (arma::uword i = 0; i < blob.fSize; i++) {
            blob.aux2.at(i).ComEx= stop;
         }

         if (crash == task) {
            MASlog << "Unknown condition. Exiting program" << std::endl;
         }

         if (quit == task) {
            const time_t Ctime = std::time(0);
            MASlog << "\n\n===============the simplx has converged=====================\n" 
                   << "The final parameters are:\n";

            for (arma::uword i = 0; i < blob.n_vert; i++) {
               for (arma::uword j = 0; j < inP.n_cols; j++) {
                  MASlog << blob.vertex.at(i).coord(j) << "\t";
               }
               MASlog << blob.vertex.at(i).value << '\t' << blob.vertex.at(i).error;
               MASlog << std::endl;
            }

            MASlog << "The lowest value is from vertex " << blob.vertex.at(0).VertID << '\n'
                         << "The parameters are: " << std::endl;
            for (arma::uword j = 0; j < inP.n_cols; j++) {
               MASlog << blob.vertex.at(0).coord(j) << "\t";
            }
            MASlog << std::endl;
            MASlog << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
         }
         MASlog << "\n\n=====================end of program============================="
                << std::endl;
         break;
      }
      default:
      {
         break;
      }
      }  // switch(task)

      if (DoIprint1) {
         MASlog << "=====Vert swap or process starting/restarting" << std::endl;
         printsimplex(blob, MASlog);
         std::ostringstream SCWstr;
         const time_t Ctime = std::time(0);
         SCWstr << ' ' << __LINE__ << " Printed above is before send/recv "
                << "=========== "
                << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
         PrintDebug(SCWstr.str(), MASlog, 9);

         std::ostringstream Mstring;
         Mstring << "\n===========ComEx list before sending===========" << std::endl;
         for (arma::uword j = 0; j < blob.n_vert; j++) {
            Mstring << "vertex" << j << " on proc " << blob.vertex.at(j).ProcID
                    << " has ComEx = " << blob.vertex.at(j).ComEx << std::endl;
         }
         for (arma::uword k = 0; k < blob.fSize; k++) {
            Mstring << "auxR" << k << " on proc " << blob.auxR.at(k).ProcID
                    << " has ComEx = " << blob.auxR.at(k).ComEx << std::endl;
         }
         for (arma::uword k = 0; k < blob.fSize; k++) {
            Mstring << "aux2" << k << " on proc " << blob.aux2.at(k).ProcID
                    << " has ComEx = " << blob.aux2.at(k).ComEx << std::endl;
         }
         Mstring << "===================================================" << std::endl;
         PrintDebug(Mstring.str(), MASlog, 9);

//         DoIprint1 = false;
      }

      // sends and receives for optimize happens here
      for (int Pj = 1; Pj < Proc; Pj++) {

         WORKER_status PREV_wstat = getPrevWstat(blob, Pj);  // previous status: before recv. Check one by one

         Mcomm = GetMcomm(Pj, blob);

         {//XXX test message
            std::ostringstream Mstring;
            Mstring << "Current Proc number = " << Pj << '\n'
                       << "Mcomm = " << Mcomm << std::endl;
            PrintDebug(Mstring.str(), MASlog, 9);
         }

         MPI_Send(&Mcomm, 1, MPI_INT, Pj, 7, MPI_COMM_WORLD);
         PrintDebug("Done sending above", MASlog, 9);

         if(GMTV == Mcomm) {
            {
               std::ostringstream Mstring;
               Mstring << "proc" << Pj << " is active";
               PrintDebug(Mstring.str(), MASlog, 9);
            }

            {
               std::ostringstream Mstring;
               Mstring << "--recieving g_VALUE from proc " << Pj;
               PrintDebug(Mstring.str(), MASlog, 9);
            }

            MPI_Recv(&g_VALUE, 1, MPI_DOUBLE, Pj, MPI_ANY_TAG, MPI_COMM_WORLD, &StatMPI);
            MPI_Get_count(&StatMPI, MPI_DOUBLE, &CountMPI);
            {
               std::ostringstream Mstring;
               Mstring << "*MPI_Status check:\n"
                       << "Recieved " << CountMPI << " count from Proc "
                       << StatMPI.MPI_SOURCE << "\n tag = " << StatMPI.MPI_TAG;
               PrintDebug(Mstring.str(), MASlog, 9);
            }

            PrintDebug("--and now recieving g_UNC", MASlog, 9);
            MPI_Recv(&g_UNC, 1, MPI_DOUBLE, Pj, MPI_ANY_TAG, MPI_COMM_WORLD, &StatMPI);
            MPI_Get_count(&StatMPI, MPI_DOUBLE, &CountMPI);
            {
               std::ostringstream Mstring;
               Mstring << "*MPI_Status check:\n"
                       << "Recieved " << CountMPI << " count from Proc "
                       << StatMPI.MPI_SOURCE << "\n tag = " << StatMPI.MPI_TAG;
               PrintDebug(Mstring.str(), MASlog, 9);
            }

            {
               std::ostringstream Mstring;
               Mstring << "proc " << Pj << std::endl;
               Mstring << "Recv g_VALUE = " << g_VALUE << std::endl;
               Mstring << "Recv g_UNC = " << g_UNC;
               PrintDebug(Mstring.str(), MASlog, 7);
            }
            updateVALUEs(Pj, g_VALUE, g_UNC, blob);
         }

         if(reinitiate == Mcomm || start == Mcomm) {
            {
               std::ostringstream Mstring;
               Mstring << "proc" << Pj << " needs to start/reinitiate";
               PrintDebug(Mstring.str(), MASlog, 1);
            }
            // send starting parameter
            assert(temp_dim > 0);
            arma::vec tempvec(static_cast<arma::uword>(temp_dim));
            SetVec(tempvec, Pj, blob);

            {
               std::ostringstream Mstring;
               Mstring << "About to send starting parameter to proc "
                       << Pj << "-tag6";
               PrintDebug(Mstring.str(), MASlog, 9);
            }

            MPI_Send(tempvec.memptr(), temp_dim, MPI_DOUBLE, Pj, 6, MPI_COMM_WORLD);

            {
               std::ostringstream Mstring;
               Mstring << "Done sending starting parameter to proc" << Pj;
               PrintDebug(Mstring.str(), MASlog, 9);
            }
         }

         // Every time MPI communication occurs, the worker sends back the status.
         // So this is here to avoid going out of sync.

         {
            std::ostringstream Mstring;
            Mstring << "--Recieving response from proc " << Pj;
            PrintDebug(Mstring.str(), MASlog, 9);
         }
         MPI_Recv(&Wtemp, 1, MPI_INT, Pj, MPI_ANY_TAG, MPI_COMM_WORLD, &StatMPI);
         MPI_Get_count(&StatMPI, MPI_INT, &CountMPI);
         {
            std::ostringstream Mstring;
            Mstring << "*MPI_Status check:\n"
                    << "Recieved " << CountMPI << " count from Proc "
                    << StatMPI.MPI_SOURCE << "\n tag = " << StatMPI.MPI_TAG
                    << "Received wstat = " << Wtemp << '\n';
            PrintDebug(Mstring.str(), MASlog, 9);
         }

         updateMWstat(Wtemp, Pj, blob);
         
         {//XXX test message
            std::ostringstream Mstring;
            Mstring << "Current Proc number = " << Pj << '\n'
                       << "Wtemp = " << Wtemp << std::endl;
            PrintDebug(Mstring.str(), MASlog, 9);
         }

         DoIprint2 = DetectChangeWstat(blob, PREV_wstat, Pj);
         if (quit == task) {
            DoIprint2 = false;
         }
      }  // end for Pj

      // print current g and g_UNC and status
      if (DoIprint1) {
         PrintDebug("Printed below is simplex after send/recv", MASlog, 1);
         printsimplex(blob, MASlog);//print to MASlog
         std::ostringstream Mstring2;//SIMP output

         for (arma::uword j = 0; j < blob.n_vert; j++) {
            Mstring2 << blob.vertex.at(j).ProcID << ' ';
            for (arma::uword q = 0; q < inP.n_cols; q++) {
               Mstring2 << std::setw(8) << blob.vertex.at(j).coord(q) << ' ';
            }
            Mstring2 << std::setw(8) << blob.vertex.at(j).value << ' '
                     << std::setw(8) << blob.vertex.at(j).error << ' '
                     << std::endl;
         }
         PrintDebug(Mstring2.str(), SIMPlog, 1);

         DoIprint1 = false;
      }

      if (DoIprint2) {
         std::ostringstream SCWstr;
         const time_t Ctime = std::time(0);
         SCWstr << ' ' << __LINE__ << " Detected change in Wstat "
                << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
         PrintDebug(SCWstr.str(), MASlog, 1);
         printsimplex(blob, MASlog);
         DoIprint2 = false;
      }

      // don't use 100% CPU in a busy loop
      std::this_thread::sleep_for(std::chrono::seconds(1));
   }  // end while

   PrintDebug("Closing file...", MASlog, 5);
   PrintDebug("Closing file...", SIMPlog, 5);
   MASlog.close();
   SIMPlog.close();
}  // end MASTER

//Debug functions
void
CheckProcID(const Notsimplex blob){

   //check ProcIDs between vertex and auxR
   for(arma::uword i = 0; i < blob.n_vert; i++){
      for(arma::uword j = 0; j < blob.fSize; j++){
         if(blob.vertex.at(i).ProcID == blob.auxR.at(j).ProcID){
            std::cerr << "Identical ProcID between vertex" << i 
                      << " Proc" << blob.vertex.at(i).ProcID
                      << " and auxR" << j << " Proc" << blob.auxR.at(j).ProcID << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);   
         }
      }
   }

   //check ProcIDs between vertex and aux2
   for(arma::uword i = 0; i < blob.n_vert; i++){
      for(arma::uword j = 0; j < blob.fSize; j++){
         if(blob.vertex.at(i).ProcID == blob.aux2.at(j).ProcID){
            std::cerr << "Identical ProcID between vertex" << i 
                      << " Proc" << blob.vertex.at(i).ProcID
                      << " and aux2" << j << " Proc" << blob.aux2.at(j).ProcID << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);   
         }
      }
   }


   //check ProcIDs between all the vertices
   for(arma::uword i = 0; i < blob.n_vert; i++){
      for(arma::uword j = i+1; j < blob.n_vert; j++){
         if(blob.vertex.at(i).ProcID == blob.vertex.at(j).ProcID){
            std::cerr << "Identical ProcID between vertex" << i 
                      << " Proc" << blob.vertex.at(i).ProcID
                      << " and vertex" << j << " Proc" << blob.vertex.at(j).ProcID << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);   
         }
      }
   }

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
