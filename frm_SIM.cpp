/*
 * This file is a source code of the simpelx algorithm, which 99% of the
 * code is written by my advisor long time ago.
 */


#include "frame.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>


// function used for std::sort. compares two vertex points
bool
POINTcompare(point a, point b)
{
   return a.value < b.value;
}

/*
  reflect inverts the worst point of a simplex through its base,
  storing the resulting point in aux[0]....sjs 6/25/09
*/
void
reflect(simplex *blob)
{
   int idim, status;

   // no one should call reflect() without knowing which point is worst,
   // and where the base is, but we'd better be sure
   if (stale == blob->sort_status) {
      //XXX
      GlobalPrintDebug("Order2 called reflect", 1);
      status = order2(blob, true);
      if (status < 0) {  // means some vertices do not have values yet
         std::fprintf(stderr, "reflect: order returned %d", status);
//         exit(1);
      }
      getbase(blob);
   }

   blob->aux[0].dim = blob->vertex[blob->highest].dim;  // just in case

   for (idim = 0; idim < blob->aux[0].dim; ++idim) {
      blob->aux[0].coord[idim] = blob->vertex[blob->highest].coord[idim]
         + 2.0 * (blob->base[idim] - blob->vertex[blob->highest].coord[idim]);
   }
   blob->aux[0].value = -1.0;
   blob->aux[0].error = -1.0;
}

/*
  extend inverts the worst point of a simplex through its base and
  extends it farther, storing the resulting point in aux[1]....sjs 6/25/09
*/
void
extend(simplex *blob)
{
   // no one should call extend() without knowing which point is worst,
   // and where the base is, but we'd better be sure
   if (stale == blob->sort_status) {
      //XXX
      GlobalPrintDebug("Order2 called extend", 1);
      int status = order2(blob, true);
      if (status < 0) {  // means some vertices do not have values yet
         std::fprintf(stderr, "extend: order returned %d", status);
//         exit(1);
      }
      getbase(blob);
   }

   blob->aux[1].dim = blob->vertex[blob->highest].dim;  // just in case

   for (int idim = 0; idim < blob->aux[1].dim; ++idim) {
      blob->aux[1].coord[idim] = blob->vertex[blob->highest].coord[idim]
         + 3.0 * (blob->base[idim] - blob->vertex[blob->highest].coord[idim]);
   }
   blob->aux[1].value = -1.0;
   blob->aux[1].error = -1.0;
   blob->aux_holds = extension;
}

/*
  contract moves the worst point of a simplex halfway towards its base,
  storing the resulting point in aux[1]....sjs 6/26/09
*/
void
contract(simplex *blob)
{
   // no one should call contract() without knowing which point is worst,
   // and where the base is, but we'd better be sure
   if (stale == blob->sort_status) {
      //XXX
      GlobalPrintDebug("Order2 called contract", 1);
      int status = order2(blob, true);
      if (status < 0) {  // means some vertices do not have values yet
         std::fprintf(stderr, "contract: order returned %d", status);
//         exit(1);
      }
      getbase(blob);
   }

   blob->aux[1].dim = blob->vertex[blob->highest].dim;  // just in case

   for (int idim = 0; idim < blob->aux[1].dim; ++idim) {
      blob->aux[1].coord[idim] = blob->vertex[blob->highest].coord[idim]
         + 0.5 * (blob->base[idim] - blob->vertex[blob->highest].coord[idim]);  // XXX original value = 0.5 new = 0.8
   }
   blob->aux[1].value = -1.0;
   blob->aux[1].error = -1.0;
   blob->aux_holds = contraction;
}

/*
  collapse moves all of the points of a simplex halfway towards the
  best point...sjs 6/26/09
*/
void
collapse(simplex *blob)
{
   // no one should call collapse() without knowing which point is best,
   // but we'd better be sure
   if (stale == blob->sort_status) {
      //XXX
      GlobalPrintDebug("Order2 called collapse", 1);
      int status = order2(blob, true);
      if (status < 0) {  // means some vertices do not have values yet
         std::fprintf(stderr, "collapse: order returned %d", status);
//         exit(1);
      }
      // we don't need the base here, but the protocol is to always
      // calculate the base when resorting; base status is assumed
      // elsewhere to be identical to sort_status
      getbase(blob);
   }

   for (int ivert = 0; ivert < blob->vertices; ++ivert) {
      if (ivert != blob->lowest) {
         for (int idim = 0; idim < blob->vertex[ivert].dim; ++idim) {
            blob->vertex[ivert].coord[idim] = blob->vertex[ivert].coord[idim]
               + 0.5 * (blob->vertex[blob->lowest].coord[idim]  //XXX original value = 0.5 new = 0.8
                        - blob->vertex[ivert].coord[idim]);
         }
         blob->vertex[ivert].value = -1.0;
         blob->vertex[ivert].error = -1.0;
      }
   }

   blob->sort_status = stale;
   blob->reflect_status = stale;
   blob->aux_holds = none;
}

/*
  slide will return a point that has been translated from start through
  target by distance multiples of the distance between target and start.
  ...sjs 6/23/09
*/
point
slide(point start, point target, double distance)
{
   point finish;

   if (start.dim != target.dim) {
      std::fprintf(stderr, "slide: dimensionality of points do not match\n");
//      exit(1);
   }
   finish.dim = start.dim;

   for (int idim = 0; idim < finish.dim; ++idim) {
      finish.coord[idim] = start.coord[idim]
         + distance * (target.coord[idim] - start.coord[idim]);
   }
   finish.value = 0;
   finish.error = 0;
//   finish.status = undef;

   return finish;
}

/*
  order sorts the vertices of a simplex well enough to find out a partial
  ordering.
  if the function values are not available at some vertices, the return
  status is negative.
  ...sjs 6/23/09
  sort by average. Does NOT account for any error when ordering
*/

int
order1(simplex *blob)
{
   for (int ivert = 0; ivert < blob->vertices; ++ivert) {
      if (undef == blob->vertex[ivert].status) {
         return -1;
      }
   }
   // look for lowest value

   blob->lowest = 0;
   double low = blob->vertex[blob->lowest].value;
   for (int ivert = 1; ivert < blob->vertices; ++ivert) {
      if (blob->vertex[ivert].value < low) {
         blob->lowest = ivert;
         low = blob->vertex[blob->lowest].value;
      }
   }

   // look for highest value

   blob->highest = 0;
   double high = blob->vertex[blob->highest].value;
   for (int ivert = 1; ivert < blob->vertices; ++ivert) {
      if (blob->vertex[ivert].value > high) {
         blob->highest = ivert;
         high = blob->vertex[blob->highest].value;
      }
   }

   // look for second-highest value

   blob->sechigh = blob->highest == 0 ? 1 : 0;
   double sechigh = blob->vertex[blob->sechigh].value;
   for (int ivert = blob->sechigh + 1; ivert < blob->vertices; ++ivert) {
      if (blob->vertex[ivert].value > sechigh && ivert != blob->highest) {
         blob->sechigh = ivert;
         sechigh = blob->vertex[blob->sechigh].value;
      }
   }

   blob->sort_status = current;
//   std::cerr << "should verify order change before invalidating auxes" << std::endl;
   blob->reflect_status = stale;
   blob->aux_holds = none;

   return 0;
}

int
order2(simplex *blob, bool reset)
{
   // sort each point in order of blob->vertex from low to high
   // lowest index = lowest value (I think)
   std::cout << "FRAME SIM TEST: blob->vertex[] BEFORE SORT: ";
   for (int i = 0; i < blob->vertices; i++) {
      std::cout << blob->vertex[i].value << ' ';
   }
   std::cout << std::endl;

   std::cout << "FRAME SIM TEST: blob->vertex[] BEFORE SORT: ";
   for (int i = 0; i < blob->vertices; i++) {
      std::cout << blob->vertex[i].ProcID << ' ';
   }
   std::cout << std::endl;

   std::sort(blob->vertex, blob->vertex + blob->vertices, POINTcompare);
   // look for lowest value
   blob->lowest = 0;
   // look for highest value
   blob->highest = blob->vertices - 1;
   // look for second-highest value
   blob->sechigh = blob->vertices - 2;

   std::cout << "FRAME SIM TEST  blob->vertex[] AFTER SORT: ";
   for (int i = 0; i < blob->vertices; i++) {
      std::cout << blob->vertex[i].value << ' ';
   }
   std::cout << std::endl;

   std::cout << "FRAME SIM TEST  blob->vertex[] AFTER SORT: ";
   for (int i = 0; i < blob->vertices; i++) {
      std::cout << blob->vertex[i].ProcID << ' ';
   }
   std::cout << std::endl;

   std::cout << "FRAME SIM: Vert highest = " << blob->highest << std::endl;
   std::cout << "FRAME SIM: Proc highest = " << blob->vertex[blob->highest].ProcID << std::endl;
   std::cout << "FRAME SIM: Vert lowest = " << blob->lowest << std::endl;
   std::cout << "FRAME SIM: Proc lowest = " << blob->vertex[blob->lowest].ProcID << std::endl;
   std::cout << "FRAME SIM: Vert sechigh = " << blob->sechigh << std::endl;
   std::cout << "FRAME SIM: Proc sechigh = " << blob->vertex[blob->sechigh].ProcID << std::endl;

   if (reset) {
      blob->sort_status = current;
      blob->reflect_status = stale;
      blob->aux_holds = none;
   }

   return 0;
}

// getbase finds the centroid of all points in a simplex except the worst
void
getbase(simplex *blob)
{
   if (blob->sort_status != current) {
      //XXX
      GlobalPrintDebug("Order2 called getbase",1);
      int status = order2(blob, true);
      if (status < 0) {
         std::fprintf(stderr, "getbase: not all vertices have been evaluated.\n");
//         exit(1);
      }
   }

   for (int idim = 0; idim < blob->vertex[0].dim; ++idim) {
      blob->base[idim] = 0.0;
   }

   // sum all of the points except the worst

   for (int ivert = 0; ivert < blob->vertices; ++ivert) {
      if (ivert != blob->highest) {
         for (int idim = 0; idim < blob->vertex[ivert].dim; ++idim) {
            blob->base[idim] += blob->vertex[ivert].coord[idim];
         }
      }
   }

   // and divide by the number of vertices summed

   std::ostringstream message("base is at ");
   for (int idim = 0; idim < blob->vertex[0].dim; ++idim) {
      blob->base[idim] /= (blob->vertices - 1);
      message << blob->base[idim] << " ";
   }
   message << std::endl;
   report(message.str(), simplexStatus);
}

void
report(std::string message, verbosity level)
{
   extern verbosity outputLevel;
   if (level <= outputLevel) {
      std::cout << message;
   }
}

int
compare(point PntA, point PntB, double mult)
{
   // A > B w/o overlapping
   if (PntA.value - mult*PntA.error > PntB.value + mult*PntB.error) {
      return 1;
   }

   // A < B w/o overlapping
   if (PntA.value + mult*PntA.error < PntB.value - mult*PntB.error) {
      return 2;
   }

   // return this if error bar is overlapping
   return 0;
}

/*
 vN = vertex number to be swapped
 aN = aux number that is replacing vN
 Once the replacements happen, all of the aux points are not worth keeping
 so both is set to infected any time swap happens
*/
void
swapPoint(simplex &blob, int vN, int aN)
{
   point temp = blob.vertex[vN];
   blob.vertex[vN] = blob.aux[aN];
   blob.vertex[vN].Pkind = verte;  // change Pkind to vert
   blob.aux[aN] = temp;
   blob.aux[aN].Pkind = auxil;

   std::cout << "FRAME Point vertex " << vN << " is replaced with aux " << aN << std::endl;
}


//XXX THE ONLY FXN THAT MODIFIES ComEx & status OUTSIDE OF MASTER
void
resetRound(simplex &blob)
{
   if (blob.aux[1].status != undef) {
      blob.aux[1].Pkind = auxil;  // change Pkind to aux
      blob.aux[1].status = infected;  // keep it infected until new simulation starts
      blob.aux[1].ComEx = contagious;  // keep it infected until new simulation starts
   }

   blob.aux[0].Pkind = auxil;  // change Pkind to aux->fail safe
   blob.aux[0].status = infected;  // keep it infected until new simulation starts
   blob.aux[0].ComEx = contagious;  // keep it infected until new simulation starts
   blob.sort_status = stale;   // probably, but maybe not
   blob.reflect_status = stale;
   blob.aux_holds = none;
   ++blob.age;
}
///////////////////////////simplex action code//////////////////////////

// core logic for a deterministic simplex.  simplex is evaluated here, and
// action is returned to be performed externally....sjs 6/24/09

enum simplex_action
optimize(simplex *blob, point **target, bool &SIMPswap, double mult)
{
   // quit when all simplex edges are shorter than this tolerance
   // convergence criteria
   const double tolerance = 1.e-4;
   // or when the simplex has made this many moves
   const int old_age = 10000;

   while (1) {  // we'll come back here after vertex replacements
      // first, check to make sure we're not due to be euthanized

      if (blob->age >= old_age) {
         const time_t Ctime = std::time(0);
         std::stringstream messString;
         messString << __FILE__ << ' ' << __LINE__ << "reached maximum "
                    << old_age << " generations: Quitting... "
                    << std::asctime( std::localtime(&Ctime) ) << std::endl;
         GlobalPrintDebug(messString.str(),5 );

         return quit;
      }

      // first, check to make sure evaluation has begun at all vertices

      for (int ivert = 0; ivert < blob->vertices; ++ivert) {
         if (undef == blob->vertex[ivert].status) {
            *target = &blob->vertex[ivert];
            std::stringstream messString;
            messString << __FILE__ << ' ' << __LINE__ << " initiated vertex "
                       << ivert;
            GlobalPrintDebug(messString.str(),5 );
            return initiate;
         }
      }

      // next, make sure we have values at each vertex
      for (int ivert = 0; ivert < blob->vertices; ++ivert) {
         if (pending == blob->vertex[ivert].status) {
//            std::cout <<"FRAME  "  << __FILE__ << ' ' << __LINE__ << "vertex " 
//                      << ivert << " is pending" << std::endl;
            return simplex_wait;
         }
      }

      // next, make sure we have VALID values at each vertex
      for (int ivert = 0; ivert < blob->vertices; ++ivert) {
         if (blob->vertex[ivert].value < 0) {
            const time_t Ctime = std::time(0);
            std::cout <<"FRAME  " << __FILE__ << ' ' << __LINE__
                      << "Proc " << blob->vertex[ivert].ProcID
                      << " -Vertex " << ivert
                      << " Does not have valid value yet "
                      << std::asctime( std::localtime(&Ctime) ) << std::endl;
            return simplex_wait;
         }
      }

      // make sure the vertices are sorted
      if (stale == blob->sort_status) {
         //XXX
         GlobalPrintDebug("Order2 called optimize",1);
         int O_status = order2(blob, true);
         if (O_status < 0) {  // means status is undefined for some
            const time_t Ctime = std::time(0);
            std::cerr << "optimize: order returned " << O_status << std::endl
                      << std::asctime( std::localtime(&Ctime) );
//            exit(1)
            return crash;
         }
//	 std::cerr << "should have checked whether order changed, but didn't" << std::endl;
         // are we done yet?
         bool simplex_converged = true;
         for (int ivert = 0; ivert < blob->vertices; ++ivert) {
            for (int jvert = ivert+1; jvert < blob->vertices; ++jvert) {
            //XXX put distance back to original place when debug is over!
               double distance = 0.;
               for (int idim = 0; idim < blob->vertex[0].dim; ++idim) {
                  distance += std::pow( (blob->vertex[ivert].coord[idim]
                                  - blob->vertex[jvert].coord[idim])/blob->base[idim], 2.0);
               }
               if (distance > std::pow(tolerance, 2.0)) {
                  simplex_converged = false;
               }
         //XXX these are debug statements to track the termination criteria/distance
         //TODO debug and delete these print statements
            std::stringstream messString;
            messString << "distance = " << distance  << std::endl;
            GlobalPrintDebug(messString.str(),5 );
            }
         }
         if (simplex_converged) {
            // shut down all simulations
            const time_t Ctime = std::time(0);
            std::stringstream messString;
            messString << __FILE__ << ' ' << __LINE__ << " simplex converged! "
                       << std::asctime( std::localtime(&Ctime) );
            GlobalPrintDebug(messString.str(),5 );
            SIMPswap = false;  // privents extra printing
            return quit;
         }
         // calculate the base of the simplex
         getbase(blob);
         // calculate the reflected point
         reflect(blob);
         blob->reflect_status = current;
         *target = &blob->aux[0];
         if (undef == blob->aux[0].status || done == blob->aux[0].status) {
            const time_t Ctime = std::time(0);
            std::stringstream messString;
            messString << __FILE__ << ' ' << __LINE__ << " initiated Refletion "
                       << std::asctime( std::localtime(&Ctime) );
            GlobalPrintDebug(messString.str(),5 );

            return initiate;
         }
         else {  // something is already running in the reflection slot
            const time_t Ctime = std::time(0);
            std::stringstream messString;
            messString << __FILE__ << ' ' << __LINE__ << " REinitiated Refletion "
                       << std::asctime( std::localtime(&Ctime) );
            GlobalPrintDebug(messString.str(),5 );

            return restart;
         }
      }  // end if (blob->sort_status == stale)

      // at this point, we should have a valid reflection
      if (blob->reflect_status != current) {
         std::cerr << "optimize: why don't we have a valid reflection?" << std::endl;
//         exit(1);
      }
      // and the calculation should be underway
      if (undef == blob->aux[0].status) {
         std::cerr << "optimize: why hasn't reflection started?" << std::endl;
//         exit(1);
      }

      // but we can't do anything until we have a value for the reflection
      if (pending == blob->aux[0].status) {
         const time_t Ctime = std::time(0);
         std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                   << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
         return simplex_wait;
      }
      if (blob->aux[0].value < 0) {  //ensure the valid value on aux[0] before anything
         const time_t Ctime = std::time(0);
         std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                   << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
         return simplex_wait;
      }

      //------------------------start comparisons---------------------//
      /*
        list of cases:
        -NONOVERLAPPING for reflected points -> make decision
        1. reflect + error < lowest - error
        2. reflect - error > lowest + error && reflect + error < sechigh - error
        3. reflect - error > highest + error
        4. reflect - error > sechigh + error && reflect + error < highest - error

        -overlapping case for reflected points -> return wait
        1. overlapping w/ lowest but NOT sechigh -> probably lower than sechigh
        2. overlapping w/ sechigh but NOT highest -> between sechigh and highest
        3. overlapping w/ highest ONLY -> probably highest???
        4. overpapping w/ all important points -> indeterminate

        TODO !!!!! This is the major additon!
        -important points with everything else -> return wait
        1. lowest overlapping w/ aggirgate
        2. highest overlappign w/ aggrigate
        3. sechigh overlapping w/ aggrigate
      */

      //TODO TODO TODO XXX move to order or new fxn, eventually

      // compare highest and sechigh
      int chk1 = compare(blob->vertex[blob->sechigh], blob->vertex[blob->highest], mult);
      if (0 == chk1) {
         std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                   << ' ' << blob->vertex[blob->sechigh].value << " +/- " 
                   << blob->vertex[blob->sechigh].error << ' '
                   << " ?< " << blob->vertex[blob->highest].value << " +/- "
                   << blob->vertex[blob->highest].error << ' ';
         const time_t Ctime = std::time(0);
         std::cout << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
         return simplex_wait;
      }

      // Now chekc important points with non-important (middle) values
      int chk2 = -1;  // used inside the for loop below
      int chk3 = -1;  // used inside the for loop below
      int chk4 = -1;  // used inside the for loop below
      bool chkT = false;  // used to indicate order change on lowest & sechigh
      bool chkT2 = false;  // used to indicate order change on highest point
      for (int i = 0; i < blob->vertices; i++) {
         // check lowest point
         if (i != blob->lowest) {
            chk3 = compare(blob->vertex[blob->lowest], blob->vertex[i], mult);
            if (1 == chk3) {
               chkT = true;
            }
         }

         // check highest and sechigh point
         if (i != blob->sechigh && i != blob->highest) {
            chk2 = compare(blob->vertex[blob->sechigh], blob->vertex[i], mult);
            chk4 = compare(blob->vertex[blob->highest], blob->vertex[i], mult);
            if (chk2 == 2) {
               chkT = true;
            }
            if (chk4 == 2) {
               chkT2 = true;
            }
         }

         // XXX if any of chk returns 0, wait
         /*
           if(chk2 == 0 || chk3 ==0 || chk4 == 0){
           const time_t Ctime = std::time(0);
           std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
           << "chk2 = " << chk2 << " : chk3 = " << chk3 << " : chk4 = " << chk4
           << ' ' << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
           return simplex_wait;
           }
         */

         //XXX DEBUG VERSION of above
         if (0 == chk2) {
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait -- "
                      << "chk2 == 0 : " << blob->vertex[blob->sechigh].value << " +/- " << blob->vertex[blob->sechigh].error
                      << " <? " << blob->vertex[i].value << " +/- " << blob->vertex[i].error
                      << ' ' << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
         if (0 == chk3) {
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait -- "
                      << "chk3 == 0 : " << blob->vertex[blob->lowest].value << " +/- " << blob->vertex[blob->lowest].error
                      << " >? " << blob->vertex[i].value << " +/- " << blob->vertex[i].error
                      << ' ' << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
         if (0 == chk4) {
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait -- "
                      << "chk4 == 0 : " << blob->vertex[blob->highest].value << " +/- " << blob->vertex[blob->highest].error
                      << " <? " << blob->vertex[i].value << " +/- " << blob->vertex[i].error
                      << ' ' << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
      }  // end for(int i = 0; i < blob->vertices; i++)

      if (chk1 == 1 || chkT2) {
         //XXX IMPORTANT! value changed and sechigh is no
         // longer sechigh. since aux[0] is ready at this point,
         // go ahead and reorder. Then compare with new highest.
         // reorder and compare R with NEW-highest
         order2(blob, false);
         if (compare(blob->aux[0], blob->vertex[blob->highest], mult) == 2) {
            GlobalPrintDebug("keeping reflected point-0", 5);
            SIMPswap = true;//tells calling fxn that simplex has moved
            swapPoint(*blob, blob->highest, 0);
         }
         else {
            GlobalPrintDebug("chk1 == 1, and R is worse tha NEW-high", 5);
         }
         resetRound(*blob);

         const time_t Ctime = std::time(0);
         std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                   << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
         return simplex_wait;
      }

      if (chkT) {
         // Case where lowest is no longer lowest and sechigh is no longer sec high
         // Since these two only changes the comparison coming up, and not new point(s)
         // ie, point of reflection stays the same.Therefore, order it without
         // restarting the reflection
         // reorder and go on to next steps
         order2(blob, false);
      }

      if (-1 == chk2 || -1 == chk3 || -1 == chk4) {
         std::cout << "FRAME SIM ERROR: MOOSE" << std::endl;
         const time_t Ctime = std::time(0);
         std::cout << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
         // no actioin yet
      }

      // NOW ALL POINTS ARE READY TO BE COMPARED WITH REFLECTED POINTS NORMALLY
      //
      //
      // 2. reflected point is greater than lowest vale, but lower than 2nd highest
      if (compare(blob->aux[0], blob->vertex[blob->lowest], mult) == 1 &&
          compare(blob->aux[0], blob->vertex[blob->sechigh], mult) == 2) {

         // keep the reflected point.
         GlobalPrintDebug("keeping relfected point-1",5 );
         SIMPswap = true;  // tells calling fxn that simplex has moved
         swapPoint(*blob, blob->highest, 0);
         resetRound(*blob);

         const time_t Ctime = std::time(0);
         std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                   << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
         return simplex_wait;
      }
      // 1. reflected value is less than or equal to lowest value
      else if (compare(blob->aux[0], blob->vertex[blob->lowest], mult) == 2) {
         std::stringstream whatever;
         whatever << blob->aux[0].value << " +/- " << blob->aux[0].error << " < "
                  << blob->vertex[blob->lowest].value << " +/- " << blob->vertex[blob->lowest].error
                  << " Exension? branch" << std::endl;
         GlobalPrintDebug(whatever.str(), 5);
         // start a calculation of the extended point, if needed
         if (blob->aux_holds != extension) {
            extend(blob);
            *target = &blob->aux[1];
            if (undef == blob->aux[1].status || done == blob->aux[0].status) {
               std::stringstream messString;
               messString << __FILE__ << ' ' << __LINE__ << " Initiated extension";
               GlobalPrintDebug(messString.str(),5 );
               return initiate;
            }
            else {  // a simulation is currently running
               std::stringstream messString;
               messString << __FILE__ << ' ' << __LINE__ << " REinitiated extension";
               GlobalPrintDebug(messString.str(),5 );
               return restart;
            }
         }

         // at this point aux[1] should hold an extended point
         // and the calculation should be started
         if (undef == blob->aux[1].status) {
            std::cerr << "optimize: why hasn't extension started?" << std::endl;
//            exit(1);
         }

         // wait for a value, if needed
         if (pending == blob->aux[1].status) {
            const time_t Ctime = std::time(0);
            std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }

         if (blob->aux[1].value < 0) {
            const time_t Ctime = std::time(0);
            std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }

         // at this point we know the VALID value at the extended point
         if (compare(blob->aux[1], blob->aux[0], mult) == 2) {
            // keep extended point
            GlobalPrintDebug("keeping extended point", 5);
            SIMPswap = true;  // tells calling fxn that simplex has moved
            swapPoint(*blob, blob->highest, 1);
            resetRound(*blob);
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  " << __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
         else if (compare(blob->aux[1], blob->aux[0], mult) == 1) {
            // keep reflected point
            GlobalPrintDebug("keeping reflected point-2", 5);
            SIMPswap = true;  // tells calling fxn that simplex has moved
            swapPoint(*blob, blob->highest, 0);
            resetRound(*blob);
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
         else {  // maybe branch
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
      }  // end else if (blob->aux[0].value <= blob->vertex[blob->lowest].value)
      //3.  if reflected point is higher than worst
      else if (compare(blob->aux[0], blob->vertex[blob->highest], mult) == 1) {  //TODO 2 types of contraction??? see drawing
         // start a calculation of the contracted point if needed
         if (blob->aux_holds != contraction) {
            contract(blob);
            *target = &blob->aux[1];
            if (undef == blob->aux[1].status || done == blob->aux[1].status) {
               std::stringstream messString;
               messString << __FILE__ << ' ' << __LINE__ << " initiated Contraction";
               GlobalPrintDebug(messString.str(),5 );
               return initiate;
            }
            else {  // simulation is already running at this vertex
               std::stringstream messString;
               messString << __FILE__ << ' ' << __LINE__ << " REinitiated Contraction";
               GlobalPrintDebug(messString.str(),5 );
               return restart;
            }
         }

         // at this point the contraction should be valid, and running
         if (undef == blob->aux[1].status) {
            std::cerr << "optimize: why is contraction not running?" << std::endl;
//            exit(1);
         }

         // but we may have to wait for a value
         if (pending == blob->aux[1].status) {
            const time_t Ctime = std::time(0);
            std::stringstream messString;
            messString <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                       << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            GlobalPrintDebug(messString.str(),5 );
            return simplex_wait;
         }
         // but we may have to wait for a VALID value
         if (blob->aux[1].value < 0) {
            const time_t Ctime = std::time(0);
            std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }

         // at this point we know the value at the contraction
         int OOI = compare(blob->aux[1], blob->aux[0], mult);
         int OOT = compare(blob->aux[1], blob->vertex[blob->highest], mult);
         if (0 == OOI || 0 == OOT) {
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  " << __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            std::stringstream whatever;
            whatever << blob->aux[0].value << " +/- " << blob->aux[0].error << " >  "
                     << blob->vertex[blob->highest].value << " +/- " << blob->vertex[blob->highest].error
                     << " Contraction1? branch"<< std::endl;
            GlobalPrintDebug(whatever.str(),5 );
            return simplex_wait;
         }

         if (2 == OOI && 2 == OOT) {
            // keep contracted
            // replace highests with contracted point
            GlobalPrintDebug("keeping contracted point 1",5 );
            SIMPswap = true;//tells calling fxn that simplex has moved
            swapPoint(*blob, blob->highest, 1);
            resetRound(*blob);
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
         else {  // contracted point was no better than the anemic reflection
            // when all else fails, collapse towards the best point
            GlobalPrintDebug("collapsing-0", 5);
            collapse(blob);
            SIMPswap = true;  // tells calling fxn that simplex has moved
            resetRound(*blob);
            return shrink;
         }
      }
      // 4. reflected point is worse than the second worst, but better than highest
      else if (compare(blob->aux[0], blob->vertex[blob->sechigh], mult) == 1
               && compare(blob->aux[0], blob->vertex[blob->highest], mult) == 2) {
        // 2nd type of contraction??? see drawing
        // replace the highest point with reflected point
         if (blob->aux_holds != contraction) {
            GlobalPrintDebug("Half-step in contraction", 5);
            // swap highest and R for breif moment while contraction point is calculated
            point temp = blob->vertex[blob->highest];
            blob->vertex[blob->highest] = blob->aux[0];
            // Now find contracted point
            contract(blob);
            *target = &blob->aux[1];
            // swap back to original points so that simplex can find 
            // way back here on next cycle-Will swap for last time if new cycle comes back here
            // and attempt to replace/compare contracted point
            blob->vertex[blob->highest] = temp;

            if (undef == blob->aux[1].status || done == blob->aux[1].status) {
               std::stringstream messString;
               messString << __FILE__ << ' ' << __LINE__ << " initiated Contraction";
               GlobalPrintDebug(messString.str(), 5);
               return initiate;
            }
            else {
               // simulation is already running at this vertex
               std::stringstream messString;
               messString << __FILE__ << ' ' << __LINE__ << " REinitiated Contraction";
               GlobalPrintDebug(messString.str(), 5);
               return restart;
            }
         }

         // at this point the contraction should be valid, and running
         if (undef == blob->aux[1].status) {
            std::cerr << "optimize: why is contraction not running?" << std::endl;
//            exit(1);
         }

         // but we may have to wait for a value
         if (pending == blob->aux[1].status) {
            const time_t Ctime = std::time(0);
            std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
         // but we may have to wait for a VALID value
         if (blob->aux[1].value < 0) {
            const time_t Ctime = std::time(0);
            std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }

         int OOI = compare(blob->aux[1], blob->aux[0], mult);
         if (0 == OOI) {
            // error bar is overlapping
            const time_t Ctime = std::time(0);
            std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            std::stringstream whatever;
            whatever << blob->vertex[blob->sechigh].value << blob->vertex[blob->sechigh].error
                     << " < " << blob->aux[0].value << " +/- " << blob->aux[0].error << " < " 
                     <<  blob->vertex[blob->highest].value << " +/- "<< blob->vertex[blob->highest].error
                     << " \nContraction2? branch"<< std::endl;
            GlobalPrintDebug(whatever.str(), 5);
            return simplex_wait;
         }

         // Put this here, so that proper comparison can happen. If this is not here,
         // then reflection will infinitely repeat (theoretically)
         // NOW replacing old highest with NEW highest--old H > new H
         // Also, the error bar comparison already happened on OOI
         swapPoint(*blob, blob->highest, 0);
         int OOT = compare(blob->aux[1], blob->vertex[blob->highest], mult);
         // NOW old high is replaced by aux[0]--VALUE OF aux[0] is the NEW high

         // at this point we know the value at the contraction
         // Also, COMPARING "NEW High" point with others
         if (2 == OOT) {
            // keep contracted C < NEW H
            GlobalPrintDebug("keeping contracted point 2", 5);
            SIMPswap = true;  // tells calling fxn that simplex has moved
            swapPoint(*blob, blob->highest, 1);
            resetRound(*blob);

            const time_t Ctime = std::time(0);
            std::cout << "FRAME  " << __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
         else {
            //Probably C > NEW H
            // when all else fails, collapse towards the best point
            // Also, this means something is not working properly
            GlobalPrintDebug("collapsing-1", 5);
            collapse(blob);
            resetRound(*blob);
            SIMPswap = true;  // tells calling fxn that simplex has moved
            return shrink;
         }
      }
      else {  // other cases-probably overlapping. So wait.
         const time_t Ctime = std::time(0);
         std::cout << "FRAME  " << __FILE__ << ' ' << __LINE__ << " simplex wait "
                   << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
         return simplex_wait;
      }
   }  // while(1)

   const time_t Ctime = std::time(0);
   std::cout << "FRAME  " << __FILE__ << ' ' << __LINE__ << " simplex wait "
             << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
   return simplex_wait;
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
