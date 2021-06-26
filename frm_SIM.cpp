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
reflect(Notsimplex *blob)
{
   int status;
   getbase(blob); //call getbase anyways to be safe

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

   //index for fHigh = blob->n_vert-blob.fSize
   //P(high) + 2*( P(bar) - P(High) )
   std::cout << "Reflected coord: " << '\t';

   for(arma::uword i = 0; i < blob->fSize; i++){
      for (arma::uword idim = 0; idim < blob->n_dim; ++idim) {
         blob->auxR.at(i).coord(idim) = blob->vertex.at(i + blob->n_vert - blob->fSize).coord(idim)
            + 2.0 * (blob->baseCent(idim) - blob->vertex.at(i + blob->n_vert - blob->fSize).coord(idim));

         std::cout << blob->auxR.at(i).coord(idim) << '\t';
      }
      std::cout << std::endl;

      blob->auxR.at(i).value = -1.0;
      blob->auxR.at(i).error = -1.0;
   }


}

/*
  extend inverts the worst point of a simplex through its base and
  extends it farther, storing the resulting point in aux[1]....sjs 6/25/09
*/
void
extend(Notsimplex *blob)
{
   // no one should call extend() without knowing which point is worst,
   // and where the base is, but we'd better be sure

   getbase(blob); //call getbase anyways to be safe for calculations

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

   for(arma::uword i = 0; i < blob->fSize; i++){
      for (arma::uword idim = 0; idim < blob->n_dim; ++idim) {
         blob->aux2.at(i).coord[idim] = blob->vertex.at(i + blob->n_vert - blob->fSize).coord(idim)
            + 3.0 * (blob->baseCent(idim) - blob->vertex.at(i + blob->n_vert - blob->fSize).coord(idim));

      }
      blob->aux2.at(i).value = -1.0;
      blob->aux2.at(i).error = -1.0;
   }

   blob->aux_holds = extension;
}

/*
  contract moves the worst point of a simplex halfway towards its base,
  storing the resulting point in aux[1]....sjs 6/26/09
*/
void
contract(Notsimplex *blob)
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

   getbase(blob);//call anyways to be safe

   //P(High) + 0.5*( P(bar) - P(High) )
   for(arma::uword i = 0; i < blob->fSize; i++){
      std::cout << "contract(): PHighs are: ";
      
      for (arma::uword idim = 0; idim < blob->n_dim; ++idim) {
         blob->aux2.at(i).coord(idim) = blob->vertex.at(i + blob->n_vert - blob->fSize).coord(idim)
            + 0.5 * (blob->baseCent(idim) - blob->vertex.at(i + blob->n_vert - blob->fSize).coord(idim));
         std::cout << blob->vertex.at(i + blob->n_vert - blob->fSize).coord(idim) << '\t';
      }     
      
      std::cout << std::endl;
      std::cout << "contract(): baseCent are: ";
      for (arma::uword idim = 0; idim < blob->n_dim; ++idim) {
         std::cout << blob->baseCent(idim) << '\t';
      
      }
      std::cout << std::endl;

      blob->aux2.at(i).value = -1.0;
      blob->aux2.at(i).error = -1.0;
   }

   blob->aux_holds = contraction;
}

/*
  collapse moves all of the points of a simplex halfway towards the
  best point...sjs 6/26/09
*/
void
collapse(Notsimplex *blob)
{
   // no one should call collapse() without knowing which point is best,
   // but we'd better be sure
   //
   getbase(blob); //call getbase anyways to be safe

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

   //these index is all the index other than lowest fraction
   for (arma::uword ivert = blob->fSize; ivert < blob->n_vert; ++ivert) {
         for (arma::uword idim = 0; idim < blob->n_dim; ++idim) {
            blob->vertex.at(ivert).coord(idim) = blob->vertex.at(ivert).coord(idim)
               + 0.5 * (blob->vertex.at(ivert - blob->fSize).coord(idim)
                        - blob->vertex.at(ivert).coord(idim));
         }
         blob->vertex.at(ivert).value = -1.0;
         blob->vertex.at(ivert).error = -1.0;
   }

   blob->sort_status = stale;
   blob->reflect_status = stale;
   blob->aux_holds = none;
}

/*
  order sorts the vertices of a simplex well enough to find out a partial
  ordering.
  if the function values are not available at some vertices, the return
  status is negative.
  ...sjs 6/23/09
  sort by average. Does NOT account for any error when ordering
*/
/*
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
}*/

int
order2(Notsimplex *blob, bool reset)
{
   // sort each point in order of blob->vertex from low to high
   // lowest index = lowest value (I think)
   std::cout << "FRAME SIM TEST: blob->vertex[].value BEFORE SORT: ";
   for (arma::uword i = 0; i < blob->n_vert; i++) {
      std::cout << blob->vertex[i].value << ' ';
   }
   std::cout << std::endl;

   std::cout << "FRAME SIM TEST: blob->vertex[].ProcID BEFORE SORT: ";
   for (arma::uword i = 0; i < blob->n_vert; i++) {
      std::cout << blob->vertex[i].ProcID << ' ';
   }
   std::cout << std::endl;

//   std::sort(blob->vertex, blob->vertex + blob->n_vert, POINTcompare);
   std::sort(blob->vertex.begin(), blob->vertex.end(), POINTcompare);

   std::cout << "FRAME SIM TEST  blob->vertex[].value AFTER SORT: ";
   for (arma::uword i = 0; i < blob->n_vert; i++) {
      std::cout << blob->vertex[i].value << ' ';
   }
   std::cout << std::endl;

   std::cout << "FRAME SIM TEST  blob->vertex[].ProcID AFTER SORT: ";
   for (arma::uword i = 0; i < blob->n_vert; i++) {
      std::cout << blob->vertex[i].ProcID << ' ';
   }
   std::cout << std::endl;

   if (reset) {
      blob->sort_status = current;
      blob->reflect_status = stale;
      blob->aux_holds = none;
   }

   return 0;
}

// getbase finds the centroid of all fractions
void
getbase(Notsimplex *blob)
{
   if (blob->sort_status != current) {
      GlobalPrintDebug("Order2 called getbase",1);
      int status = order2(blob, true);
      if (status < 0) {
         std::fprintf(stderr, "getbase: not all vertices have been evaluated.\n");
//         exit(1);
      }
   }

//XXX TODO take the propagation of error for the AVERAGE value for each fraction

//reinitialize the components of each of the vectors, for obvious reasons
         blob->baseCent.zeros();
         blob->baseLow.zeros();
         blob->baseHigh.zeros();
         blob->baseSecHigh.zeros();
         blob->yfCent = 0;
         blob->yfLow = 0;
         blob->yfSecHigh = 0;
         blob->yfHigh = 0;


// calculate the bases of the swarm that is not the fHigh. These are just the sum
   for (arma::uword i = 0;i < blob->n_vert - blob->fSize; i++){
      for (arma::uword j = 0; j < blob->n_dim; j++){
         blob->baseCent.at(j) += blob->vertex.at(i).coord(j);
      }
      blob->yfCent += blob->vertex.at(i).value;
   }
// must divide by n_dim to get the actual base value
// base of all vertex
   std::cout << "The baseCent is: ";
   for(arma::uword j = 0; j < blob->n_dim; j++){
      blob->baseCent(j) = blob->baseCent(j) / static_cast<double>(blob->n_vert - blob->fSize);
   }
// calculate the base of the fLow, along with sum of the values.
   for (arma::uword i = 0;i < blob->fSize; i++){
      for (arma::uword j = 0; j < blob->n_dim; j++){
         blob->baseLow(j) += blob->vertex.at(i).coord(j);
      }
      blob->yfLow += blob->vertex.at(i).value;
   }

   //calculate the average of baseLow
// must divide by n_dim to get the actual base value
// base of all vertex
   for(arma::uword j = 0; j < blob->n_dim; j++){
      blob->baseLow(j) = blob->baseLow(j) / static_cast<double>(blob->fSize);
   }

   //then calculate the average of yLow
   blob->yfLow = blob->yfLow / static_cast<double>(blob->fSize);

//sechigh is from blob.vertex[index] n_vert-1-fSize (upper) to n_vert-1-fSize-fSize
//        //calculate the center and y(sum) of fSecHigh
   for (arma::uword i = blob->n_vert-2*blob->fSize; i < blob->n_vert - blob->fSize;i++){
      for (arma::uword j = 0; j < blob->n_dim; j++){
         blob->baseSecHigh(j) += blob->vertex[i].coord(j);
      }
      blob->yfSecHigh += blob->vertex.at(i).value;
   }
   //then calculate the average of ySecHigh
   blob->yfSecHigh = blob->yfSecHigh / static_cast<double>(blob->fSize);

   //calculate the average of SecHigh
// must divide by n_dim to get the actual base value
// base of all vertex
   for(arma::uword j = 0; j < blob->n_dim; j++){
      blob->baseSecHigh(j) = blob->baseSecHigh(j) / static_cast<double>(blob->fSize);
   }
   
   //highest is from indx n_vert-1 (upper) to n_vert-1-fSize
   //calculate the center and y(sum) of fHigh
   for (arma::uword i = blob->n_vert - blob->fSize; i < blob->n_vert; i++){
      for (arma::uword j = 0; j < blob->n_dim; j++){
         blob->baseHigh.at(j) += blob->vertex.at(i).coord(j);
      }
      blob->yfHigh += blob->vertex.at(i).value;
   }
   //then calculate the average of yHigh
   blob->yfHigh = blob->yfHigh / static_cast<double>(blob->fSize);

// must divide by n_dim to get the actual base value
// base of all vertex
   for(arma::uword i = 0; i < blob->n_dim; i++){
      blob->baseHigh.at(i) = blob->baseHigh.at(i) / static_cast<double>(blob->fSize);
   }
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

int
compare2(double PntA, double ErrA, double PntB, double ErrB, double mult)
{
   //XXX TODO FIX ErrA and ErrB is set to 0 for debugging purpose XXX TODO
   ErrA = 0;
   ErrB = 0;

   // A > B w/o overlapping
   if (PntA - mult*ErrA > PntB + mult*ErrB) {
      return 1;
   }

   // A < B w/o overlapping
   if (PntA + mult*ErrA < PntB - mult*ErrB) {
      return 2;
   }

   // return this if error bar is overlapping
   return 0;
}


void
swapPoint(Notsimplex &blob, int aN)
{

   for(arma::uword i = 0; i<blob.fSize; i++ ){
   point temp = blob.vertex.at(i + blob.n_vert - blob.fSize);
      //if aN = 0, aux is auxR
      if(aN == 0){
         blob.vertex.at(i + blob.n_vert - blob.fSize) = blob.auxR.at(i);
         blob.vertex.at(i + blob.n_vert - blob.fSize).Pkind = verte;  // change Pkind to vert
         blob.auxR.at(i) = temp;
         blob.auxR.at(i).Pkind = auxil;
      }
   
      //if aN = 1, aux is aux2
      if(aN == 1){
         blob.vertex.at(i + blob.n_vert - blob.fSize) = blob.aux2.at(i);
         blob.vertex.at(i + blob.n_vert - blob.fSize).Pkind = verte;  // change Pkind to vert
         blob.aux2.at(i) = temp;
         blob.aux2.at(i).Pkind = auxil;
      }

   }

//   std::cout << "FRAME Point vertex " << vN << " is replaced with aux " << aN << std::endl;
}


void
resetRound(Notsimplex &blob)
{
   for(arma::uword i = 0; i < blob.fSize; i++){
      if (blob.aux2.at(i).status != undef) {
         blob.aux2.at(i).Pkind = auxil;  // change Pkind to aux
         blob.aux2.at(i).status = infected;  // keep it infected until new simulation starts
         blob.aux2.at(i).ComEx = contagious;  // keep it infected until new simulation starts
      }

      if (blob.auxR.at(i).status != undef) {
         blob.auxR.at(i).Pkind = auxil;  // change Pkind to aux->fail safe
         blob.auxR.at(i).status = infected;  // keep it infected until new simulation starts
         blob.auxR.at(i).ComEx = contagious;  // keep it infected until new simulation starts
      }
   }
   blob.sort_status = stale;   // probably, but maybe not
   blob.reflect_status = stale;
   blob.aux_holds = none;
   ++blob.age;
}
///////////////////////////simplex action code//////////////////////////

// core logic for a deterministic simplex.  simplex is evaluated here, and
// action is returned to be performed externally....sjs 6/24/09

enum simplex_action
optimize(Notsimplex *blob, std::vector<int> &target, bool &SIMPswap, double mult)
{
   //clear the target vector
   target.clear();

   //set the size of target to fSize
   //The only time the size of target is not fSize
   //is during initializaiton. Hence, this vector will be 
   //resized accordingly within the scope
   target.resize(blob->fSize);

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

      //change the size of the target vector for loop index in master()
      //target IS NOW A VECTOR OF INTEGER(S)
      //the size of target here is 1, just because this is targeting individual
      //point in the vertex, not any of the aux.
      for (arma::uword ivert = 0; ivert < blob->n_vert; ++ivert) {
         if (undef == blob->vertex.at(ivert).status) {
            target.resize(1);
            target.at(0) = blob->vertex.at(ivert).ProcID;
            
            std::stringstream messString;
            messString << __FILE__ << ' ' << __LINE__ << " initiated vertex "
                       << ivert;
            GlobalPrintDebug(messString.str(),5 );
            return initiate;
         }
      }

      // next, make sure we have values at each vertex
      for (arma::uword ivert = 0; ivert < blob->n_vert; ++ivert) {
         if (pending == blob->vertex.at(ivert).status) {
//            std::cout <<"FRAME  "  << __FILE__ << ' ' << __LINE__ << "vertex " 
//                      << ivert << " is pending" << std::endl;
            return simplex_wait;
         }
      }

      // next, make sure we have VALID values at each vertex
      for (arma::uword ivert = 0; ivert < blob->n_vert; ++ivert) {
         if (blob->vertex.at(ivert).value < 0) {
            const time_t Ctime = std::time(0);
            std::cout <<"FRAME  " << __FILE__ << ' ' << __LINE__
                      << "Proc " << blob->vertex.at(ivert).ProcID
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
         for (arma::uword ivert = 0; ivert < blob->n_vert; ++ivert) {
            for (arma::uword jvert = ivert+1; jvert < blob->n_vert; ++jvert) {
            //XXX put distance back to original place when debug is over!
               double distance = 0.;
               for (arma::uword idim = 0; idim < blob->n_dim; ++idim) {
                  distance += std::pow( (blob->vertex.at(ivert).coord(idim)
                                  - blob->vertex.at(jvert).coord(idim))/blob->baseCent(idim), 2.0);
               }
               if (distance > std::pow(tolerance, 2.0)) {
                  simplex_converged = false;
               }
         //XXX these are debug statements to track the termination criteria/distance
         //TODO debug and delete these print statements
            std::stringstream messString;
            messString << "distance = " << distance  << std::endl;
            GlobalPrintDebug(messString.str(),9);
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

         //flags for return statemens in the coming loop
         bool FLGinit = false;
         bool FLGreinit = false;

         //Only loop over upto fSize, since this is about one of the aux
         for(arma::uword i = 0; i < blob->fSize;i++){
            if (undef == blob->auxR.at(i).status || done == blob->auxR.at(i).status) {
               const time_t Ctime = std::time(0);
               std::stringstream messString;
               messString << __FILE__ << ' ' << __LINE__ << " Proc" << blob->auxR.at(i).ProcID
                          << " initiated Refletion " << std::asctime( std::localtime(&Ctime) );
               GlobalPrintDebug(messString.str(),5 );
               target.at(i) = blob->auxR.at(i).ProcID;
               FLGinit = true;
            }
            else {  // something is already running in the reflection slot
               const time_t Ctime = std::time(0);
               std::stringstream messString;
               messString << __FILE__ << ' ' << __LINE__ << " Proc" << blob->auxR.at(i).ProcID
                       << " REinitiated Refletion "
                       << std::asctime( std::localtime(&Ctime) );
               GlobalPrintDebug(messString.str(),5 );
               target.at(i) = blob->auxR.at(i).ProcID;
               FLGreinit = true;
            }
         }
         
         if(FLGinit){
            return initiate;
         }
         if(FLGreinit){
            return restart;
         }

      }  // end if (blob->sort_status == stale)

      // at this point, we should have a valid reflection
      if (blob->reflect_status != current) {
         std::cerr << "optimize: why don't we have a valid reflection?" << std::endl;
//         exit(1);
      }
      // and the calculation should be underway
//      if (undef == blob->aux[0].status) {
//         std::cerr << "optimize: why hasn't reflection started?" << std::endl;
//         exit(1);
//      }

      // but we can't do anything until we have a value for the reflection
      for(arma::uword i = 0; i < blob->fSize; i++){
         if (pending == blob->auxR.at(i).status) {
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                   << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }    
         if (blob->auxR.at(i).value < 0) {  //ensure the valid value on aux[0] before anything
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                   << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
      }

//TODO calculate all the average value of the fractions here TODO

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
      

      //call calcFrac so that all the yfrac values are ready
      //should have valid value for all vertex at this point
      
      calcFrac(*blob);

      // compare highest and sechigh
      int chk1 = compare2(blob->yfSecHigh, blob->ufSecHigh, blob->yfHigh, blob->ufHigh, mult);
      if (0 == chk1) {
         std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                   << ' ' << blob->yfSecHigh << " +/- " 
                   << blob->ufSecHigh << ' '
                   << " ?< " << blob->yfHigh << " +/- "
                   << blob->ufHigh << ' ';
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
         // check lowest point
         chk3 = compare2(blob->yfLow, blob->ufLow, blob->yfCent, blob->ufCent, mult);
         if (1 == chk3) {
               chkT = true;
         }
         // check highest and sechigh point
         chk2 = compare2(blob->yfSecHigh, blob->ufSecHigh, blob->yfCent, blob->ufCent, mult);
         chk4 = compare2(blob->yfHigh, blob->ufHigh, blob->yfCent, blob->ufCent, mult);
         if (chk2 == 2) {
            chkT = true;
         }
         if (chk4 == 2) {
            chkT2 = true;
         }

         // XXX if any of chk returns 0, wait
           if(chk2 == 0 || chk3 ==0 || chk4 == 0){
           const time_t Ctime = std::time(0);
           std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
           << "chk2 = " << chk2 << " : chk3 = " << chk3 << " : chk4 = " << chk4
           << ' ' << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
           return simplex_wait;
           }

         //XXX DEBUG VERSION of above
         if (0 == chk2) {
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait -- "
                      << "chk2 == 0 : " << blob->yfSecHigh
                      << " +/- " << blob->ufSecHigh
                      << " <? " << blob->yfCent << " +/- " << blob->ufCent
                      << ' ' << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
         if (0 == chk3) {
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait -- "
                      << "chk3 == 0 : " << blob->yfLow
                      << " +/- " << blob->ufLow
                      << " >? " << blob->yfCent << " +/- " << blob->ufCent
                      << ' ' << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
         if (0 == chk4) {
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait -- "
                      << "chk4 == 0 : " << blob->yfHigh << " +/- " 
                      << blob->ufHigh
                      << " <? " << blob->yfCent << " +/- " << blob->ufCent
                      << ' ' << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
 
     if (chk1 == 1 || chkT2) {
         //XXX IMPORTANT! value changed and sechigh is no
         // longer sechigh. since aux[0] is ready at this point,
         // go ahead and reorder. Then compare with new highest.
         // reorder and compare R with NEW-highest
         order2(blob, false);
         if (compare2(blob->yfR, blob->ufR, blob->yfHigh, blob->ufHigh, mult) == 2) {
            GlobalPrintDebug("keeping reflected point-0", 5);
            SIMPswap = true;//tells calling fxn that simplex has moved
            swapPoint(*blob, 0);
            std::cerr << "before CheckProcID "<< __FILE__ << ' ' << __LINE__ << std::endl;
            CheckProcID(*blob);
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
      if (compare2(blob->yfR, blob->ufR, blob->yfLow, blob->ufLow, mult) == 1 &&
          compare2(blob->yfR, blob->ufR, blob->yfSecHigh, blob->ufSecHigh, mult) == 2) {

         // keep the reflected point.
         GlobalPrintDebug("keeping relfected point-1",5 );
         SIMPswap = true;  // tells calling fxn that simplex has moved
         
         swapPoint(*blob, 0);
         std::cerr << "before CheckProcID "<< __FILE__ << ' ' << __LINE__ << std::endl;
         CheckProcID(*blob);

         resetRound(*blob);
         std::cerr << "before CheckProcID "<< __FILE__ << ' ' << __LINE__ << std::endl;
         CheckProcID(*blob);

         const time_t Ctime = std::time(0);
         std::cout << "FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                   << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
         return simplex_wait;
      }
      // 1. reflected value is less than or equal to lowest value
      else if (compare2(blob->yfR, blob->ufR, blob->yfLow, blob->ufLow, mult) == 2) {
         std::stringstream whatever;
         whatever << blob->yfR << " +/- " << blob->ufR << " < "
                  << blob->yfLow << " +/- " << blob->ufLow
                  << " Extension? branch" << std::endl;
         GlobalPrintDebug(whatever.str(), 5);
         // start a calculation of the extended point, if needed
         if (blob->aux_holds != extension) {
            extend(blob);
            bool FLGinit = false;
            bool FLGreinit = false;

            for(arma::uword i = 0; i < blob->fSize; i++){
               target.at(i) = blob->aux2.at(i).ProcID;
               if (undef == blob->aux2.at(i).status || done == blob->aux2.at(i).status) {
                  std::stringstream messString;
                  messString << __FILE__ << ' ' << __LINE__ << " Initiated extension";
                  GlobalPrintDebug(messString.str(),5 );
                  FLGinit = true;
               }
               else {  // a simulation is currently running
                  std::stringstream messString;
                  messString << __FILE__ << ' ' << __LINE__ << " REinitiated extension";
                  GlobalPrintDebug(messString.str(),5 );
                  FLGreinit = true;
               }
            }

            if(FLGinit){return initiate;}
            if(FLGreinit){return restart;}
         }
         // at this point aux[1] should hold an extended point
         // and the calculation should be started
         for(arma::uword i = 0; i < blob->fSize; i++){
            if (undef == blob->aux2.at(i).status) {
               std::cerr << "optimize: why hasn't extension started?" << std::endl;
//            exit(1);
            }
         
            // wait for a value, if needed
            if (pending == blob->aux2.at(i).status) {
               const time_t Ctime = std::time(0);
               std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                         << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
               return simplex_wait;
            }

            if (blob->aux2.at(i).value < 0) {
               const time_t Ctime = std::time(0);
               std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                         << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
               return simplex_wait;
            }
         }

         // at this point we know the VALID value at the extended point
         if (compare2(blob->yf2, blob->uf2, blob->yfR, blob->ufR, mult) == 2) {
            // keep extended point
            GlobalPrintDebug("keeping extended point", 5);
            SIMPswap = true;  // tells calling fxn that simplex has moved
            
            swapPoint(*blob, 1);
            std::cerr << "before CheckProcID "<< __FILE__ << ' ' << __LINE__ << std::endl;
            CheckProcID(*blob);
            
            resetRound(*blob);
            std::cerr << "before CheckProcID "<< __FILE__ << ' ' << __LINE__ << std::endl;
            CheckProcID(*blob);

            const time_t Ctime = std::time(0);
            std::cout << "FRAME  " << __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            return simplex_wait;
         }
         else if (compare2(blob->yf2, blob->uf2, blob->yfR, blob->ufR, mult) == 1) {
            // keep reflected point
            GlobalPrintDebug("keeping reflected point-2", 5);
            SIMPswap = true;  // tells calling fxn that simplex has moved

            swapPoint(*blob, 0);
            std::cerr << "before CheckProcID "<< __FILE__ << ' ' << __LINE__ << std::endl;
            CheckProcID(*blob);
            
            resetRound(*blob);
            std::cerr << "before CheckProcID "<< __FILE__ << ' ' << __LINE__ << std::endl;
            CheckProcID(*blob);

            
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
      else if (compare2(blob->yfR, blob->ufR, blob->yfHigh, blob->ufHigh, mult) == 1) {  
         //TODO 2 types of contraction??? see drawing
         // start a calculation of the contracted point if needed
         if (blob->aux_holds != contraction) {
            contract(blob);
            bool FLGinit = false;
            bool FLGreinit = false;

            for(arma::uword i = 0; i < blob->fSize; i++){
            target.at(i) = blob->aux2.at(i).ProcID;
               if (undef == blob->aux2.at(i).status || done == blob->aux2.at(i).status) {
                  blob->aux2.at(i).Pstat = not_running;
                  std::stringstream messString;
                  messString << __FILE__ << ' ' << __LINE__ << " initiated Contraction";
                  GlobalPrintDebug(messString.str(),5 );
                  FLGinit = true;
               }
               else {  // simulation is already running at this vertex
                  blob->aux2.at(i).Pstat = running;
                  std::stringstream messString;
                  messString << __FILE__ << ' ' << __LINE__ << " REinitiated Contraction";
                  GlobalPrintDebug(messString.str(),5 );
                  FLGreinit = true;
               }
            }

            if(FLGinit){return initiate;}
            if(FLGreinit){return restart;}
         }
         // at this point the contraction should be valid, and running
         for(arma::uword i = 0; i < blob->fSize; i++){
            if (undef == blob->aux2[i].status) {
               std::cerr << "optimize: why is contraction not running?" << std::endl;
//            exit(1);
            }
         }
         // but we may have to wait for a value
         for(arma::uword i = 0; i < blob->fSize; i++){
            if (pending == blob->aux2.at(i).status) {
               const time_t Ctime = std::time(0);
               std::stringstream messString;
               messString <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                        << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
               GlobalPrintDebug(messString.str(),5 );
               return simplex_wait;
            }
            // but we may have to wait for a VALID value
            if (blob->aux2.at(i).value < 0) {
               const time_t Ctime = std::time(0);
               std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                         << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
               return simplex_wait;
            }
         }

         // at this point we know the value at the contraction
         int OOI = compare2(blob->yf2, blob->uf2, blob->yfR, blob->ufR, mult);
         int OOT = compare2(blob->yf2, blob->uf2, blob->yfHigh, blob->ufHigh, mult);
         if (0 == OOI || 0 == OOT) {
            const time_t Ctime = std::time(0);
            std::cout << "FRAME  " << __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            std::stringstream whatever;
            whatever << blob->yfR << " +/- " << blob->ufR << " >  "
                     << blob->yfHigh << " +/- " << blob->ufHigh
                     << " Contraction1? branch"<< std::endl;
            GlobalPrintDebug(whatever.str(),5 );
            return simplex_wait;
         }

         if (2 == OOI && 2 == OOT) {
            // keep contracted
            // replace highests with contracted point
            GlobalPrintDebug("keeping contracted point 1",5 );
            SIMPswap = true;//tells calling fxn that simplex has moved
  
            swapPoint(*blob, 1);
            std::cerr << "before CheckProcID "<< __FILE__ << ' ' << __LINE__ << std::endl;
            CheckProcID(*blob);

            resetRound(*blob);
            std::cerr << "before CheckProcID "<< __FILE__ << ' ' << __LINE__ << std::endl;
            CheckProcID(*blob);

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
            std::cerr << "before CheckProcID "<< __FILE__ << ' ' << __LINE__ << std::endl;
            CheckProcID(*blob);
            return shrink;
         }
      }
      // 4. reflected point is worse than the second worst, but better than highest
      else if (compare2(blob->yfR, blob->ufR, blob->yfSecHigh, blob->ufSecHigh, mult) == 1
               && compare2(blob->yfR, blob->ufR, blob->yfHigh, blob->ufHigh, mult) == 2) {
        // 2nd type of contraction??? see drawing
        // replace the highest point with reflected point
         if (blob->aux_holds != contraction) {
            GlobalPrintDebug("Half-step in contraction", 5);
            // swap highest and R for breif moment while contraction point is calculated
            std::vector<point> temp;
            temp.resize(blob->fSize);
            
            for(arma::uword i = 0; i < blob->fSize; i++){
                  temp.at(i) = blob->vertex.at(i + blob->n_vert - blob->fSize);
                  blob->vertex.at(i + blob->n_vert - blob->fSize) = blob->auxR.at(i);   
                  blob->auxR.at(i) = temp.at(i);
               }
            std::cerr << __FILE__ << ' ' << __LINE__ << " before CheckProcID" << std::endl;
            CheckProcID(*blob);

            // Now find contracted point
            contract(blob);

            for(arma::uword i = 0; i < blob->fSize; i++){
               target.at(i) = blob->aux2.at(i).ProcID;
            }
   
            std::cerr << __FILE__ << ' ' << __LINE__ << " before CheckProcID" << std::endl;
            CheckProcID(*blob);

            // swap back to original points so that simplex can find 
            // way back here on next cycle-Will swap for last time if new cycle comes back here
            // and attempt to replace/compare contracted point
            for(arma::uword i = 0; i < blob->fSize; i++){
//               blob->vertex.at(i + blob->n_vert - blob->fSize) = temp.at(i);
               temp.at(i) = blob->auxR.at(i);
               blob->auxR.at(i) = blob->vertex.at(i + blob->n_vert - blob->fSize);   
               blob->vertex.at(i + blob->n_vert - blob->fSize) = temp.at(i);
               std::cerr << __FILE__ << ' ' << __LINE__ << " before CheckProcID" << std::endl;
               CheckProcID(*blob);
            
               if (undef == blob->aux2.at(i).status || done == blob->aux2.at(i).status) {
                  blob->aux2.at(i).Pstat = not_running;
                  std::stringstream messString;
                  messString << __FILE__ << ' ' << __LINE__ << " initiated Contraction";
                  GlobalPrintDebug(messString.str(), 5);
                  return initiate;
               }
               else {
                  // simulation is already running at this vertex
                  blob->aux2.at(i).Pstat = running;
                  std::stringstream messString;
                  messString << __FILE__ << ' ' << __LINE__ << " REinitiated Contraction";
                  GlobalPrintDebug(messString.str(), 5);
                  return restart;
               }
            }
         
         }

         // at this point the contraction should be valid, and running
         for(arma::uword i = 0; i < blob->fSize; i++){
            if (undef == blob->aux2.at(i).status) {
               std::cerr << "optimize: why is contraction not running?" << std::endl;
//            exit(1);
            }
         // but we may have to wait for a value
            if (pending == blob->aux2.at(i).status) {
               const time_t Ctime = std::time(0);
               std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                         << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
               return simplex_wait;
               }
               // but we may have to wait for a VALID value
               if (blob->aux2.at(i).value < 0) {
               const time_t Ctime = std::time(0);
               std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                         << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
               return simplex_wait;
            }
         }

         int OOI = compare2(blob->yf2, blob->uf2, blob->yfR, blob->ufR, mult);
         if (0 == OOI) {
            // error bar is overlapping
            const time_t Ctime = std::time(0);
            std::cout <<"FRAME  "<< __FILE__ << ' ' << __LINE__ << " simplex wait "
                      << std::asctime( std::localtime(&Ctime) ) << ' ' << std::endl;
            std::stringstream whatever;
            whatever << blob->yfSecHigh << " +/- "  << blob->ufSecHigh
                     << " < " << blob->yf2 << " +/- " << blob->uf2 << " < " 
                     <<  blob->yfR << " +/- "<< blob->ufR
                     << " \nContraction2? branch"<< std::endl;
            GlobalPrintDebug(whatever.str(), 5);
            return simplex_wait;
         }

         // Put this here, so that proper comparison can happen. If this is not here,
         // then reflection will infinitely repeat (theoretically)
         // NOW replacing old highest with NEW highest--old H > new H
         // Also, the error bar comparison already happened on OOI
         swapPoint(*blob, 0);
         int OOT = compare2(blob->yf2, blob->uf2, blob->yfHigh, blob->ufHigh, mult);
         // NOW old high is replaced by aux[0]--VALUE OF aux[0] is the NEW high

         // at this point we know the value at the contraction
         // Also, COMPARING "NEW High" point with others
         if (2 == OOT) {
            // keep contracted C < NEW H
            GlobalPrintDebug("keeping contracted point 2", 5);
            SIMPswap = true;  // tells calling fxn that simplex has moved
            swapPoint(*blob, 1);
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
