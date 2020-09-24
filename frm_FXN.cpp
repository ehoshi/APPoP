#include "frame.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <signal.h>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

#include <armadillo>


// This function is function to print debug statements.
// Dlevel is set and passed from the calling function, and is compared
//       to SetL value to print or not print
// Dlevel = 0 -always print
// Dlevel = 9 -only print on HardCore debug
//
// SetL =0 -Do NOT print ANY debug statment, except for "always print" statements
// SetL =1 -print least amount
// SetL =9 -HardCore debug. Print everything-this will clutter ALL of the log files!
void
PrintDebug(const std::string &StrToPrint, std::ofstream &outfile, const int Dlevel)
{
//   const time_t Ctime = std::time(0);
   int SetL = 5;
   if (Dlevel <= SetL) {
//      time_t jikan = std::time(NULL);
      outfile << StrToPrint << ' ' << std::endl;
//      outfile << std::asctime( std::localtime(&Ctime) ) << std::endl;
   }
}

//TODO commments
void
GlobalPrintDebug(std::string StrToPrint, const int Dlevel)
{
   PrintDebug(StrToPrint, GLOBAL_debuglog, Dlevel);
}

void
printsimplex(const simplex &blob, std::ostream& SSt)
{
   const time_t Ctime = std::time(0);
   SSt << "\n===========Current g, g_UNC values, and status===========" << std::endl;
   for (int j = 0; j < blob.vertices; j++) {
      SSt << "vertex" << j << " on proc " << blob.vertex[j].ProcID
          << " = " << blob.vertex[j].value << " +/- "
          << blob.vertex[j].error << " : Status = "
          << blob.vertex[j].status << std::endl;
   }
   for (int k = 0; k < MAXAUX; k++) {
      SSt << "aux" << k << " on proc " << blob.aux[k].ProcID
          << " = " << blob.aux[k].value << " +/- " << blob.aux[k].error
          << " : Status = " << blob.aux[k].status << std::endl;
   }

   SSt << "==================================================="
           << std::asctime( std::localtime(&Ctime) ) << '\n' << std::endl;
}

/*
   This function below contains the calculation functions used by workers to calculate
   objective function, g. NOTE that g and g(r) are two different values, as g is
   objective function and g(r) is pair-correlation function.
   return values for gcalc:
   0  = success
   1  = props.test exists, but it does not have any readable values
   2  = there is at least one NaN, inf or other limit type (useless) value in the results
   -1 = can NOT find experimental or simulation pair correlation function file(s)
   -2 = failed at cat | awk
*/

//TODO make a fuction to prepare input to gcalc.
/*
   arma::mat in_var;
   arma::mat unc_var;
   std::ostringstream Vin;
   std::ostringstream Uin;

   Vin << "var.txt";
   Uin << "unc.txt";

//preparing the input files
   int C = system("cat props.prod | awk '/^<(PE|P|D)>/{print $3}' > var.txt");
   if (C == -1 || WEXITSTATUS(C) !=0)
   {
      outfile << "error in conversion to gcalc input--C" << std::endl;
      WStatus = error;
      return -2;
   }

   PrintDebug("test message between cat-inside gcalc", outfile, 9);

   int D = system("cat props.prod | awk '/^<(PE|P|D)>/{print $5}' > unc.txt");
   if (D == -1 || WEXITSTATUS(D) !=0)
   {
            outfile << "error in conversion to gcalc input--D" << std::endl;
            WStatus = error;
            return -2;
   }

   PrintDebug("PRINT AFTER cat inside gcalc", outfile, 9);

//check for nan, inf, or any other limit type value in the file
   if( is_file_useful(Vin) && is_file_useful(Uin) ){
         if( !in_var.load("var.txt") || !unc_var.load("unc.txt") ){
            PrintDebug("props.prod exits but does not have readable value\n", outfile, 8);
            WStatus = pending;
            return 1;
         }
   PrintDebug("PRINT INSIDE is_file_useful-gcalc", outfile, 9);
   }
   else{
      PrintDebug("At least one NaN value exitsts\n", outfile, 8);
      WStatus = pending;
      return 2;
   }

*/
int
gcalc(double &grFinal, double &guFinal, WORKER_status &WStatus,
      std::ofstream &outfile, unsigned long int cnt, double prev_g,
      double prev_u)
{
   PrintDebug("This is test message inside gcalc", outfile, 9);

   arma::mat in_var;
   arma::mat unc_var;
   std::ostringstream Vin;
   std::ostringstream Uin;

   Vin << "var.txt";
   Uin << "unc.txt";

//preparing the input files (remove for generalized version)
/*
   int C = system("cat props.prod | awk '/^<(PE|P|D)>/{print $3}' > var.txt");
   if (C == -1 || WEXITSTATUS(C) !=0)
   {
      outfile << "error in conversion to gcalc input--C" << std::endl;
      WStatus = error;
      return -2;
   }

   PrintDebug("test message between cat-inside gcalc", outfile, 9);

   int D = system("cat props.prod | awk '/^<(PE|P|D)>/{print $5}' > unc.txt");
   if (D == -1 || WEXITSTATUS(D) !=0)
   {
            outfile << "error in conversion to gcalc input--D" << std::endl;
            WStatus = error;
            return -2;
   }

   PrintDebug("PRINT AFTER cat inside gcalc", outfile, 9);
*/

   // check for nan, inf, or any other limit type value in the file
   if (is_file_useful(Vin) && is_file_useful(Uin)) {
      if ( !in_var.load("var.txt") || !unc_var.load("unc.txt") ) {
         PrintDebug("var & unc does not have readable value\n", outfile, 8);
         WStatus = pending;
         return 1;
      }
      PrintDebug("PRINT INSIDE is_file_useful-gcalc", outfile, 9);
   }
   else {
      PrintDebug("At least one NaN value exitsts\n", outfile, 8);
      WStatus = pending;
      return 2;
   }

   {
      std::ostringstream iVstr;
      std::ostringstream iUstr;

      iVstr << "PRINTED FROM gcalc in_var = \n"  << in_var << std::endl;
      iUstr << "PRINTED FROM gcalc unc_var = \n" << unc_var << std::endl;

      PrintDebug(iVstr.str(), outfile, 7);
      PrintDebug(iUstr.str(), outfile, 7);
   }
/*For generalized version
  Read a file (LitV.in) that contains
  exp values on col1
  unc values on col2
  A   values on col3
  to one matrix

  arma::mat LitV
  arma::vec Exp(LitV.n_rows)
  arma::vec Eunc(LitV.n_rows)
  arma::vec A(LitV.n_rows)

  if( !LitV.load("LitV.in") )
  {
      return -1;
  }

  for(amra::uword z = 0; z < LitV.n_col; z++){
      for(arma::uword y = 0; y < LitV.n_row;y++){
        if(z == 0){
           Exp(y) = LitV(y,z);
        }
        if(z == 1){
           Eunc(y) = LitV(y,z);

        }
        if(z == 2){
           A(y) = LitV(y,z);

        }
      }
  }


*/


   // variables, vectors, and matricies used for NGOR-not paircorelation fxn
   arma::vec Exp(3);  // vector to store the experimental value
   arma::vec Eunc(3);  // vector to store experiental uncertainty
   arma::vec A(3);  // vector to store the weighing factor

   Exp(0) = -9.918;  // dH(vap) =  -43.98 kJ/mol DIVIDE THIS BY 4.184 TO GET KCAL/MOL
                     // then subtract 0.592 kcal/mol to get internal energy
   Exp(1) = 1.0; // Pressure atm
   Exp(2) = 2.299e-5;  // D-self diffusion coeffiecent cm^2 / s

   Eunc(0) = 0.01164;
   Eunc(1) = 17.6638;
   Eunc(2) = Exp(2) * 0.001;

   // weighing factor can be thought as tolerance
   // calculated as 1/A
   A(0) = 1E0;
   A(1) = 5E2;
   A(2) = 4e-6;

   // temp value to store the obj.fxn value for non-pair correlation function
   double gvalue;
   double uncertainty;

   // calculating the uncertainty for the above:
   // 1. subtraction
   arma::mat temp_mat_VAL1 = in_var - Exp;
   arma::mat temp_mat_unc1 = arma::sqrt(unc_var % unc_var + Eunc % Eunc);

   // 2. square
   // multiply unc (relative) by 2 to account for accumulation of error
   // for taking a square of a value
   arma::mat temp_mat_unc2 = 2 * temp_mat_unc1 / temp_mat_VAL1;  // relative error
   arma::mat temp_mat_VAL2 = temp_mat_VAL1 % temp_mat_VAL1;

   // 3. mulitply by weighing factor squared
   // error is skipped since relative error stays same, A_error = 0
   arma::mat temp_mat_VAL3 = temp_mat_VAL2 / (A % A);


   // NEW 4. sum
   double temp_VAL4 = arma::accu( temp_mat_VAL3 );
   arma::mat temp_mat_unc3 = temp_mat_unc2 % temp_mat_VAL3;  //UNDO rel error
   double temp_unc4 = std::sqrt( arma::accu(temp_mat_unc3 % temp_mat_unc3) );

   // new 5. square root
   //divide relative error by 2 to get the uncertainty.
   //Then undo at the end
   gvalue = std::sqrt( temp_VAL4 );
   double temp_unc5 = temp_unc4 / temp_VAL4 / 2; //relative error / 2
   uncertainty = temp_unc5 * gvalue; //undo rel error
   /*
     variable to store the obj.fxn value for g(r) and associated uncertainties:
     gr = obj.fxn value for g(r)
     gu = uncertainty associated with g(r)
     label with 1 => O-O pair correlation function
     label with 2 => O-H pair correlation function
     label with 3 => H-H pair correlation function
   */
   double gr1 = 0;
   double gr2 = 0;
   double gr3 = 0;
   double gu1 = 0;
   double gu2 = 0;
   double gu3 = 0;

   // these are string to assign file names,which are files generated by simulation
   // Used for reading g(r) data from file
   std::string eFile1    = "gOO_Soper.dat";
   std::string eFile2    = "gOH_Soper.dat";
   std::string eFile3    = "gHH_Soper.dat";
   std::string sFile1    = "g0808.prod";
   std::string sFile2    = "g0801.prod";
   std::string sFile3    = "g0101.prod";

   // OO piar correlation function calculation
   int gr_calcChk = gr_calc(eFile1, sFile1, gr1, gu1);
   if (-1 == gr_calcChk) {
      outfile << "ERROR: cannot load one or more of the pair-correlation function"
              << " files.\n This can be caued by simulaion program failed to produce "
              << "pair correlation function file\n or experimental results are missing"
              << " from tempdir" << std::endl;

      WStatus = error;
      return -1;
   }

   // OH pair correlation function calculation
   gr_calcChk = gr_calc(eFile2, sFile2, gr2, gu2);
   if (-1 == gr_calcChk) {
      outfile << "ERROR: cannot load one or more of the pair-correlation function"
              << " files.\n This can be caued by simulaion program failed to produce "
              << "pair correlation function file\n or experimental results are missing"
              << " from tempdir" << std::endl;

      WStatus = error;
      return -1;
   }

   // HH pair correlation function calculation
   gr_calcChk = gr_calc(eFile3, sFile3, gr3, gu3);
   if (-1 == gr_calcChk) {
      outfile << "ERROR: cannot load one or more of the pair-correlation function"
              << " files.\n This can be caued by simulaion program failed to produce "
              << "pair correlation function file\n or experimental results are missing"
              << " from tempdir" << std::endl;

      WStatus = error;
      return -1;
   }

   // calculate the objfxn for g(r)
   //TODO find somewhere to store the A, weighing values for gr
   grFinal =  1.0*gr1 + 1.0*gr2 + 1.0*gr3 + gvalue;

   // now calculate the final value for uncertainty.
   double temp1 = gu1 * gu1 * 1.0 * 1.0;
   double temp2 = gu2 * gu2 * 1.00 * 1.00;
   double temp3 = gu3 * gu3 * 1.0 * 1.0;
   double temp4 = uncertainty * uncertainty;

   guFinal = std::sqrt(temp1 + temp2 +temp3 + temp4);

   if ( (prev_g > 0 && prev_u > 0)
        && (std::abs(prev_g - grFinal) >= prev_g * 0.01
            || std::abs(prev_u - guFinal) >= prev_u * 0.01) ) {

      const time_t Ctime = std::time(0);

      std::ostringstream terms;
//      terms << "# LOOP grF guF gvalue unc gr1 gr2 gr3 gu1 gu2 gu3 E P D Eu Pu Du time" << std::endl;

      terms << std::setw(8) << cnt << ' '
            << std::setw(12) << grFinal << ' '
            << std::setw(12) << guFinal << ' '
            << std::setw(12) << gvalue << ' '
            << std::setw(12) << uncertainty << ' '
            << std::setw(12) << gr1 << ' '
            << std::setw(12) << gu1 << ' '
            << std::setw(12) << gr2 << ' '
            << std::setw(12) << gu2 << ' '
            << std::setw(12) << gr3 << ' '
            << std::setw(12) << gu3 << ' ';
      for (arma::uword i = 0; i < in_var.n_rows; i++) {
         terms << std::setw(8) << in_var(i) << ' '
               << std::setw(8) << unc_var(i) << ' ';
      }
      terms << std::asctime( std::localtime(&Ctime) );
      PrintDebug(terms.str(), outfile, 1);
   }

   // returns 0 with gcalc success
   WStatus = active;
   return 0;
}


int
gcalc2(double &grFinal, double &guFinal, WORKER_status &WStatus,
      std::ofstream &outfile, unsigned long int cnt, double prev_g,
      double prev_u)
{
   PrintDebug("This is test message inside gcalc", outfile, 9);

   arma::mat in_var;
   arma::mat unc_var;
   std::ostringstream Vin;
   std::ostringstream Uin;

   Vin << "var.txt";
   Uin << "unc.txt";

   // check for nan, inf, or any other limit type value in the file
   if (is_file_useful(Vin) && is_file_useful(Uin)) {
      if ( !in_var.load("var.txt") || !unc_var.load("unc.txt") ) {
         PrintDebug("var & unc does not have readable value\n", outfile, 8);
         WStatus = pending;
         return 1;
      }
      PrintDebug("PRINT INSIDE is_file_useful-gcalc", outfile, 9);
   }
   else {
      PrintDebug("At least one NaN value exitsts\n", outfile, 8);
      WStatus = pending;
      return 2;
   }

   {
      std::ostringstream iVstr;
      std::ostringstream iUstr;

      iVstr << "PRINTED FROM gcalc in_var = \n"  << in_var << std::endl;
      iUstr << "PRINTED FROM gcalc unc_var = \n" << unc_var << std::endl;

      PrintDebug(iVstr.str(), outfile, 7);
      PrintDebug(iUstr.str(), outfile, 7);
   }
/*For generalized version
  Read a file (LitV.in) that contains
  exp values on col1
  unc values on col2
  A   values on col3
  to one matrix

  arma::mat LitV
  arma::vec Exp(LitV.n_rows)
  arma::vec Eunc(LitV.n_rows)
  arma::vec A(LitV.n_rows)

  if( !LitV.load("LitV.in") )
  {
      return -1;
  }

  for(amra::uword z = 0; z < LitV.n_col; z++){
      for(arma::uword y = 0; y < LitV.n_row;y++){
        if(z == 0){
           Exp(y) = LitV(y,z);
        }
        if(z == 1){
           Eunc(y) = LitV(y,z);

        }
        if(z == 2){
           A(y) = LitV(y,z);

        }
      }
  }


*/


   // variables, vectors, and matricies used for NGOR-not paircorelation fxn
   arma::vec Exp(3);  // vector to store the experimental value
   arma::vec Eunc(3);  // vector to store experiental uncertainty
   arma::vec A(3);  // vector to store the weighing factor

   Exp(0) = -9.918;  // dH(vap) =  -43.98 kJ/mol DIVIDE THIS BY 4.184 TO GET KCAL/MOL
                     // then subtract 0.592 kcal/mol to get internal energy
   Exp(1) = 1.0; // Pressure atm
   Exp(2) = 2.299e-5;  // D-self diffusion coeffiecent cm^2 / s

   Eunc(0) = 0.01164;
   Eunc(1) = 17.6638;
   Eunc(2) = Exp(2) * 0.001;

   //Tolerance ratio. The original constant vector
   //made to relative error to match val2
//original values. Unitless
/*
   A(0) = 1E-3 / -9.918;
   A(1) = 3.5E2 / 1.0;
   A(2) = 5e-7 / 2.299E-5;         
 make the values to be 1/A
   A(0) = 1/1E-3;
   A(1) = 1/3.5E2;
   A(2) = 1/5e-7;         
*/
   //testing tolerance ratio. It is 1/A
   //Units are also 1/A (relative version)
   A(0) = 1/(1E-3 / -9.918);
   A(1) = 1/(3.5E2 / 1.0);
   A(2) = 1/(5e-7 / 2.299E-5);         

   // temp value to store the obj.fxn value for non-pair correlation function
   double gvalue;
   double uncertainty;

//   arma::mat tempC1 = Eunc / Exp;//relative error of the A (unitless A values)

   //normalize A. Does not have any unc.

   arma::mat NormC = normalise(A);
 
//changed the objective function
   //1. subtraction-keep the same
   //XXX VAL1: change this, if needed. The difference is divided by a constant 
   //(same value as experimental value and units) with no error to get rid of units
//   arma::mat temp_mat_VAL1 = in_var - Exp;
   arma::mat temp_mat_VAL1 = (in_var - Exp) / Exp;
   arma::mat temp_mat_unc1 = arma::sqrt(unc_var % unc_var + Eunc % Eunc);
   //XXX original unc2
//   arma::mat temp_mat_unc2 = temp_mat_unc1 / temp_mat_VAL1;//relative error of the unc1. Original

//Since division is considered constant with no error, relative error is calculated from subtraction only
   arma::mat temp_mat_unc2 = temp_mat_unc1 / (in_var - Exp);//relative error of the unc1, subtraction only


   //2. change val1 to rel. omitted due to large propagation of error
//   arma::mat temp_mat_VAL2 = temp_mat_VAL1 / Exp;
//   arma::mat temp_mat_unc2 = temp_mat_unc1 / temp_mat_VAL1;//relative error of the unc1

   //relative error after division (rel.error of val2)
//   arma::mat temp_mat_unc3 = arma::sqrt( temp_mat_unc2 % temp_mat_unc2 + Eunc % Eunc);

   //3. determine parallel part. val2 is the P vector
   double parallel_val = 0.0;
   double parallel_unc = 0.0;

   //take dot product of P vector and NormC (w/uncertainty)
   //Calculates A_parallel. Gives abs.unc.
   double tempDunc = 0.0;
   for(arma::uword i = 0; i < NormC.n_rows; i++){
      double test_val = temp_mat_VAL1(i) * NormC(i);
      //rel.error for the mult.step of the dot product
      parallel_val += test_val; 
      //convert rel.error to abs.error, 
      //square, and sum to calculate the error of the summation of the dot product
      tempDunc += (temp_mat_unc2(i) * test_val) * (temp_mat_unc2(i) * test_val)  ;

   }

   parallel_unc = std::sqrt(tempDunc);

   //4. find the projection of the parallel vector
   double perp_unc = 0.0;
   double perp_val = arma::norm( temp_mat_VAL1 - (parallel_val*NormC) );

   //get the rel.unc of parallel 
   double tempUnc1 = parallel_unc/parallel_val;
   
   //This is to match the size of temps w/NormC
   arma::mat tempUnc2 = NormC;
   arma::mat tempVAL1 = NormC;//this is used for calculating errors in normalization
   //get abs.unc of the multiplication (tempUnc1 is rel.unc)
   for(arma::uword i = 0; i < NormC.n_rows; i++){
      double tempC = parallel_val * NormC(i);
      tempUnc2(i) = tempC * tempUnc1;
      tempVAL1(i) = temp_mat_VAL1(i) - tempC;
   }
    
    
   //abs.unc of VAL1
   arma::mat temp_mat_unc4 = temp_mat_unc2 % temp_mat_VAL1;
   arma::mat temp_mat_unc5 = NormC;
   
   //get error for the subtraction (abs.unc)
   for(arma::uword i = 0; i < NormC.n_rows; i++){
   temp_mat_unc5(i) = std::sqrt( tempUnc2(i) * tempUnc2(i) + temp_mat_unc4(i) * temp_mat_unc4(i) );
   }

   //find error for calculating the magnitude
   arma::mat tempUnc3 = temp_mat_unc5 / tempVAL1; //rel.error of VAL2 (subtraction)
   //square of each term (rel.unc) then convert to abs.unc
//   arma::mat tempUnc4 = tempVAL1 % tempVAL1 % arma::sqrt(tempUnc3 % tempUnc3 + tempUnc3 % tempUnc3);   
   arma::mat tempUnc4 = tempVAL1 % tempVAL1 % tempUnc3 * 2;//uses power rule
   //sum (abs.unc) then take the square root of the sum unc. Gives rel.unc
   double tempUnc5 = std::sqrt( arma::accu( tempUnc4 % tempUnc4 ) ) /arma::accu(tempVAL1%tempVAL1) / 2;
   //get abs.unc of above
   perp_unc = tempUnc5 * perp_val;

   gvalue = (parallel_val * parallel_val) + (100 * perp_val * perp_val);

   //calculate the unc for gvalue
   double temparaUNC1 = parallel_unc / parallel_val;
   double temparaUNC2 = 0;
   double temparaUNC3 = 0;
   double temperpUNC1 = 0;
   double temperpUNC2 = 0;

   //unc of multiplication of the gvalue (square). Rel.unc
//   temparaUNC2 = std::sqrt( temparaUNC1 * temparaUNC1 + temparaUNC1 * temparaUNC1 ) ;
   temparaUNC2 = 2 * temparaUNC1;//again, power rule
//   temperpUNC1 = std::sqrt(tempUnc5 * tempUnc5 + tempUnc5 * tempUnc5) ;
   temperpUNC1 = 2 * tempUnc5;//same

   //comvert above error to abs. error
   temparaUNC3 = temparaUNC2 * parallel_val * parallel_val;
   temperpUNC2 = temperpUNC1 * perp_val * perp_val * 100;
   //error of sum
   
   uncertainty = std::sqrt(temparaUNC3 * temparaUNC3 + temperpUNC2 * temperpUNC2);

/*
     variable to store the obj.fxn value for g(r) and associated uncertainties:
     gr = obj.fxn value for g(r)
     gu = uncertainty associated with g(r)
     label with 1 => O-O pair correlation function
     label with 2 => O-H pair correlation function
     label with 3 => H-H pair correlation function
   */
   double gr1 = 0;
   double gr2 = 0;
   double gr3 = 0;
   double gu1 = 0;
   double gu2 = 0;
   double gu3 = 0;

   // these are string to assign file names,which are files generated by simulation
   // Used for reading g(r) data from file
   std::string eFile1    = "gOO_Soper.dat";
   std::string eFile2    = "gOH_Soper.dat";
   std::string eFile3    = "gHH_Soper.dat";
   std::string sFile1    = "g0808.prod";
   std::string sFile2    = "g0801.prod";
   std::string sFile3    = "g0101.prod";

   // OO piar correlation function calculation
   int gr_calcChk = gr_calc(eFile1, sFile1, gr1, gu1);
   if (-1 == gr_calcChk) {
      outfile << "ERROR: cannot load one or more of the pair-correlation function"
              << " files.\n This can be caued by simulaion program failed to produce "
              << "pair correlation function file\n or experimental results are missing"
              << " from tempdir" << std::endl;

      WStatus = error;
      return -1;
   }

   // OH pair correlation function calculation
   gr_calcChk = gr_calc(eFile2, sFile2, gr2, gu2);
   if (-1 == gr_calcChk) {
      outfile << "ERROR: cannot load one or more of the pair-correlation function"
              << " files.\n This can be caued by simulaion program failed to produce "
              << "pair correlation function file\n or experimental results are missing"
              << " from tempdir" << std::endl;

      WStatus = error;
      return -1;
   }

   // HH pair correlation function calculation
   gr_calcChk = gr_calc(eFile3, sFile3, gr3, gu3);
   if (-1 == gr_calcChk) {
      outfile << "ERROR: cannot load one or more of the pair-correlation function"
              << " files.\n This can be caued by simulaion program failed to produce "
              << "pair correlation function file\n or experimental results are missing"
              << " from tempdir" << std::endl;

      WStatus = error;
      return -1;
   }

   // calculate the objfxn for g(r)
   //TODO find somewhere to store the A, weighing values for gr
   grFinal =  1.0*gr1 + 1.0*gr2 + 1.0*gr3 + gvalue;

   // now calculate the final value for uncertainty.
   double temp1 = gu1 * gu1 * 1.0 * 1.0;
   double temp2 = gu2 * gu2 * 1.00 * 1.00;
   double temp3 = gu3 * gu3 * 1.0 * 1.0;
   double temp4 = uncertainty * uncertainty;

   guFinal = std::sqrt(temp1 + temp2 +temp3 + temp4);

   if ( (prev_g > 0 && prev_u > 0)
        && (std::abs(prev_g - grFinal) >= prev_g * 0.01
            || std::abs(prev_u - guFinal) >= prev_u * 0.01) ) {

      const time_t Ctime = std::time(0);

      std::ostringstream terms;
//      terms << "# LOOP grF guF gvalue unc gr1 gr2 gr3 gu1 gu2 gu3 E P D Eu Pu Du time" << std::endl;

      terms << std::setw(8) << cnt << ' '
            << std::setw(12) << grFinal << ' '
            << std::setw(12) << guFinal << ' '
            << std::setw(12) << gvalue << ' '
            << std::setw(12) << uncertainty << ' '
            << std::setw(12) << parallel_val << ' '
            << std::setw(12) << parallel_unc << ' '
            << std::setw(12) << perp_val << ' '
            << std::setw(12) << perp_unc << ' '
            << std::setw(12) << gr1 << ' '
            << std::setw(12) << gu1 << ' '
            << std::setw(12) << gr2 << ' '
            << std::setw(12) << gu2 << ' '
            << std::setw(12) << gr3 << ' '
            << std::setw(12) << gu3 << ' ';
      for (arma::uword i = 0; i < in_var.n_rows; i++) {
         terms << std::setw(8) << in_var(i) << ' '
               << std::setw(8) << unc_var(i) << ' ';
      }
      terms << std::asctime( std::localtime(&Ctime) );
      PrintDebug(terms.str(), outfile, 1);
   }

   // returns 0 with gcalc success
   WStatus = active;

   return 0;
}
////////////////////////////////////////////////////////////////////////
/*
   gr_calc caluclates the uncertainty and obj.fxn value of pair correlation funciton.
   Here are the return values and meaning:
   -1 = unable to find or read the files with experimental data
    0 = normal execution
*/
int
gr_calc(std::string eStr, std::string sStr, double &gr, double &gu)
{
   double erMax;  // laergest value of r found in experimental data
   double erMin;  // smallest value of r found in experimental data
   double SoS = 0;  // sum of stuff-temporary place to store the calculated value
   arma::uword i;  // counter for loop
   arma::uword j;  // counter for loop

   /*
     variables, strings, vectors, and matricies used for GOR-paircorrelation fxn
     e -> experimental value
     s -> simulation value
   */
   arma::mat e;  // experimental matrix for pair correlation fxn
   arma::mat s;  // simulation matrix for pair correlation fxn

   // load the input files
   if ( !e.load(eStr) ) {
      return -1;
   }

   if ( !s.load(sStr) ) {
      return -1;
   }

   /*
     convert input matrix to a vector of:
     col 0 = x-axis or r
     col 1 = y-axix or g(r)
     col 2 = uncertainty associated with g(r)
   */
   arma::vec er    = e.col(0);  // r values for experimental data
   arma::vec eGofr = e.col(1);  // g(r) values for experimental data
   arma::vec eEor  = e.col(2);  // uncertainty in g(r) for experimental data
   arma::vec sr    = s.col(0);  // r values in simulation data
   arma::vec sGofr = s.col(1);  // g(r) values for simulaiton
   arma::vec sEor  = s.col(2);  // uncertainty associated with g(r) simulaiton

   /*
     vector to store only the values within the range of exiperimental x-axis or r value
     stores the "chopped" simulation outputs
   */
   arma::vec sr_Chop;
   arma::vec sGofr_Chop;
   arma::vec sEor_Chop;

   /*
     the next few lines perform the following tasks:
     set the erMax to maximum value of r found in the experimental value
     set the erMin to minimum value of r found in the experimental value
     "chop" original sr so that all sr value is within the er
     do similar for sGofr and sEor
     this is to ensure that all the simulation values fit within the range of experimental values.
     IMPORTANT NOTE: this whole steps rely on all the e files sorted from lowest to highest r
   */

   erMax = er(er.n_rows - 1);
   erMin = er(0);
   double delta_E = erMax - erMin;

   sr_Chop = sr.elem(find (sr <= erMax && sr >= erMin) );
   sGofr_Chop = sGofr.elem(find (sr <= erMax && sr >= erMin) );
   sEor_Chop  = sEor.elem( find (sr <= erMax && sr >= erMin) );

   arma::vec eInterp(sr_Chop.n_rows);  // stores interpolated values
   arma::vec eorInterp(sr_Chop.n_rows);  // stores error accumulation
   arma::vec delta_x(sr_Chop.n_rows);  // stores values of dx between two r.
   arma::vec diffSqMult_error(sr_Chop.n_rows);  // error for subtraction, 
                                                // square, and mult in objfxn.

   // interpolate using "chopped off" simulaion results and experimental value.  

   for (i = 0; i < sr_Chop.n_rows; i++) {
      //this j loop is for interpolation and piece of error calculation.
      for (j = 0; j < er.n_rows - 1; j++) {
         if (sr_Chop(i) > er(j) && sr_Chop(i) < er(j+1)) {
            // interpolate. This equation is special case of general Lagrange interpolation.
            // From Numerical Recepies P.86
            // OR fancy way of saying linear interpolation
            double A = ( er(j+1) - sr_Chop(i) ) / ( er(j+1) - er(j) );
            double B = 1 - A;

            eInterp(i) = A * eGofr(j) + B * eGofr(j+1);

            // error accumilation
            // these values are temporary values for calculating random/statistical error.
            double eA;
            double eB;

            // A and B has error accumilation of 0
            // eA = abs.error of operation A * eGofr(j) 
            // eB = abs.error of operation B * eGofr(j+1) 
            eA = A * eEor(j);
            eB=  B * eEor(j+1);

            // storeing the interp.error values into eorInterp vector.
            eorInterp(i) = std::sqrt( (eA * eA) + (eB * eB) );

            break;
         }
         else if (sr_Chop(i) == er(j)) {
            eInterp(i) = eGofr(j);
            eorInterp(i) = eEor(j);
            break;
         }
         else if (sr_Chop(i) == er(er.n_rows - 1)) {
            eInterp(i) = eGofr(er.n_rows - 1);
            eorInterp(i) = eEor(er.n_rows - 1);
            break;
         }
      }  //end j loop

      // these statment ensures to avoid divide by 0 in the steps to come
      if (i != sr_Chop.n_rows-1) {
         delta_x(i) = sr_Chop(i+1) - sr_Chop(i);
      }
      else {
         delta_x(i) = 0.0;
      }

      // calculates piece of obj_fxn. Used again outside of this loop.
      double diff;
      double subt_square;
      double subt_X_deltaX;
      diff  = sGofr_Chop(i) - eInterp(i);
      subt_square = diff * diff;
      subt_X_deltaX = ( subt_square * delta_x(i) ) ;
      SoS = SoS + subt_X_deltaX;

      // calculates piece of propagated error. Calc up to multiplying by delta x.
      // if statement is there to avoid divide by 0 on the next step.
      if (diff != 0.0) {
         //uncertanty for subtraction, square, and mult by delta_x
         double temp_sq1;
         double temp_sq2;
         double sqrt_sos;  // square root of sum of square
         double tpcnt_y;  // temp percent for y
         double rel_pcnt_error;
         double abs_error_y;

         // error for subtraction
         temp_sq1 = sEor_Chop(i) * sEor_Chop(i);
         temp_sq2 = eorInterp(i) * eorInterp(i);
         sqrt_sos = std::sqrt(temp_sq1 + temp_sq2);

         // get relative error for square error calc and calculate error of square
         tpcnt_y = sqrt_sos / std::abs(diff);
         rel_pcnt_error = 2 * tpcnt_y;

         // change previous rel_pcnt to abs_pcnt, (after square of diff).
         // delta_x is included here since it has error of 0.
         abs_error_y = rel_pcnt_error * (diff * diff * delta_x(i) );
         diffSqMult_error(i) = abs_error_y;
      }
      else {
         diffSqMult_error(i) = 0.0;
      }
   }  // end i loop

   // calculate obj fxn value
   gr = std::sqrt(SoS / delta_E);

   // error accumilation for objfxn
   // NOTE: delta_E has error of 0.
   double temp_sum_e1;
   double temp_sum_e2;
   double temp_sum_e3;
   temp_sum_e1 = std::sqrt( arma::accu (diffSqMult_error % diffSqMult_error) );  // summation
   temp_sum_e2 = temp_sum_e1 / SoS;  // gives rel.error from previous line
   temp_sum_e3 = temp_sum_e2 / 2;  // rel.error for taking sqrt
   gu  = temp_sum_e3 * gr;

   return 0;
}

////////////////////////////////////////////////////////////////////////
/*
  Here on down is auxillariy functions that is used for multiple reasons.
  1. convinence
  2. not to cram the main code
  3. because it is needed
*/
////////////////////////////////////////////////////////////////////////

// --CALLED BY WORKER ONLY
// FXN to check for NaN
bool
has_nan(const arma::mat m)
{
   for(arma::uword i = 0; i < m.n_elem; i++) {
      if (std::isnan(m[i])) {
         return true;
      }
   }

   return false;
}

// FXN to test for wether the file is useful. Armadillo does not
// read non-numbers (words/letters) other than NaN.
bool
is_file_useful(std::ostringstream &FileName)
{
   std::fstream DatFile;
   DatFile.open(FileName.str().c_str(), std::ios::in);
   double tempD = 0;
   int count = 0;
   while (DatFile >> tempD) {
      count++;
   }

   if ( !DatFile.eof() ) {
      DatFile.close();
      return false;
   }

   DatFile.close();
   return true;
}

// FXN to check for current child status for provided PID
int
StatChk(const pid_t pid, WORKER_status &Wstatus, const int Rank)
{
   std::ostringstream Dstring;

   pid_t ws;  // status checking variable
   int childExitStatus;

   ws = ::waitpid(pid, &childExitStatus, WNOHANG);
   if (-1 == ws && EINTR == errno) {
      std::cerr << "waitpid interrupted on proc " << Rank << '.' << std::endl;
   }

   if (WIFEXITED(childExitStatus) && ws > 0) {
      std::ostringstream Fstring;
      if (WEXITSTATUS(childExitStatus) == 0) {
         // normal termination of child process
         // This may or may not be an error dependin on the situation.
         // Overwrite of default Wstatus/Wstat happens outside of this fxn.
         Fstring << "Normal termination of the child process" << std::endl;
         GlobalPrintDebug(Fstring.str(), 7);

         Wstatus = finished;
         return 1;
      }
      else {
         // child process did not execute properly
         GlobalPrintDebug("ERROR: The child procss did not execute properly", 1);

         Wstatus = error;
         return 2;
      }
   }
   else if (WIFSIGNALED(childExitStatus) && ws > 0) {
      std::ostringstream astring;
      // external child killed with status of Ret
      astring << "ERROR: external process (child) killed with status of: "
              << WTERMSIG(childExitStatus) << std::endl;

      GlobalPrintDebug(astring.str(), 1);
      Wstatus = error;
      return 3;
   }
   else if (WIFSTOPPED(childExitStatus) && ws > 0) {
      std::ostringstream astring;
      // cchild process stopped with status of Ret
      astring << "Child process stopped with status of: "
              << WSTOPSIG(childExitStatus) << std::endl;
      GlobalPrintDebug(astring.str(), 1);

      Wstatus = stall;
      return 4;
   }
   else if (ws < 0) {
      //failed to call or run waitpid
      GlobalPrintDebug("ERROR: failed to call or run waitpid", 1);
      std::cerr << "ERROR: failed to call or run waitpid on proc " << Rank << std::endl;
      std::perror("waitpid");

      Wstatus = error;
      return 5;
   }

   // no abnormal/error as of now
   Dstring << "No error so far in external process -- parent" << std::endl;
   GlobalPrintDebug(Dstring.str(), 9);

   Wstatus = pending;
   return 0;
}

// FXN used to kill job of the provided PID
void
killJob(const pid_t cid)
{
   std::ostringstream Kstring;
   Kstring << "Testing from killJob fxn: the pid killed is: " << cid;
   GlobalPrintDebug(Kstring.str(), 5);
   int killReturn = ::kill(cid, SIGTERM);
   // wait here for process cid to kill???
   if (-1 == killReturn) {
      if (ESRCH == errno) {
         GlobalPrintDebug("pid does not exist. Proceed to new setup ", 1);
      }
      else if (EPERM == errno) {
         GlobalPrintDebug("no permission to send signal. Forcefully kill process", 1);
         // no permission to send signal. TODO forcfully kill process
         ::kill(cid, SIGKILL);  // I don't think I really want to do this...
      }
   }
}

// FXN used to print the worker situation.
void
PrintWsit(const WRK_situation Wsit, std::ofstream &outfile)
{
   if (0 == Wsit) {
      outfile << "Wsit is NoStart" << std::endl;
   }
   if (1 == Wsit) {
      outfile << "Wsit is EqStart" << std::endl;
   }
   if (2 == Wsit) {
      outfile << "Wsit is ProdStart" << std::endl;
   }
   if (3 == Wsit) {
      outfile << "Wsit is ProdAvail" << std::endl;
   }
   if (4 == Wsit) {
      outfile << "Wsit is ProdFin" << std::endl;
   }
}

// --ONLY CALLED BY MASTER.
// takes in Wstat and decides what Mcomm to give
// in-Wstat, return-Mcomm
// CPN = current process number that is looping over in MAS
MASTER_comand
MASwork(const WORKER_status WSTAT, const int CPN, const simplex blob)
{
   // really, there is no purpose to doing this, other than identifying the process
   // but this is also used for failsafe.
   bool IsCPN_aux = false;

   for (int i = 0; i < blob.vertices; i++) {
      if (blob.vertex[i].ProcID == CPN && blob.vertex[i].Pkind == auxil) {
         IsCPN_aux = true;
         break;
      }
   }

   for(int i = 0; i < MAXAUX; i++){
      if (blob.aux[i].ProcID == CPN && blob.aux[i].Pkind == auxil) {
         IsCPN_aux = true;
         break;
      }
   }

   if (undef == WSTAT && !IsCPN_aux) {
      std::ostringstream gnirts;
      gnirts << "Problem:Wstat=UNDEF in Proc " << CPN << ", which is NOT AUX\n."
             << "This should NOT be happeneing after initialization, but Mcomm is set to 'start'.";
      GlobalPrintDebug(gnirts.str(), 1);

      return start;
   }

   if (active == WSTAT) {
      std::ostringstream gnirts;
      gnirts << "Proc " << CPN << " is active. Mcomm is set to 'GMTV'.";
      GlobalPrintDebug(gnirts.str(), 7);

      return GMTV;
   }

   if (finished == WSTAT || stall == WSTAT) {
      std::ostringstream gnirts;
      gnirts << "Problem:FINISHED/STALL in Proc " << CPN << " is finished or stalled. "
             << "\nMcomm is set to 'reinit'-" << WSTAT;
      GlobalPrintDebug(gnirts.str(), 1);
      return reinitiate;
   }
 
   if (error == WSTAT) {
      std::ostringstream gnirts;
      gnirts << "ERROR: Proc " << CPN << " is has some error. Mcomm is set to 'reinit'.";
      GlobalPrintDebug(gnirts.str(), 1);

      return reinitiate;
   }

   if (infected == WSTAT) {
      std::ostringstream gnirts;
      gnirts << "Proc " << CPN << " is now zombified"; 
      GlobalPrintDebug(gnirts.str(), 9);

      return contagious;
   }
   // same as if WSTAT == pending
   return nothing;
}

// updateWstat
//
// inputs: 
// WST-worker status that was received in MASloop
// CPN-Current Process Number in MASloop
// PN - total number of pararell processes
// outputs(pass by refernce):
// blob-Wstat and ComEx of proc CPN is written
void
updateMWstat(const WORKER_status WST, const int CPN, simplex &blob)
{
   MASTER_comand MCM = MASwork(WST, CPN, blob);

   for (int i = 0; i < blob.vertices; i++) {
      if (blob.vertex[i].ProcID == CPN) {
         blob.vertex[i].status = WST;
         blob.vertex[i].ComEx = MCM;
         return;
      }
   }

   for (int i = 0; i < MAXAUX; i++) {
      if (blob.aux[i].ProcID == CPN) {
         blob.aux[i].status = WST;
         blob.aux[i].ComEx = MCM;
         return;
      }
   }
   GlobalPrintDebug("OH Crap. updateMWstat could not find the matching process number", 1);
}

// updateVALUEs
//
// inputs:
// CPN-Current Process Number in MASloop
// PN - total number of pararell processes
// gVALUE- well, obj fxn value
// gUNC-ehh, uncertainty for corresponding gVALUE
//
// outputs(pass by refernce):
// blob-the simplex. value and error is updated here
void
updateVALUEs(const int CPN, const double gVALUE, const double gUNC, simplex &blob)
{
   for (int i = 0; i < blob.vertices; i++) {
      if (blob.vertex[i].ProcID == CPN) {
         blob.vertex[i].value = gVALUE;
         blob.vertex[i].error = gUNC;
         return;
      }
   }

   for (int i = 0; i < MAXAUX; i++) {
      if (blob.aux[i].ProcID == CPN) {
         blob.aux[i].value = gVALUE;
         blob.aux[i].error = gUNC;
         return;
      }
   }
   GlobalPrintDebug("OH Crap. updateVALUE could not find the matching process number", 1);
}

void
SetVec(arma::vec &Tvec, const int dim, const int CPN, const simplex blob)
{
   for (int i = 0; i < blob.vertices; i++) {
      if (blob.vertex[i].ProcID == CPN) {
         for (int j = 0; j < dim; j++) {
            Tvec(static_cast<arma::uword>(j)) = blob.vertex[i].coord[j];
         }
         return;
      }
   }

   for (int i = 0; i < MAXAUX; i++ ) {
      if (blob.aux[i].ProcID == CPN) {
         for (int j = 0; j < dim; j++) {
            Tvec(static_cast<arma::uword>(j)) = blob.aux[i].coord[j];
         }
         return;
      }
   }
   GlobalPrintDebug("OH Crap. SetVec could not find the matching process number", 1);
}

// this function checks the status of each vertex and
// determines if ALL of the processors have returned error
bool
ChkErr(const simplex blob, bool III)
{
   bool FAIL = true;

   for (int i = 0; i < blob.vertices; i++) {
      if (blob.vertex[i].status != error) {
         FAIL = false;
         break;
      }
   }

   // loop over AUX
   if (!III) {
      for (int i = 0; i < MAXAUX; ++i) {
         if (blob.aux[i].status != error) {
            FAIL = false;
            break;
         }
      }
   }
   return FAIL;
}

MASTER_comand
GetMcomm(const int CPN, const simplex blob)
{
   for (int i = 0; i < blob.vertices; i++) {
      if (blob.vertex[i].ProcID == CPN) {
         return blob.vertex[i].ComEx;
      }
   }

   for (int i = 0; i < MAXAUX; i++) {
      if (blob.aux[i].ProcID == CPN) {
         return blob.aux[i].ComEx;
      }
   }
   GlobalPrintDebug("OH Crap. GetMcomm could not find the matching process number", 1);
   return question;
}

WORKER_status
getPrevWstat(const simplex blob, int CPN) {
   for (int i = 0; i < blob.vertices; i++) {
      if (blob.vertex[i].ProcID == CPN) {
         return blob.vertex[i].status;
      }
   }

   for (int i = 0; i < MAXAUX; i++) {
      if (blob.aux[i].ProcID == CPN) {
         return blob.aux[i].status;
      }
   }

   GlobalPrintDebug("OH Crap. getPrevWstat could not find the matching process number", 1);
   return undef;
}

bool
DetectChangeWstat(const simplex blob, const WORKER_status PWstat, int CPN)
{
   for (int i = 0; i < blob.vertices; i++) {
      if (blob.vertex[i].ProcID == CPN && blob.vertex[i].status != PWstat) {
         return true;
      }
   }

   for (int i = 0; i < MAXAUX; i++) {
      if (blob.aux[i].ProcID == CPN && blob.aux[i].status != PWstat) {
         return true;
      }
   }

   return false;
}

///////////////////Only other function that calls vfork() & exec
// Used for running shell script without using system(). Also
// used for generalizing the code
int
DoScript(std::string scriptStr){
   pid_t SCPTpid;
//   pid_t SCPTparent;
   int SCPTstatus;

//   SCPTparent = getpid();//parent pid
   SCPTpid = ::vfork();  // child pid

   if (SCPTpid < 0) {
      std::cout << "something is really wrong" << std::endl;
      return -1;
   }
   else if (SCPTpid == 0) {  // child process
      ::execl(scriptStr.c_str(), "sh", (char *)0);
      std::cout << "FRAME ERROR: failed at execl in DoScript" << std::endl;
      std::perror("execl");
      ::_exit(1);
   }
   else {
/*
      std::cout << "The pID of the child process is "  << pid << std::endl;
      std::cout << "The pID of the parent process is " << parent << '\n' << std::endl;
      std::cout << "waiting for process - parent" << std::endl;
*/
      ::waitpid(SCPTpid, &SCPTstatus, 0);
   }

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
