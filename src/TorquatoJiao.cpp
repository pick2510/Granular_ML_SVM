#include "TorquatoJiao.h"

//requires Boost_filesystem, glpk, and gurobi
#ifdef USE_TJ
#include <ctime>
#include <iostream>
#include <glpk.h>
#include "gurobi_c.h"
#include "gurobi_c++.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <sstream>
#include <boost/filesystem.hpp> 

using namespace std;


void TJJammer::TimerClass::Start()
{
  startTime = time(NULL);
  return;
}

void TJJammer::TimerClass::End()
{
  endTime = time(NULL);
  timePassed = endTime - startTime;
  return;
}

void TJJammer::TimerClass::Display()
{
  cout << "Elapsed time = " << timePassed << endl;
  return;
}


class LPClass 
{
public:
	LPClass(TJJammer & src)
	{
		this->pTJJammer=&src;
	}
	int  SolveLP_V00(int maxSize); //The solve function from Version0.0 (for comparison of results)
	void LPRoutine(int lpIters,
		int maxSize);    //Calls the whole process of generating, importing, solving, exporting, and cleaning up

private:
	TJJammer * pTJJammer;
	double *tempVec1;
	double *tempVec2;
	int    *ia;       //Row    index (i.e. which constraint)
	int    *ja;       //Column index (i.e. which variable)
	double *ar;
	int    *ineq;
	double *bounds;
	double coeffVal;
	double dist;
	int countCons;
	int countRows;
	int epsilonColTrack;
	int success;
	int epsilonNum;
	int ifResize;
	int j;
	gsl_matrix *updateLambdas;
	gsl_matrix *updateTemp;
	double *maxMoveVals;
	double *oldLambdasNewLocals;
	double *InternalDeltaLocalCoords;
	double maxMoveVal;
	double greatestMaxMoveVal;
	double epsilonVal;
	double updateTempVal;
	double latVol;
	int solverStatus;
	bool warmUpParam;
	int warmUpReturn;
	//GLPK variables:
	glp_prob *lp;
	glp_smcp solverParams;
	glp_bfcp otherParams;
	//Gurobi variables:
	GRBEnv     *gEnv;          //The Gurobi environment
	GRBModel   *gModel;        //The LP model
	GRBVar     *gVar;          //The variables
	GRBConstr  *gConstr;       //The constraints
	GRBLinExpr *gExpr;         //The (linear) expressions for the constraints
	double      grbExtraScale; //To push the Feasibility tolerances beyond what Gurobi usually will do

	//Functions:

	//Initialize the variables to starting values:
	void Initialize();

	//Generating constraints:
	void Generate(int maxSize); //Generates the LP constraints and imports them into the solver of choice (maxSize is the maximum number of constraints)
	void Generate_NNL(int maxSize);    //If there's a near-neighbor list to use
	void Generate_NoNNL(int maxSize);  //No near-neighbor list.
	void ApplyBulkStrainConstraint();  //Prevents volume-increasing strains (to 1st order)

	//Importing the problem into a solver:
	void Import();
	void Import_GLPK();
	//void Import_HOPDM();
	//void Import_BPMPD();
	void   Import_Gurobi();
	void   SetupEnvironment_Gurobi();
	void   AddEpsilon_Gurobi(int i,
		int l,
		int varIndex); //Functions for adding variables to gVar
	void   AddDelta_Gurobi(int InternalDeltaIndex,
		int varIndex);
	void   AddConstrs_Gurobi();
	double AddConstr_Gurobi(int constrIndex,
		int numTerms,
		int offset);
	//Solving using the LP software:
	void Solve();
	void Solve_GLPK();
	//int Solve_HOPDM();
	//int Solve_BPMPD();
	void Solve_Gurobi();

	//Exporting the solve's results:
	void Export(int lpIters);
	void ExportStrain_GLPK();
	//void ExportStrain_HOPDM();
	//void ExportStrain_BPMPD();
	void ExportStrain_Gurobi();
	void ApplyStrain();
	void ExportMovements();
	void GetMaxDisp(gsl_matrix  *globalTotDisp,
		gsl_matrix **globalDisp); //if the SLP termination criterion is maxDisp, this is called
	void ReviewResults(); //To see if any variables hit limits (especially if the non-positive trace constraint is being invoked!)
	void ExportLP_Gurobi(); //Export the Lp's solution to file
	void UpdateLPFile(int lpIters,
		double strainTrace,
		double freeSpheres);  //Write the interesting qwuantities to file

	//Cleanup: deallocation of arrays, etc.
	void Cleanup();

};
double TJJammer::getGlobalLength(double *inputVec) 
{
  // variables
  double dist = 0.0;

  for (int l=0; l<dim; l++) {
    distTempG[l] = 0.0;
    for (int d=0; d<dim; d++) {
      distTempG[l] = distTempG[l] + gsl_matrix_get(lambdas,l,d)*distTempL[d];    // calculate Lambda*localVec = global vec
    }
    dist = dist + distTempG[l]*distTempG[l];    // calculate square of global component
  }
  dist = sqrt(dist);

  // returns
  return dist;
}
gsl_matrix *getInverse(gsl_matrix *matrixToInvert) {
  // check to see if matrix is square
  if (matrixToInvert->size1 != matrixToInvert->size2) {
    printf("This matrix is not square, exiting\n");
    exit(1);
  }

  // variables needed
  int matrixSize = matrixToInvert->size1;
  int permInt = 0;                                             		       // necessary for gsl inverse calculation
  int successLU = 0;                                           			   // returns the success of the LU decomposition and inverse calc
  gsl_permutation *inverseP = gsl_permutation_calloc(matrixSize);            // necessary for gsl inverse calculation
  gsl_matrix *LUdecomp = gsl_matrix_calloc(matrixSize,matrixSize);           // necessary for the gsl inverse calculation
  gsl_matrix *inverseMatrix = gsl_matrix_calloc(matrixSize,matrixSize);      // the inverse matrix

  // calculate inverse
  gsl_matrix_memcpy(LUdecomp,matrixToInvert);

  successLU = gsl_linalg_LU_decomp(LUdecomp,inverseP,&permInt);
  if (successLU == 1) {
    cerr << "Inverse calculation of the lambdas matrix failed at LU decomposition\n";
  }
  successLU = gsl_linalg_LU_invert(LUdecomp,inverseP,inverseMatrix);
  if (successLU == 1) {
    cerr << "Inverse calculation of the lambdas matrix failed at inverse calculation\n";
  }

  // memory free
  gsl_matrix_free(LUdecomp);
  gsl_permutation_free(inverseP);

  // returns
  return inverseMatrix;
}


int TJJammer::resizeIfOverlap(bool randOver) 
{
  // variables
  int resizeSuccess = 0;           // a "success" (1) means that a resize was performed
  double dist = 0.0;
  double resize = 0.0;
  double maxResize = 0.0;
  double maxOverlap = 0.0;
  //double storeResizeDist = 0.0;
  int i,j = 0;
  NeighborNode *addressTempNNL = 0;
  double shift = 0.0;
  int quotient = 0;
  int remain = 0;
  int indexNum = int(pow(2.0*double(overBoxes) +1.0,dim));
  int selfIndexSkip = (indexNum -1)/2;              // this number is the self-image in the same unit cell
  double *randGlobalMovement = 0;
  double *globalPos = 0;

  if (useNNL) {
    // use the L/2 method if an NNL is in use
    for (i=0; i<N; i++) {
      addressTempNNL = neighborListOverlap[i];
      while (addressTempNNL) {
	j = addressTempNNL->index;
	if (i > j) {
	  addressTempNNL = addressTempNNL->next;
	  continue;
	}
	for (int d=0; d<dim; d++) {
	  distTempL[d] = localCoords[dim*i + d] - localCoords[dim*j + d];
	  if (distTempL[d] > 0.5) {
	    distTempL[d] = distTempL[d] - 1.0;
	  }
	  else if (distTempL[d] < -0.5) {
	    distTempL[d] = distTempL[d] + 1.0;
	  }
	}

	dist = getGlobalLength(distTempL);
	double overlap = radii[i] + radii[j] - (dist + resizeTol);
	// if the distance is less than the sum of the radii, then calculate and store the resize factor (resize)
	//if (dist + resizeTol < radii[i] + radii[j]) {
	if (overlap > 0) {
	  resize = ((radii[i] + radii[j])*(1.0 + resizeSpace))/dist;
	  if (resize > maxResize) {
	    maxResize = resize;
	    //storeResizeDist = dist;
	  }
	  if (overlap > maxOverlap) {maxOverlap = overlap;}
	}
	addressTempNNL = addressTempNNL->next;
      }
    }
  }//end of UseNNL if-statement
  else {//no NNL
	// calculate all (2*overBoxes +1)^dim unit cell images in the overBox method
    for (i=0; i<N; i++) {
      for (j=i; j<N; j++) {
	for (int k=0; k < indexNum; k++) {
	  if (j==i && k == selfIndexSkip) {    // skip self image in same box
	    continue;
	  }
	  quotient = k;
	  for (int d=0; d<dim; d++) {
	    remain = quotient % (2*overBoxes +1);
	    shift = double(remain - overBoxes);        // this is the image that is shift lattice cells over
	    distTempL[d] = localCoords[dim*i + d] - localCoords[dim*j + d] + shift;
	    quotient = quotient/(2*overBoxes +1);
	  }

	  dist = getGlobalLength(distTempL);
	  double overlap = radii[i] + radii[j] - (dist + resizeTol);
	  //if (dist + resizeTol < radii[i] + radii[j]) {// note that this is for additive diameters!
	  if (overlap > 0) {
	    resize = ((radii[i] + radii[j])*(1.0 + resizeSpace))/dist;
	    if (resize > maxResize) {
	      maxResize = resize;
	      //storeResizeDist = dist;
	    }
	    if (overlap > maxOverlap) {maxOverlap = overlap;}
	  }
	}
      }
    }
  }

  // rescale the lambdas matrix lattice vectors by maxResize so that no detected overlap is present
  if (maxResize > 0.0) {
	  if(Verbosity>7)
		  cout << "  resizeTol requests pack resizing (maxOverlap = " << maxOverlap << ")" << endl;
	  resizeSuccess = 1;

	  // if no random overlap is desired, then simply rescale so that farthest overlap is now contacting
	  if (!randOver) {
		  maxResize /= (1.0 + resizeSpace);
    }

    // rescale the lambdas matrix
    //gsl_matrix_scale(lambdas,maxResize + termTol); //Original--CAN'T USE termTol!!!
    //gsl_matrix_scale(lambdas,maxResize + resizeTol);
    gsl_matrix_scale(lambdas,maxResize);

    // if stochastic moves are used, then move spheres around a little after resize
    if (randOverlapMove && randOver) {
      // variables
      randGlobalMovement = new double[dim];
      globalPos = new double[dim];
      double randLength = 0.0;
      //gsl_matrix *inverseLambdas = gsl_matrix_calloc(dim,dim);

      // get the inverse of lambdas
      inverseLambdas = getInverse(lambdas);

      // get the random movement
      for (i=0; i<N; i++) {
	randLength = 0.0;
	for (int d=0; d<dim; d++) {
	  randGlobalMovement[d] = 2.0*(gsl_rng_uniform_pos(RNG) - 0.5);
	  randLength = randLength + randGlobalMovement[d]*randGlobalMovement[d];
	}
	randLength = sqrt(randLength);

	// get the new global position
	for (int d=0; d<dim; d++) {
	  globalPos[d] = 0.0;
	  for (j=0; j<dim; j++) {
	    globalPos[d] = globalPos[d] + gsl_matrix_get(lambdas,d,j)*localCoords[dim*i +j];
	  }
	  globalPos[d] = globalPos[d] + ((radii[i]*resizeSpace)/randLength)*randGlobalMovement[d];
	}

	// replace local coordinates with converted global coordinates
	for (int d=0; d<dim; d++) {
	  localCoords[dim*i +d] = 0.0;
	  for (j=0; j<dim; j++) {
	    localCoords[dim*i +d] = localCoords[dim*i +d] + gsl_matrix_get(inverseLambdas,d,j)*globalPos[j];
	  }
	  if (localCoords[dim*i +d] >= 1.0) {
	    localCoords[dim*i +d] = localCoords[dim*i +d] -1.0;
	  }
	  else if (localCoords[dim*i +d] < 0.0) {
	    localCoords[dim*i +d] = localCoords[dim*i +d] +1.0;
	  }
	}
      }

      // deletions
      delete [] randGlobalMovement;
      delete [] globalPos;
    }
  }

  // returns
  return resizeSuccess;
}

double TJJammer::getVol() {//Gets the volume of the fundamental cell
  // variables
  double volume = 0.0;
  int permInt = 0;                                             		       // necessary for gsl LU decomp and gsl determinant
  int successLU = 0;                                           			   // returns the success of the LU decomposition
  gsl_permutation *inverseP = gsl_permutation_calloc(dim);            	   // necessary for gsl LU decomp
  gsl_matrix *LUdecomp = gsl_matrix_calloc(dim,dim);          			   // necessary for gsl LU decomp

  gsl_matrix_memcpy(LUdecomp,lambdas);

  successLU = gsl_linalg_LU_decomp(LUdecomp,inverseP,&permInt);
  if (successLU == 1) {
    cerr << "Inverse calculation of the lambdas matrix failed at LU decomposition\n";
  }

  volume = gsl_linalg_LU_det(LUdecomp, permInt);

  // memory free
  gsl_matrix_free(LUdecomp);
  gsl_permutation_free(inverseP);

  // returns
  return fabs(volume);
}
double TJJammer::getSphereVol()
//Gets the total volume of the spheres, sphereVol
//phi = sphereVol / latVol
{
  double sphereVol = 0.0;
  double sphereVolConst = pow(pi,dim/2.0)/gsl_sf_gamma(dim/2.0 +1.0);

  for (int i = 0 ; i < N ; i++) {
    sphereVol += sphereVolConst * pow(radii[i],dim);
  }

  return(sphereVol);
}


//////////////////////////////////////////////////////
//THE OLD FUNCTION:

//int LPClass::SolveLP_V00(int maxSize)
////The LP-solving routine from Version 0.0
//{
//  // organization of variables is as follows: strain matrix (11, 12, 13, 14, 22, 23, 24, 33, 34, 44) then particle nums. This is an N*pTJJammer->dim + (pTJJammer->dim*(pTJJammer->dim+1)/2) column constraint matrix
//  // variables
//  double *tempVec1 = new double[pTJJammer->dim];
//  double *tempVec2 = new double[pTJJammer->dim];
//  int *ia;
//  int *ja;
//  double *ar;
//  double *bounds;
//  double coeffVal = 0.0;
//  double dist = 0.0;
//  int countCons = 1;
//  int countRows = 0;
//  int epsilonColTrack = 0;
//  int success = 0;
//  int epsilonNum = (pTJJammer->dim*(pTJJammer->dim +1))/2;
//  int ifResize = 0;
//  int j = 0;
//  gsl_matrix *updateLambdas = gsl_matrix_calloc(pTJJammer->dim,pTJJammer->dim);
//  gsl_matrix *updateTemp = gsl_matrix_calloc(pTJJammer->dim,pTJJammer->dim);
//  double *maxMoveVals = new double[pTJJammer->dim];
//  double *oldLambdasNewLocals = new double[pTJJammer->dim];
//  double *InternalDeltaLocalCoords = new double[pTJJammer->dim];
//  double maxMoveVal = 0.0;
//  double greatestMaxMoveVal = 0.0;
//  double epsilonVal = 0.0;
//  double updateTempVal = 0.0;
//  double latVol = 0.0;
//  int solverStatus = 0;
//  bool warmUpParam = false;
//  int warmUpReturn = 0;
//  glp_prob *lp;
//  glp_smcp solverParams;
//  glp_bfcp otherParams;
//
//  for (int i=0; i<pTJJammer->dim; i++) {
//    maxMoveVals[i] = 0.0;
//    oldLambdasNewLocals[i] = 0.0;
//    InternalDeltaLocalCoords[i] = 0.0;
//  }
//
//  // this segment calculates the constraints (coefficients of the "row" variables in the LP solver)
//  if (pTJJammer->useNNL) {
//    // NNL variables
//    TJJammer::NeighborNode *addressTempNNL = 0;
//
//    // LP solver variables
//    ia = new int[(2*pTJJammer->dim + epsilonNum)*maxSize +1];
//    ja = new int[(2*pTJJammer->dim + epsilonNum)*maxSize +1];
//    ar = new double[(2*pTJJammer->dim + epsilonNum)*maxSize +1];
//    bounds = new double[(2*pTJJammer->dim + epsilonNum)*maxSize +1];
//    for (int i=0; i<(2*pTJJammer->dim + epsilonNum)*maxSize +1; i++) {
//      ia[i] = 0;
//      ja[i] = 0;
//      ar[i] = 0.0;
//      bounds[i] = 0.0;
//    }
//
//    // solution; use L/2 method to calc distance
//    for (int i=0; i<pTJJammer->N; i++) {
//      addressTempNNL = pTJJammer->neighborListDelta[i];
//      while (addressTempNNL) {
//	j = addressTempNNL->index;
//
//	// continue if j < i, i.e., only perform NNL calcs once for each contact
//	if (j < i) {
//	  addressTempNNL = addressTempNNL->next;
//	  continue;
//	}
//
//	for (int d=0; d<pTJJammer->dim; d++) {
//	  pTJJammer->distTempL[d] = pTJJammer->localCoords[pTJJammer->dim*i + d] - pTJJammer->localCoords[pTJJammer->dim*j + d];
//	  if (pTJJammer->distTempL[d] > 0.5) {
//	    pTJJammer->distTempL[d] = pTJJammer->distTempL[d] - 1.0;
//	  }
//	  else if (pTJJammer->distTempL[d] < -0.5) {
//	    pTJJammer->distTempL[d] = pTJJammer->distTempL[d] + 1.0;
//	  }
//	}
//
//	dist = pTJJammer->getGlobalLength(pTJJammer->distTempL);
//
//	if (dist <= pTJJammer->radii[i] + pTJJammer->radii[j] + pTJJammer->InternalDelta) {
//	  countRows = countRows +1;
//
//	  // calc tempVec1 (Lambdas*pTJJammer->distTempL)
//	  for (int l=0; l<pTJJammer->dim; l++) {
//	    tempVec1[l] = 0.0;
//	    for (int d=0; d<pTJJammer->dim; d++) {
//	      tempVec1[l] = tempVec1[l] + gsl_matrix_get(pTJJammer->lambdas,l,d)*pTJJammer->distTempL[d];
//	    }
//	  }
//
//	  // input the epsilon coeffs
//	  epsilonColTrack = 1;
//	  for (int d=0; d<pTJJammer->dim; d++) {
//	    ar[countCons] = tempVec1[d]*tempVec1[d];
//	    ia[countCons] = countRows;
//	    ja[countCons] = epsilonColTrack;
//	    countCons = countCons +1;
//	    epsilonColTrack = epsilonColTrack +1;
//	    for (int l=d +1; l<pTJJammer->dim; l++) {
//	      ar[countCons] = 2*tempVec1[d]*tempVec1[l];
//	      ia[countCons] = countRows;
//	      ja[countCons] = epsilonColTrack;
//	      countCons = countCons +1;
//	      epsilonColTrack = epsilonColTrack +1;
//	    }
//	  }
//
//	  // calc tempVec2 (Lambdas_T*tempVec1)
//	  for (int l=0; l<pTJJammer->dim; l++) {
//	    tempVec2[l] = 0.0;
//	    for (int d=0; d<pTJJammer->dim; d++) {
//	      tempVec2[l] =  tempVec2[l] + gsl_matrix_get(pTJJammer->lambdas,d,l)*tempVec1[d];
//	    }
//	  }
//
//	  // input the pTJJammer->InternalDeltaR coeffs (order is: x1,x2,y1,y2,z1,z2...)
//	  for (int d=0; d<pTJJammer->dim; d++) {
//	    ar[countCons] = tempVec2[d];
//	    ia[countCons] = countRows;
//	    ja[countCons] = epsilonNum + i*pTJJammer->dim +d +1;
//	    countCons = countCons +1;
//	    ar[countCons] = -1.0*tempVec2[d];
//	    ia[countCons] = countRows;
//	    ja[countCons] = epsilonNum + j*pTJJammer->dim +d +1;
//	    countCons = countCons +1;
//	  }
//
//	  // calc bounds (i.e., RHS of constraint)
//	  coeffVal = ((pTJJammer->radii[i] + pTJJammer->radii[j])*(pTJJammer->radii[i] + pTJJammer->radii[j]));    // note this is for additive diameters!
//	  for (int d=0; d<pTJJammer->dim; d++) {
//	    coeffVal = coeffVal - pTJJammer->distTempL[d]*tempVec2[d];
//	  }
//	  bounds[countRows] = coeffVal/2.0;
//	}
//
//	addressTempNNL = addressTempNNL->next;
//      }
//    }
//  }
//  else {
//    // box variables
//    int indexNum = int(pow(2.0*double(pTJJammer->overBoxes) +1.0,pTJJammer->dim));
//    int selfIndexSkip = (indexNum -1)/2;             // this number is the self-image in the same unit cell
//    double shift = 0.0;
//    int quotient = 0;
//    int remain = 0;
//
//    // LP solver variables
//    maxSize = ((pTJJammer->N*(pTJJammer->N+1))/2)*indexNum -1;             // overkill if pTJJammer->InternalDelta is not large. But if it is large, this is barely sufficient!
//    ia = new int[(2*pTJJammer->dim + epsilonNum)*maxSize +1];
//    ja = new int[(2*pTJJammer->dim + epsilonNum)*maxSize +1];
//    ar = new double[(2*pTJJammer->dim + epsilonNum)*maxSize +1];
//    bounds = new double[(2*pTJJammer->dim + epsilonNum)*maxSize +1];
//    for (int i=0; i<(2*pTJJammer->dim + epsilonNum)*maxSize +1; i++) {
//      ia[i] = 0;
//      ja[i] = 0;
//      ar[i] = 0.0;
//      bounds[i] = 0.0;
//    }
//
//    for (int i=0; i<pTJJammer->N; i++) {
//      for (j=i; j<pTJJammer->N; j++) {
//	for (int k=0; k < indexNum; k++) {
//	  if (j==i && k == selfIndexSkip) {    // skip self image in same box
//	    continue;
//	  }
//	  quotient = k;
//	  for (int d=0; d<pTJJammer->dim; d++) {
//	    remain = quotient % (2*pTJJammer->overBoxes +1);
//	    shift = double(remain - pTJJammer->overBoxes);        // this is the image that is shift lattice cells over
//	    pTJJammer->distTempL[d] = pTJJammer->localCoords[pTJJammer->dim*i + d] - pTJJammer->localCoords[pTJJammer->dim*j + d] + shift;
//	    quotient = quotient/(2*pTJJammer->overBoxes +1);
//	  }
//
//	  dist = pTJJammer->getGlobalLength(pTJJammer->distTempL);
//
//	  if (dist <= pTJJammer->radii[i] + pTJJammer->radii[j] + pTJJammer->InternalDelta) {
//	    countRows = countRows +1;
//
//	    // calc tempVec1 (Lambdas*pTJJammer->distTempL)
//	    for (int l=0; l<pTJJammer->dim; l++) {
//	      tempVec1[l] = 0.0;
//	      for (int d=0; d<pTJJammer->dim; d++) {
//		tempVec1[l] = tempVec1[l] + gsl_matrix_get(pTJJammer->lambdas,l,d)*pTJJammer->distTempL[d];
//	      }
//	    }
//
//	    // input the epsilon coeffs
//	    epsilonColTrack = 1;
//	    for (int d=0; d<pTJJammer->dim; d++) {
//	      ar[countCons] = tempVec1[d]*tempVec1[d];
//	      ia[countCons] = countRows;
//	      ja[countCons] = epsilonColTrack;
//	      countCons = countCons +1;
//	      epsilonColTrack = epsilonColTrack +1;
//	      for (int l=d +1; l<pTJJammer->dim; l++) {
//		ar[countCons] = 2*tempVec1[d]*tempVec1[l];
//		ia[countCons] = countRows;
//		ja[countCons] = epsilonColTrack;
//		countCons = countCons +1;
//		epsilonColTrack = epsilonColTrack +1;
//	      }
//	    }
//
//	    // calc tempVec2 (Lambdas_T*tempVec1)
//	    for (int l=0; l<pTJJammer->dim; l++) {
//	      tempVec2[l] = 0.0;
//	      for (int d=0; d<pTJJammer->dim; d++) {
//		tempVec2[l] =  tempVec2[l] + gsl_matrix_get(pTJJammer->lambdas,d,l)*tempVec1[d];
//	      }
//	    }
//
//	    // input the pTJJammer->InternalDeltaR coeffs (order is: x1,x2,y1,y2,z1,z2...)
//	    if (i!=j) {              // if i==j, then pTJJammer->InternalDeltaR coeffs are immaterial b/c coeff*(pTJJammer->InternalDeltaR1 - pTJJammer->InternalDeltaR1) = 0
//	      for (int d=0; d<pTJJammer->dim; d++) {
//		ar[countCons] = tempVec2[d];
//		ia[countCons] = countRows;
//		ja[countCons] = epsilonNum + i*pTJJammer->dim +d +1;
//		countCons = countCons +1;
//		ar[countCons] = -1.0*tempVec2[d];
//		ia[countCons] = countRows;
//		ja[countCons] = epsilonNum + j*pTJJammer->dim +d +1;
//		countCons = countCons +1;
//	      }
//	    }
//
//	    // calc bounds (i.e., RHS of constraint)
//	    coeffVal = ((pTJJammer->radii[i] + pTJJammer->radii[j])*(pTJJammer->radii[i] + pTJJammer->radii[j]));    // note this is for additive diameters!
//	    for (int d=0; d<pTJJammer->dim; d++) {
//	      coeffVal = coeffVal - pTJJammer->distTempL[d]*tempVec2[d];
//	    }
//	    bounds[countRows] = coeffVal/2.0;
//	  }
//	}
//      }
//    }
//  }
//
//  // create problem and read values in
//  lp = glp_create_prob();
//  glp_set_prob_name(lp,"SLP");
//  glp_set_obj_dir(lp,GLP_MIN);
//
//  // specify control parameters: IMPORTANT!!!
//  glp_init_smcp(&solverParams);
//  solverParams.presolve = pTJJammer->GLPK_PresolveVar;
//  solverParams.tol_bnd = pTJJammer->feasibleTol;
//  solverParams.tol_dj = pTJJammer->feasibleTol;
//  solverParams.msg_lev = GLP_MSG_ERR; //Only print error and warning messages to the terminal
//
//  glp_get_bfcp(lp,&otherParams);
//  otherParams.type = pTJJammer->GLPK_BasisFactType;
//  glp_set_bfcp(lp, &otherParams);
//
//
//  // first do the objective function and structural variables (i.e., position displacement and epsilons)
//  glp_add_cols(lp,pTJJammer->dim*pTJJammer->N + epsilonNum);
//  epsilonColTrack = 1;
//  for (int i=0; i<pTJJammer->dim; i++) {
//    for (int l=i; l<pTJJammer->dim; l++) {
//      if (l==i) {
//	glp_set_obj_coef(lp,epsilonColTrack,1.0);
//	glp_set_col_bnds(lp,epsilonColTrack,GLP_DB,-1.0*pTJJammer->InternalCompMax,0.0);
//      }
//      else {
//	if (pTJJammer->strictJam == 1) {
//	  glp_set_col_bnds(lp,epsilonColTrack,GLP_DB,-1.0*pTJJammer->InternalShearMax,pTJJammer->InternalShearMax);
//	}
//	else {
//	  glp_set_col_bnds(lp,epsilonColTrack,GLP_DB,0.0,0.0);
//	}
//      }
//      epsilonColTrack = epsilonColTrack +1;
//    }
//  }
//
//  for (int i=epsilonNum +1; i<=pTJJammer->dim*pTJJammer->N + epsilonNum; i++) {
//    glp_set_col_bnds(lp,i,GLP_DB,-1.0*pTJJammer->InternalTransMax,1.0*pTJJammer->InternalTransMax);
//  }
//
//  // now do the constraints lower bounds
//  glp_add_rows(lp,countRows);
//  for (int i=1; i<=countRows; i++) {
//    glp_set_row_bnds(lp,i,GLP_LO,bounds[i],0.0);
//  }
//
//  // now load the problem and solve it
//  //#load
//  cout << "  Loading and solving the LP..." << endl;
//  glp_load_matrix(lp,countCons-1,ia,ja,ar);
//
//  // in case the basis matrix was singular on the last run
//  if (warmUpParam) {
//    warmUpReturn = glp_warm_up(lp);
//    if (warmUpReturn == 0) {
//      cout << "warm up successfully fixed the basis matrix\n";
//    }
//    else {
//      cout << "warm up indicated that the basis matrix is numerically singular within the precision specified\n";
//    }
//  }
//
//  solverStatus = glp_simplex(lp,&solverParams);
//  cout << "  Complete!" << endl;
//
//  if (solverStatus > 0) {
//    cout << "solver encountered problems, running warm up to attempt to resolve\n";
//    warmUpParam = true;
//  }
//  else {
//    warmUpParam = false;
//  }
//
//  // now update the coordinates and lattice
//  // first the lattice vecs; start by creating the matrix (I + epsilon), which is updateLambdas
//  epsilonColTrack = 1;
//  for (int i=0; i<pTJJammer->dim; i++) {
//    for (int d=i; d<pTJJammer->dim; d++) {
//      if (i==d) {
//	gsl_matrix_set(updateLambdas,i,d,1.0 + glp_get_col_prim(lp,epsilonColTrack));
//      }
//      else {
//	gsl_matrix_set(updateLambdas,i,d,glp_get_col_prim(lp,epsilonColTrack));
//	gsl_matrix_set(updateLambdas,d,i,glp_get_col_prim(lp,epsilonColTrack));
//      }
//      epsilonColTrack = epsilonColTrack +1;
//    }
//  }
//
//  // now multiply: pTJJammer->lambdas = (I + epsilon)*pTJJammer->lambdas
//  for (int l=0; l<pTJJammer->dim; l++) {
//    for (int i=0; i<pTJJammer->dim; i++) {
//      updateTempVal = 0.0;
//      for (int d=0; d<pTJJammer->dim; d++) {
//	updateTempVal = updateTempVal + gsl_matrix_get(updateLambdas,l,d)*gsl_matrix_get(pTJJammer->lambdas,d,i);
//      }
//      gsl_matrix_set(updateTemp,l,i,updateTempVal);   // now updateTemp is the new pTJJammer->lambdas matrix, (I+epsilon)*pTJJammer->lambdas
//    }
//  }
//
//  // now update the coordinates
//  greatestMaxMoveVal = 0.0;
//  for (int i=0; i<pTJJammer->N; i++) {
//    for (int d=0; d<pTJJammer->dim; d++) {
//      InternalDeltaLocalCoords[d] = glp_get_col_prim(lp,epsilonNum +1 +pTJJammer->dim*i +d);  //Translation in terms of the (new) lattice
//      pTJJammer->localCoords[pTJJammer->dim*i +d] = pTJJammer->localCoords[pTJJammer->dim*i +d] + InternalDeltaLocalCoords[d]; //Increment
//      //Apply PBCs
//      if (pTJJammer->localCoords[pTJJammer->dim*i +d] > 1.0) {
//	pTJJammer->localCoords[pTJJammer->dim*i +d] = pTJJammer->localCoords[pTJJammer->dim*i +d] - 1.0;
//      }
//      else if (pTJJammer->localCoords[pTJJammer->dim*i +d] < 0.0) {
//	pTJJammer->localCoords[pTJJammer->dim*i +d] = pTJJammer->localCoords[pTJJammer->dim*i +d] + 1.0;
//      }
//    }
//
//    // now calc the maximum movement - only if NNLs are used
//    if (pTJJammer->useNNL) {
//      maxMoveVal = 0.0;
//      for (int l=0; l<pTJJammer->dim; l++) {
//	oldLambdasNewLocals[l] = 0.0;
//	for (int d=0; d<pTJJammer->dim; d++) {
//	  oldLambdasNewLocals[l] = oldLambdasNewLocals[l] + gsl_matrix_get(pTJJammer->lambdas,l,d)*pTJJammer->localCoords[pTJJammer->dim*i +d];
//	}
//      }
//      for (int l=0; l<pTJJammer->dim; l++) {
//	maxMoveVals[l] = 0.0;
//	for (int d=0; d<pTJJammer->dim; d++) {
//	  if (l==d) {
//	    epsilonVal = gsl_matrix_get(updateLambdas,l,d) - 1.0;
//	  }
//	  else {
//	    epsilonVal = gsl_matrix_get(updateLambdas,l,d);
//	  }
//	  maxMoveVals[l] = maxMoveVals[l] + epsilonVal*oldLambdasNewLocals[d];
//	}
//	maxMoveVals[l] = maxMoveVals[l] + InternalDeltaLocalCoords[l];
//	maxMoveVal = maxMoveVal + maxMoveVals[l]*maxMoveVals[l];
//      }
//      if (maxMoveVal > greatestMaxMoveVal) {
//	greatestMaxMoveVal = maxMoveVal;
//      }
//    }
//  }
//
//  if (pTJJammer->useNNL) {
//    greatestMaxMoveVal = sqrt(greatestMaxMoveVal);       // the scaling is necessary to match InternalNNLextraDist
//    pTJJammer->maxSingleMove = pTJJammer->maxSingleMove + greatestMaxMoveVal;
//    cout << "  greatestMaxMoveVal = " << greatestMaxMoveVal << endl;
//    cout << "  maxSingleMove      = " << pTJJammer->maxSingleMove << "\n";
//  }
//
//  // now update the pTJJammer->lambdas matrix officially
//  gsl_matrix_memcpy(pTJJammer->lambdas,updateTemp);
//
//  // delete the LP problem
//  glp_delete_prob(lp);
//
//  /*cout << "New coordinates and pTJJammer->lambdas, pre-resize check:\n";
//    for (int i=0; i<pTJJammer->dim; i++) {
//    for (j=i; j<pTJJammer->dim; j++) {
//    cout << "e" << i << j << " = " << gsl_matrix_get(pTJJammer->lambdas,i,j) << "\n";
//    }
//    }
//
//    for (int i=0; i<N; i++) {
//    cout << "x = " << pTJJammer->localCoords[pTJJammer->dim*i] << ";  y = " << pTJJammer->localCoords[pTJJammer->dim*i +1] << ";  z = " << pTJJammer->localCoords[pTJJammer->dim*i +2] << "\n";
//    }*/
//
//  latVol = pTJJammer->getVol();
//  cout << "  lattice volume     = " << latVol << "\n";
//
//  // now check for overlap
//  if (pTJJammer->lpIters == pTJJammer->maxIters -1) {                 // if it's the last iteration, resize only so greatest overlap is at contact!
//    ifResize = pTJJammer->resizeIfOverlap(false);
//  }
//  else {
//    ifResize = pTJJammer->resizeIfOverlap(true);
//  }
//
//  if (ifResize) {
//    //cout << "overlap detected, resize necessary\n";
//    latVol = pTJJammer->getVol();
//    cout << "  new lattice volume is " << latVol << "\n";
//  }
//
//
//  // free memory
//  delete [] ia;
//  delete [] ja;
//  delete [] ar;
//  delete [] bounds;
//  delete [] tempVec1;
//  delete [] tempVec2;
//  delete [] maxMoveVals;
//  delete [] oldLambdasNewLocals;
//  delete [] InternalDeltaLocalCoords;
//  gsl_matrix_free(updateLambdas);
//  gsl_matrix_free(updateTemp);
//
//  // returns
//  success = 1;
//  return success;
//
//}
//

void LPClass::LPRoutine(int lpIters,
		       int maxSize)
{
  Initialize();
  Generate(maxSize);
  Import();
  Solve(); //Sets solverStatus
  Export(lpIters);
  Cleanup();

  return;
}


//////////////////////////////////////////////////////
//INITIALIZE OBJECT VARIABLES


void LPClass::Initialize()
{
  tempVec1 = new double[pTJJammer->dim];
  tempVec2 = new double[pTJJammer->dim];
  //int *ia;
  //int *ja;
  //double *ar;
  //double *bounds;
  coeffVal = 0.0;
  dist = 0.0;
  countCons = 1;
  countRows = 0;
  epsilonColTrack = 0;
  success = 0;
  epsilonNum = pTJJammer->strictJam == 1 ? (pTJJammer->dim*(pTJJammer->dim +1))/2 : 1; //epsilon is a scalar for collective jamming!
  ifResize = 0;
  j = 0;
  updateLambdas = gsl_matrix_calloc(pTJJammer->dim,pTJJammer->dim);
  updateTemp    = gsl_matrix_calloc(pTJJammer->dim,pTJJammer->dim);
  maxMoveVals = new double[pTJJammer->dim];
  oldLambdasNewLocals = new double[pTJJammer->dim];
  InternalDeltaLocalCoords = new double[pTJJammer->dim];
  maxMoveVal = 0.0;
  greatestMaxMoveVal = 0.0;
  epsilonVal = 0.0;
  updateTempVal = 0.0;
  latVol = 0.0;
  solverStatus = 0;
  warmUpParam = false;
  warmUpReturn = 0;
  //glp_prob *lp;
  //glp_smcp solverParams;
  //glp_bfcp otherParams;


  for (int i=0; i<pTJJammer->dim; i++) {
    maxMoveVals[i] = 0.0;
    oldLambdasNewLocals[i] = 0.0;
    InternalDeltaLocalCoords[i] = 0.0;
  }

  return;
}


//////////////////////////////////////////////////////
//GENERATE THE LP


void LPClass::Generate(int maxSize)
//Generates the LP constraints and imports them into the solver of choice (maxSize is the maximum number of constraints)
{
  //cout << "  Generate Constraints" << endl;
  if (pTJJammer->useNNL) {Generate_NNL  (maxSize);}
  else        {Generate_NoNNL(maxSize);}

  return;
}


void LPClass::Generate_NNL(int maxSize)
//NOTE: ia,ja,ar START AT INDEX 1 (NOT ZERO!)
{
  // NNL variables
  TJJammer::NeighborNode *addressTempNNL = 0;

  // LP solver variables
  //Each sphere pair constraint equation has the following entries:
  //  pTJJammer->dim entries for each sphere's displacement (total = 2*pTJJammer->dim)
  //  epsilonNum entries for the deforming box
  //    (If collective jamming, then this is a scalar)
  const int arraySize = pTJJammer->strictJam ? 
    (2*pTJJammer->dim + epsilonNum)*maxSize+pTJJammer->dim+1 : //+pTJJammer->dim for the bulk strain constraint (Trace <= 0), and +1 because the arrays start at index 1
    (2*pTJJammer->dim + epsilonNum)*maxSize    +1;  //Don't worry--epsilonNum=1 if strictJam==0 

  //cout << "Allocate for " << arraySize << " entries" << endl;
  ia     = new int   [arraySize];
  ja     = new int   [arraySize];
  ar     = new double[arraySize];
  ineq   = new int   [arraySize];
  bounds = new double[arraySize];
  for (int i=0 ; i<arraySize ; i++) {
    ia[i] = 0;
    ja[i] = 0;
    ar[i] = 0.0;
    ineq[i] = 0;
    bounds[i] = 0.0;
  }

  //Add the bulk strain constraint:
  if (pTJJammer->strictJam == 1)
    ApplyBulkStrainConstraint();

  //use L/2 method to calc distance
  //int startIndex = 1;
  for (int i=0; i<pTJJammer->N; i++) {
    addressTempNNL = pTJJammer->neighborListDelta[i];
    while (addressTempNNL) {
      j = addressTempNNL->index;

      // continue if j < i, i.e., only perform NNL calcs once for each contact
      if (j < i) {
	addressTempNNL = addressTempNNL->next;
	continue;
      }

      for (int d=0; d<pTJJammer->dim; d++) {
	pTJJammer->distTempL[d] = pTJJammer->localCoords[pTJJammer->dim*i + d] - pTJJammer->localCoords[pTJJammer->dim*j + d];
	if (pTJJammer->distTempL[d] > 0.5) {
	  pTJJammer->distTempL[d] = pTJJammer->distTempL[d] - 1.0;
	}
	else if (pTJJammer->distTempL[d] < -0.5) {
	  pTJJammer->distTempL[d] = pTJJammer->distTempL[d] + 1.0;
	}
      }

      dist = pTJJammer->getGlobalLength(pTJJammer->distTempL);
      double distSq = dist*dist;

      if (dist <= pTJJammer->radii[i] + pTJJammer->radii[j] + pTJJammer->InternalDelta) {//If the two spheres are within influence sphere-range
	//cout << "   Interaction: " << i << " & " << j << endl;
	countRows++;

	//
	//
	//  ---Here's how all the "temps" stack up:---
	//
	//  distTempL * Lambdas_T * Lambdas * distTempL 
	//                          \_________________/
	//                               tempVec1
	//              \_____________________________/
	//                         tempVec2
	//  \_________________________________________/
	//        Global distance squared btwn 
	//           sphere centers i and j
	//
	//

	// calc tempVec1 (Lambdas*pTJJammer->distTempL)
	// and 
	for (int l=0; l<pTJJammer->dim; l++) {
	  tempVec1[l] = 0.0;
	  for (int d=0; d<pTJJammer->dim; d++) {
	    tempVec1[l] += gsl_matrix_get(pTJJammer->lambdas,l,d)*pTJJammer->distTempL[d];
	  }
	}

	// calc tempVec2 (Lambdas_T*tempVec1)
	for (int l=0; l<pTJJammer->dim; l++) {
	  tempVec2[l] = 0.0;
	  for (int d=0; d<pTJJammer->dim; d++) {
	    tempVec2[l] += gsl_matrix_get(pTJJammer->lambdas,d,l)*tempVec1[d];
	  }
	}

	// input the epsilon coeffs
	if (pTJJammer->strictJam == 0) {//Simple for collective jamming!
	  //#check
	  ar[countCons] = distSq; //Global distance squared between the spheres
	  ia[countCons] = countRows;
	  ja[countCons] = 1;
	  countCons++;
	}
	else if (pTJJammer->strictJam == 1) {
	  epsilonColTrack = 1;
	  for (int d=0; d<pTJJammer->dim; d++) {//First index of the epsilon
	    ar[countCons] = tempVec1[d]*tempVec1[d];
	    ia[countCons] = countRows;
	    ja[countCons] = epsilonColTrack;
	    countCons++;
	    epsilonColTrack = epsilonColTrack +1;
	    for (int l=d +1; l<pTJJammer->dim; l++) {//Second index of the epsilon
	      ar[countCons] = 2*tempVec1[d]*tempVec1[l];
	      ia[countCons] = countRows;
	      ja[countCons] = epsilonColTrack;
	      countCons++;
	      epsilonColTrack = epsilonColTrack +1;
	    }
	  }
	}



	// input the pTJJammer->InternalDeltaR coeffs (order is: x1,x2,y1,y2,z1,z2...)
	for (int d=0; d<pTJJammer->dim; d++) {
	  ar[countCons] = tempVec2[d];
	  ia[countCons] = countRows;
	  ja[countCons] = epsilonNum + i*pTJJammer->dim +d +1;
	  countCons++;
	  ar[countCons] = -1.0*tempVec2[d];
	  ia[countCons] = countRows;
	  ja[countCons] = epsilonNum + j*pTJJammer->dim +d +1;
	  countCons++;
	}

	// calc bounds (i.e., RHS of constraint)
	coeffVal = ((pTJJammer->radii[i] + pTJJammer->radii[j])*(pTJJammer->radii[i] + pTJJammer->radii[j])) - distSq;    // note this is for additive diameters!
	//for (int d=0; d<pTJJammer->dim; d++) {
	//coeffVal -= pTJJammer->distTempL[d]*tempVec2[d];
	//}
	ineq  [countRows] = GLP_LO;
	bounds[countRows] = coeffVal/2.0;

	//Check:
	//int lastIndex = countCons;
	//cout << "  Constraint " << countRows+1 << ":" << endl;
	//for (int checkIndex = startIndex ; checkIndex < lastIndex ; checkIndex++) {
	//cout << "A(" << ia[checkIndex] << "," << ja[checkIndex] << ") = " << ar[checkIndex] << endl;
	//}
	//cout << "b(" << countRows << ") = " << bounds[countRows] << endl;
	//startIndex = lastIndex;
      }

      addressTempNNL = addressTempNNL->next;
    }
  }
  if(Verbosity>7)
	  cout << "  Generated " << countRows << " constraints, " << countCons << " entries" << endl;

  return;
}


void LPClass::Generate_NoNNL(int maxSize)
{
  // box variables
  int indexNum = int(pow(2.0*double(pTJJammer->overBoxes) +1.0,pTJJammer->dim));
  int selfIndexSkip = (indexNum -1)/2;             // this number is the self-image in the same unit cell
  double shift = 0.0;
  int quotient = 0;
  int remain = 0;

  // LP solver variables
  maxSize = ((pTJJammer->N*(pTJJammer->N+1))/2)*indexNum -1;             // overkill if pTJJammer->InternalDelta is not large. But if it is large, this is barely sufficient!
  const int arraySize = pTJJammer->strictJam ? 
    (2*pTJJammer->dim + epsilonNum)*maxSize+pTJJammer->dim+1 : //+pTJJammer->dim for the bulk strain constraint (Trace <= 0), and +1 because the arrays start at index 1
    (2*pTJJammer->dim + epsilonNum)*maxSize    +1;  //Don't worry--epsilonNum=1 if strictJam==0


  if(Verbosity>7)
  cout << "  ia,ja,ar,... arraySize = " << arraySize << endl;
  ia     = new int   [arraySize];
  ja     = new int   [arraySize];
  ar     = new double[arraySize];
  ineq   = new int   [arraySize];  
  bounds = new double[arraySize];
  for (int i=0 ; i<arraySize ; i++) {
    ia[i] = 0;
    ja[i] = 0;
    ar[i] = 0.0;
    ineq[i] = 0;
    bounds[i] = 0.0;
  }

  //Add the bulk strain constraint:
  //  (Unneeded in collective jamming because the strain variable is a scalar)
  if (pTJJammer->strictJam==1)
    ApplyBulkStrainConstraint();

  for (int i=0; i<pTJJammer->N; i++) {
    for (j=i; j<pTJJammer->N; j++) {
      for (int k=0; k < indexNum; k++) {
	if (j==i && k == selfIndexSkip) {    // skip self image in same box
	  continue;
	}
	quotient = k;
	for (int d=0; d<pTJJammer->dim; d++) {
	  remain = quotient % (2*pTJJammer->overBoxes +1);
	  shift = double(remain - pTJJammer->overBoxes);        // this is the image that is shift lattice cells over
	  pTJJammer->distTempL[d] = pTJJammer->localCoords[pTJJammer->dim*i + d] - pTJJammer->localCoords[pTJJammer->dim*j + d] + shift;
	  quotient = quotient/(2*pTJJammer->overBoxes +1);
	}

	dist = pTJJammer->getGlobalLength(pTJJammer->distTempL);
	double distSq = dist*dist;

	if (dist <= pTJJammer->radii[i] + pTJJammer->radii[j] + pTJJammer->InternalDelta) {
	  //cout << "   Interaction: " << i << " & " << j << endl;
	  countRows++;

	  //
	  //
	  //  ---Here's how all the "temps" stack up:---
	  //
	  //  distTempL * Lambdas_T * Lambdas * distTempL 
	  //                          \_________________/
	  //                               tempVec1
	  //              \_____________________________/
	  //                         tempVec2
	  //  \_________________________________________/
	  //        Global distance squared btwn 
	  //           sphere centers i and j
	  //
	  //

	  // calc tempVec1 (Lambdas*pTJJammer->distTempL)
	  // and 
	  for (int l=0; l<pTJJammer->dim; l++) {
	    tempVec1[l] = 0.0;
	    for (int d=0; d<pTJJammer->dim; d++) {
	      tempVec1[l] += gsl_matrix_get(pTJJammer->lambdas,l,d)*pTJJammer->distTempL[d];
	    }
	  }

	  // calc tempVec2 (Lambdas_T*tempVec1)
	  for (int l=0; l<pTJJammer->dim; l++) {
	    tempVec2[l] = 0.0;
	    for (int d=0; d<pTJJammer->dim; d++) {
	      tempVec2[l] += gsl_matrix_get(pTJJammer->lambdas,d,l)*tempVec1[d];
	    }
	  }

	  // input the epsilon coeffs
	  if (pTJJammer->strictJam == 0) {//Simple for collective jamming!
	    //#check
	    ar[countCons] = distSq; //Global distance squared between the spheres
	    ia[countCons] = countRows;
	    ja[countCons] = 1;
	    countCons++;
	  }
	  else if (pTJJammer->strictJam == 1) {
	    epsilonColTrack = 1;
	    for (int d=0; d<pTJJammer->dim; d++) {//First index of the epsilon
	      ar[countCons] = tempVec1[d]*tempVec1[d];
	      ia[countCons] = countRows;
	      ja[countCons] = epsilonColTrack;
	      countCons++;
	      epsilonColTrack = epsilonColTrack +1;
	      for (int l=d +1; l<pTJJammer->dim; l++) {//Second index of the epsilon
		ar[countCons] = 2*tempVec1[d]*tempVec1[l];
		ia[countCons] = countRows;
		ja[countCons] = epsilonColTrack;
		countCons++;
		epsilonColTrack = epsilonColTrack +1;
	      }
	    }
	  }



	  // input the pTJJammer->InternalDeltaR coeffs (order is: x1,x2,y1,y2,z1,z2...)
	  if (i!=j) {              // if i==j, then pTJJammer->InternalDeltaR coeffs are immaterial b/c coeff*(pTJJammer->InternalDeltaR1 - pTJJammer->InternalDeltaR1) = 0
	    for (int d=0; d<pTJJammer->dim; d++) {
	      ar[countCons] = tempVec2[d];
	      ia[countCons] = countRows;
	      ja[countCons] = epsilonNum + i*pTJJammer->dim +d +1;
	      countCons++;
	      ar[countCons] = -1.0*tempVec2[d];
	      ia[countCons] = countRows;
	      ja[countCons] = epsilonNum + j*pTJJammer->dim +d +1;
	      countCons++;
	    }
	  }
	  else {              // BUT WE STILL NEED TO HOLD THEIR PLACES FOR GUROBI!
	    for (int d=0; d<pTJJammer->dim; d++) {
	      ar[countCons] = 0.0;
	      ia[countCons] = countRows;
	      ja[countCons] = epsilonNum + i*pTJJammer->dim +d +1;
	      countCons++;
	      ar[countCons] = 0.0;
	      ia[countCons] = countRows;
	      ja[countCons] = epsilonNum + j*pTJJammer->dim +d +1;
	      countCons++;
	    }
	  }

	  // calc bounds (i.e., RHS of constraint)
	  coeffVal = ((pTJJammer->radii[i] + pTJJammer->radii[j])*(pTJJammer->radii[i] + pTJJammer->radii[j])) - distSq;    // note this is for additive diameters!
	  //for (int d=0; d<pTJJammer->dim; d++) {
	  //coeffVal -= pTJJammer->distTempL[d]*tempVec2[d];
	  //}
	  ineq  [countRows] = GLP_LO;
	  bounds[countRows] = coeffVal/2.0;

	  //Check:
	  //int lastIndex = countCons;
	  //cout << "  Constraint " << countRows+1 << ":" << endl;
	  //for (int checkIndex = startIndex ; checkIndex < lastIndex ; checkIndex++) {
	  //cout << "A(" << ia[checkIndex] << "," << ja[checkIndex] << ") = " << ar[checkIndex] << endl;
	  //}
	  //cout << "b(" << countRows << ") = " << bounds[countRows] << endl;
	  //startIndex = lastIndex;
	}
      }
    }
  }

  if(Verbosity>7)
  cout << "  Generated " << countRows << " constraints, " << countCons << " entries" << endl;

  return;
}


void LPClass::ApplyBulkStrainConstraint()
//
//  (For Strict Jamming only)
//
//The trace of the diagonal strain terms must be non-positive.
//NOTE: Since all of the constraints are uniformly input as
//"GREATER THAN OR EQUAL TO" inequalities, e.g.,
//a1x1 + a2x2 + ... >= b, 
//then the correct way to pose this is as
//-Tr(E) >= 0
//instead of 
// Tr(E) <= 0
//That's why ar[countCons] = -1.0, not +1.0
{
  //Add a term for each diagonal strain component:
  countRows++;
  ineq  [countRows] = GLP_LO;
  bounds[countRows] = 0.0;
  int epsilonIndex = 1;
  for (int i = 0 ; i < pTJJammer->dim ; i++) {
    ia[countCons] = countRows;
    ja[countCons] = epsilonIndex;
    ar[countCons] = -1.0;
    epsilonIndex += pTJJammer->dim - i;
    countCons++;
  }
  return;
}

//////////////////////////////////////////////////////////
//IMPORT THE PROBLEM TO THIRD-PARTY SOFTWARE


void LPClass::Import()
{
	if(Verbosity>7)
		cout << "  Import constraints into LP software" << endl;
	if (strcmp(pTJJammer->solver.c_str(),"GLPK") == 0) {
    Import_GLPK();
  }
  else if (strcmp(pTJJammer->solver.c_str(),"Gurobi") == 0) {
    Import_Gurobi();
  }
  else {
    cerr << "WARNING: Solver " << pTJJammer->solver << " is not programmed" << endl;
  }

  return;
}


void LPClass::Import_GLPK()
{
// create problem and read values in
  cout << "  Create new GLPK Object" << endl;
  lp = glp_create_prob();
  glp_set_prob_name(lp,"SLP");
  glp_set_obj_dir(lp,GLP_MIN);

  // specify control parameters: IMPORTANT!!!
  glp_init_smcp(&solverParams);
  solverParams.presolve = pTJJammer->GLPK_PresolveVar;
  solverParams.tol_bnd = pTJJammer->feasibleTol;
  solverParams.tol_dj = pTJJammer->feasibleTol;
  solverParams.msg_lev = GLP_MSG_ERR; //Only print error and warning messages to the terminal

  glp_get_bfcp(lp,&otherParams);
  otherParams.type = pTJJammer->GLPK_BasisFactType;
  glp_set_bfcp(lp, &otherParams);


  // first do the objective function and structural VARIABLES (i.e., position displacement and epsilons)
  glp_add_cols(lp,pTJJammer->dim*pTJJammer->N + epsilonNum); //if collective jamming, then epsilonNum=1

  //Deformable box variables:
  if (pTJJammer->strictJam==0) {
    glp_set_obj_coef(lp,1,1.0);
    glp_set_col_bnds(lp,1,GLP_DB,-1.0*pTJJammer->InternalCompMax,0.0); //Nowhere to expand in collective jamming!
  }
  else if (pTJJammer->strictJam==1) {
    epsilonColTrack = 1;
    for (int i=0; i<pTJJammer->dim; i++) {
      for (int l=i; l<pTJJammer->dim; l++) {
	if (l==i) {
	  glp_set_obj_coef(lp,epsilonColTrack,1.0);
	  glp_set_col_bnds(lp,epsilonColTrack,GLP_DB,-1.0*pTJJammer->InternalCompMax,pTJJammer->InternalCompMax);
	}
	else {
	  glp_set_col_bnds(lp,epsilonColTrack,GLP_DB,-1.0*pTJJammer->InternalShearMax,pTJJammer->InternalShearMax);
	}
	epsilonColTrack = epsilonColTrack +1;
      }
    }
  }

  //Sphere displacement variables
  for (int i=epsilonNum +1; i<=pTJJammer->dim*pTJJammer->N + epsilonNum; i++) {
    glp_set_col_bnds(lp,i,GLP_DB,-1.0*pTJJammer->InternalTransMax,1.0*pTJJammer->InternalTransMax);
  }

  // now do the constraints lower bounds
  glp_add_rows(lp,countRows);
  for (int i=1; i<=countRows; i++) {
    glp_set_row_bnds(lp,i,ineq[i],bounds[i],0.0);
  }

  return;
}


void LPClass::Import_Gurobi()
//Steps:
//1) Initialize the environment and run parameters
//2) Allocate variables and constraints
//3) Add strain variables
//4) Add pTJJammer->InternalDelta (movement) variables
//5) Add constraints

{
  //Part 1: Initialize the environment/set parameters
	try {
		if(Verbosity>7)
			cout << "  Initialize environment and model" << endl;
		gEnv   = new GRBEnv;
		//Set up the environment:
		SetupEnvironment_Gurobi();
		//Now intialize the model off of that environment:
		gModel = new GRBModel(*gEnv);
		gModel->set(GRB_IntAttr_ModelSense,1); //Set the problem to be a minimization problem

		//Part 2: Allocate variables (gVar) and constraints (gConstr)
		//cout << "  Allocate vars and constrs: " << epsilonNum + pTJJammer->dim*N << " variables and " << countRows << " constraints." << endl;
		gVar    = new GRBVar   [epsilonNum + pTJJammer->dim*pTJJammer->N];
		gConstr = new GRBConstr[countRows];
		gExpr   = new GRBLinExpr[countRows];

		//Part 3: Add strain variables
		//cout << "  Add strain vars" << endl;
		epsilonColTrack = 0;  //Start at zero because this is an array location (name is adjusted to be "+1")
		if (pTJJammer->strictJam==0) {
			AddEpsilon_Gurobi(0,0,epsilonColTrack);
			epsilonColTrack++;
		}
		else if (pTJJammer->strictJam==1) {

			for (int i=0; i<pTJJammer->dim; i++) {
				for (int l=i; l<pTJJammer->dim; l++) {
					AddEpsilon_Gurobi(i,l,epsilonColTrack);
					epsilonColTrack++;
				}
			}
		}

    //Part 4: add pTJJammer->InternalDelta (movement) variables
    //cout << "  Add pTJJammer->InternalDelta vars" << flush;
    int InternalDeltaIndex = 1;  //The ACTUAL NUMBER (for the name)!  (STARTS AT ONE)
    int varIndex = epsilonColTrack;
    for (int i=epsilonNum +1; i<=pTJJammer->dim*pTJJammer->N + epsilonNum; i++) {
      AddDelta_Gurobi(InternalDeltaIndex,varIndex);
      InternalDeltaIndex++;
      varIndex++;
    }
    //cout << "...done" << endl;

    //Update the model so that the variables are loaded in:
    gModel->update();

    //Part 5: Add the constraint equations:
    AddConstrs_Gurobi();

    //...update again?
    gModel->update();
  }
  catch (GRBException e) {
    cout << "    Error code   = " << e.getErrorCode() << endl;
    cout << "    Error Message: "     << e.getMessage()   << endl;
	throw;//has error in initialization, need to inform the caller of TJIterations()
  } catch (...) {
    cout << "  Unhandled error while initializing Gurobi " << endl;
	throw;//has error in initialization, need to inform the caller of TJIterations()
  }

  //cout<<"  finished import into Gurobi"<<endl;

  return;
}

void LPClass::SetupEnvironment_Gurobi()
{

	//Set feasibility limits:
	double minTol = 1.0e-9; //The limits on the FeasibilityTol and OptimalityTol parameters for Gurobi (Simplex)
	double maxTol = 1.0e-2;
	if (pTJJammer->feasibleTol >= minTol && pTJJammer->feasibleTol <= maxTol) {//Alowed by Gurobi:
		grbExtraScale = 1.0; //Import constraints as-is
		gEnv->set(GRB_DoubleParam_FeasibilityTol,pTJJammer->feasibleTol);
		gEnv->set(GRB_DoubleParam_OptimalityTol ,pTJJammer->feasibleTol);
	}
	else if (pTJJammer->feasibleTol < minTol) {//Lower than what is allowed by Gurobi!
		grbExtraScale = minTol / pTJJammer->feasibleTol;
		gEnv->set(GRB_DoubleParam_FeasibilityTol,minTol);
		gEnv->set(GRB_DoubleParam_OptimalityTol ,minTol);
		if(Verbosity>7)
			cout << "  Gurobi: Extra-tight Feasibility tolerances--grbExtraScale = " << grbExtraScale << endl;
	}
	else if (pTJJammer->feasibleTol > maxTol) {
		if(Verbosity>7)
			cout << "  Gurobi: Feasibility toleranes must be below 0.01" << endl;
		gEnv->set(GRB_DoubleParam_FeasibilityTol,maxTol);
		gEnv->set(GRB_DoubleParam_OptimalityTol ,maxTol);
	}
	gEnv->set(GRB_DoubleParam_BarConvTol,pTJJammer->feasibleTol); //Barrier convergence doesn't care :)

	//Choose which solver will be used (Default is "Auto"
	if (strcmp(pTJJammer->GRB_Method.c_str(),"Auto") == 0) {
		if(Verbosity>7)
			cout << "  Gurobi: Use Auto Solver\n";
		gEnv->set(GRB_IntParam_Method,-1);
	}
	else if (strcmp(pTJJammer->GRB_Method.c_str(),"PrimalSimplex") == 0) {
		if(Verbosity>7)
			cout << "  Gurobi: Use Primal Simplex\n";
		gEnv->set(GRB_IntParam_Method,0);
	}
	else if (strcmp(pTJJammer->GRB_Method.c_str(),"DualSimplex") == 0) {
		if(Verbosity>7)
			cout << "  Gurobi: Use Dual Simplex\n";
		gEnv->set(GRB_IntParam_Method,1);
	}
	else if (strcmp(pTJJammer->GRB_Method.c_str(),"Barrier") == 0) {
		if(Verbosity>7)
			cout << "  Gurobi: Use Barrier Solver\n";
		gEnv->set(GRB_IntParam_Method,2);
	}
	else if (strcmp(pTJJammer->GRB_Method.c_str(),"Concurrent") == 0) {
		if(Verbosity>7)
			cout << "  Gurobi: Use Concurrent Solver\n";
		gEnv->set(GRB_IntParam_Method,3); //Concurrent (non-deterministic!)
	}
	else if (strcmp(pTJJammer->GRB_Method.c_str(),"ConcurrentDeterm") == 0) {
		if(Verbosity>7)
			cout << "  Gurobi: Use Concurrent Deterministic Solver\n";
		gEnv->set(GRB_IntParam_Method,3); //Concurrent (deterministic)
	}
//Disable crossover:
  gEnv->set(GRB_IntParam_Crossover,0);

  //Set printing level:
  gEnv->set(GRB_IntParam_OutputFlag,0); //Disable output

  //Limit the number of threads that Gurobi uses:
  gEnv->set(GRB_IntParam_Threads,pTJJammer->grbThreads);

  //Presolve level:
  gEnv->set(GRB_IntParam_Presolve,pTJJammer->GRB_Presolve);   //(-1,0,1,2) = (Auto (Default) , none , standard , aggressive)

  return;
}


void LPClass::AddEpsilon_Gurobi(int i,
				int l,
				int varIndex)
{
  //Get the name for this new variable:
  stringstream ss;
  ss << i+1 << l+1; //Start at x1...
  string varName = string("e") + ss.str();

  //cout << "  Adding variable " << varName << " to the model" << endl;

  //Get the right bounds and objective coefficient:
  double objCoef   = 0.0;
  double lowBound  = 0.0;
  double highBound = 0.0;
  if (pTJJammer->strictJam==0) {
    objCoef = grbExtraScale*pTJJammer->dim;
    lowBound = -1.0*pTJJammer->InternalCompMax;
    highBound = 0.0; //No expansion allowed in Collective Jamming
  }
  else if (pTJJammer->strictJam==1) {
    if (i == l) {//Bulk strain (matrix trace terms)
      objCoef = grbExtraScale; //1.0 when LP is not pre-scaled
      lowBound = -1.0*pTJJammer->InternalCompMax;
      highBound = pTJJammer->InternalCompMax;
    }
    else {//Shear strain
      objCoef = 0.0;
      lowBound  = -1.0*pTJJammer->InternalShearMax;
      highBound =      pTJJammer->InternalShearMax;
    }
  }

  //C++ method:
  try {
    gVar[varIndex] = gModel->addVar(lowBound, //Lower bound
				    highBound,  //Upper bound
				    objCoef,  //Coefficient in objective fcn
				    GRB_CONTINUOUS,//Type of variable (real)
				    varName.c_str()); //Name of the variable (first variable is "x1", second is "x2," etc.
  }
  catch(GRBException e) {
    cout << "    Error code   = " << e.getErrorCode() << endl;
    cout << "    Error Message: "     << e.getMessage()   << endl;
	throw;//has error, inform the caller
  } catch (...) {
    cout << "  Unhandled error while adding variable " << varName << endl;
	throw;//has error, inform the caller
  }

  return;
}


void LPClass::AddDelta_Gurobi(int InternalDeltaIndex,
			      int varIndex)
{
  //Get the name for this new variable:
  stringstream ss;
  ss << InternalDeltaIndex;
  string varName = string("d") + ss.str();

  //cout << "  Adding variable " << varName << " to the model" << endl;

  //Get the right bounds and objective coefficient:
  double objCoef   = 0.0;
  double lowBound  = -1.0*pTJJammer->InternalTransMax;
  double highBound = pTJJammer->InternalTransMax;


  //C++ method:
  try {
    gVar[varIndex] = gModel->addVar(lowBound, //Lower bound
				    highBound,  //Upper bound
				    objCoef,  //Coefficient in objective fcn
				    GRB_CONTINUOUS,//Type of variable (real)
				    varName.c_str()); //Name of the variable (first variable is "x1", second is "x2," etc.
  }
  catch(GRBException e) {
    cout << "    Error code   = " << e.getErrorCode() << endl;
    cout << "    Error Message: "     << e.getMessage()   << endl;
	throw;//has error, inform the caller
  } catch (...) {
    cout << "  Unhandled error while adding variable " << varName << endl;
	throw;//has error, inform the caller
  }

  return;
}


void LPClass::AddConstrs_Gurobi()
{
  //cout << "  Gurobi: Add constraints..." << endl;
  //Allocate row indices array:
  int numTerms = epsilonNum + 2*pTJJammer->dim; //The number of variables associated with any given nonoverlap constraint
                                     //i.e. the strain terms, plus the movements for both particles

  //Loop through the constraints:
  int offset = 0;
  int constrIndex = 0;
  //The strain trace constraint (strict jamming only)
  if (pTJJammer->strictJam==1) {
    offset = AddConstr_Gurobi(constrIndex,pTJJammer->dim,offset); 
    constrIndex++;
  }

  //The interparticle contraints
  while (constrIndex < countRows) {
    //cout << "    Constraint " << constrIndex+1 << endl;
    offset = AddConstr_Gurobi(constrIndex,numTerms,offset);
    constrIndex++;
  }//End of constraint loop
  //cout << "done" << endl;

  return;
}


double LPClass::AddConstr_Gurobi(int constrIndex,
				 int numTerms,
				 int offset)
{
  //Get the name for this new constraint:
  stringstream ss;
  ss << constrIndex+1; //Start at C1...
  string constrName = string("C") + ss.str();
  //cout << "  Adding constraint " << constrName << endl;

  //Initialize the constraint expression:
  //GRBLinExpr expr;

  //Add the terms to the expression:
  //int offset = constrIndex * numTerms; //Where to start reading from ja and ar
  //cout << "  Building..." << endl;
  //cout << "    ia,ja,ar up to " << numTerms+offset << endl;
  for (int iTerm = 1 ; iTerm <= numTerms ; iTerm++) {//ia, ja, ar start at 1, not 0.
    int varIndex = ja[iTerm+offset]-1; //Which variable we're dealing with (MUST MOVE BACK 1!)
    if (varIndex >= (epsilonNum + pTJJammer->dim*pTJJammer->N)) {
      cerr << "Tried to call varIndex = " << varIndex << endl;
    }
    double coef  = ar[iTerm+offset];
    //cout << "  (row , entry , varIndex , coef) = (" << ia[iTerm+offset]-1 << " , " << iTerm+offset << " , " << varIndex << " , " << coef << ")" << endl;
    gExpr[constrIndex] += grbExtraScale * coef * gVar[varIndex]; //Wow, yeah, this is weird to be doing -__-
  }

  //figure out what comparator to use:
  char grbSign;
  if (ineq[constrIndex+1] == GLP_LO) {
    grbSign = GRB_GREATER_EQUAL;
  }
  else if (ineq[constrIndex+1] == GLP_FX) {
    grbSign = GRB_EQUAL;
  }
  else {
    grbSign = GRB_GREATER_EQUAL;
    cerr << "Error: Unspecified comparator on constraint " << constrIndex+1 << endl;
  }


  /*
  cout.precision(2);
  cout << "  Constraint " << constrName << ": " << expr[constrIndex] << endl;
  cout << grbSign << " " << bounds[constrIndex];
  cout.precision(COUT_PRECISION);
  */  

  //Add the constraint to the model:
  //cout << "  Adding..." << endl;

  try {
    gConstr[constrIndex] = gModel->addConstr(gExpr[constrIndex],                    //Left-hand side
					     grbSign,                               //Inequality symbol
					     grbExtraScale * bounds[constrIndex+1], //Right-hand side...ADD ONE FOR STUPID GLPK CONVENTION ARTIFACT!!!
					     constrName);                           //Name of the constraint (C1,C2,C3,etc...)
  }
  catch(GRBException e) {
    cout << "    Error code   = " << e.getErrorCode() << endl;
    cout << "    Error Message: " << e.getMessage()   << endl;
	throw;//has error, inform the caller
  } catch (...) {
    cout << "  Unhandled error while adding constraint " << constrName << endl;
	throw;//has error, inform the caller
  }

  //cout << "  Done." << endl;
  
  return(offset+numTerms);//So that the next constraint knows where to start
}


//////////////////////////////////////////////////////
//SOLVE THE LP


void LPClass::Solve()
{
	if(Verbosity>7)
		cout << "  Solve the LP" << endl;
	if (strcmp(pTJJammer->solver.c_str(),"GLPK") == 0) {
		Solve_GLPK();
	}
	else if (strcmp(pTJJammer->solver.c_str(),"Gurobi") == 0) {
		Solve_Gurobi();
	}
	else {
		cerr << "WARNING: Solver " << pTJJammer->solver << " is not programmed" << endl;
	}

	return;
}


void LPClass::Solve_GLPK()
{
  glp_load_matrix(lp,countCons-1,ia,ja,ar);

  // in case the basis matrix was singular on the last run
  if (warmUpParam) {
    warmUpReturn = glp_warm_up(lp);
    if (warmUpReturn == 0) {
      cout << "warm up successfully fixed the basis matrix\n";
    }
    else {
      cout << "warm up indicated that the basis matrix is numerically singular within the precision specified\n";
    }
  }

  solverStatus = glp_simplex(lp,&solverParams);

  if (solverStatus > 0) {
    cout << "solver encountered problems, running warm up to attempt to resolve\n";
    warmUpParam = true;
  }
  else {
    warmUpParam = false;
  }

  return;
}


void LPClass::Solve_Gurobi()
{
  //cout << "  GUROBI OPTIMIZE" << endl;
  try {
    gModel->optimize();
    
    //Figure out how the optimization went:
    int optimstatus = gModel->get(GRB_IntAttr_Status);

	if(optimstatus == 13 || optimstatus == 12)
	{
		//disabling crossover is causing precision issues, re-enable it and optimize again
		gModel->getEnv().set(GRB_IntParam_Crossover,1);
        gModel->optimize();
        optimstatus = gModel->get(GRB_IntAttr_Status);
	}

    if (optimstatus == GRB_INF_OR_UNBD) {
      gModel->getEnv().set(GRB_IntParam_Presolve, 0);
      gModel->optimize();
      optimstatus = gModel->get(GRB_IntAttr_Status);
    }


    if (optimstatus == GRB_OPTIMAL) {
      //cout << "OPTIMAL SOLUTION FOUND!" << endl;
      //double objval = gModel->get(GRB_DoubleAttr_ObjVal);
      //cout << "Optimal objective: " << objval << endl;
    } 
    else if (optimstatus == GRB_INFEASIBLE) {
		throw TJJammer::error("Model is infeasible");

      // compute and write out IIS

      //gModel->computeIIS();
      //gModel->write("model.ilp");
    } 
    else if (optimstatus == GRB_UNBOUNDED) {
		throw TJJammer::error("Model is unbounded");
	} 
	else {
		std::stringstream ss;
		ss << "Optimization was stopped with status = "
			<< optimstatus ;
		throw TJJammer::error(ss.str().c_str());
	}

  } catch(GRBException e) {
	  std::stringstream ss;
	  ss << "Gurobi optimization exception; Error code = " << e.getErrorCode() << endl;
	  ss << e.getMessage() << endl;
	  throw TJJammer::error(ss.str().c_str());
  } 
  catch(std::exception & a)
  {
	  std::cout<<"std::exception, what="<<a.what();
	  throw;
  }
  catch (...) 
  {
	  cout << "Error during optimization" << endl;
	  throw;
  }

  return;
}


//////////////////////////////////////////////////////
//EXPORT THE LP'S RESULTS


void LPClass::Export(int lpIters)
//Update the packing with the results of the LP solve
{
	if(Verbosity>7)
		cout << "  Export LP results..." << flush;
	//Part 1: Export Strain matrix
	if (strcmp(pTJJammer->solver.c_str(),"GLPK") == 0) {
		ExportStrain_GLPK();
	}
	else if (strcmp(pTJJammer->solver.c_str(),"Gurobi") == 0) {
		ExportStrain_Gurobi();
	}
	else {
		cerr << "WARNING: Solver " << pTJJammer->solver << " is not programmed" << endl;
	}

	//Part 2: Apply the strain:
	ApplyStrain();

	//Part 3: Apply the movements:
	ExportMovements();

	//Part 4: Look at and store the interesting quantities to file
	ReviewResults();

	//Part 5: Export the LP solution to file (if wanted)
	if (pTJJammer->printThisIter_LP) {
		if (strcmp(pTJJammer->solver.c_str(),"GLPK") == 0) {
			cerr<<"LP export not supported for GLPK\n";
		}
		else if (strcmp(pTJJammer->solver.c_str(),"Gurobi") == 0) {
			ExportLP_Gurobi();
		}
	}



	if(Verbosity>7)
		cout<<"done."<<endl;

	return;
}


void LPClass::ExportStrain_GLPK()
{
  // first the lattice vecs; start by creating the matrix (I + epsilon), which is updateLambdas
  if (pTJJammer->strictJam==0) {
    epsilonColTrack=1;
    double dVal = glp_get_col_prim(lp,epsilonColTrack);
    epsilonColTrack++;
    gsl_matrix_set(updateLambdas,0,0,1.0 + dVal); //Only need to store a scalar
  }
  else if (pTJJammer->strictJam==1) {
    epsilonColTrack = 1;
    for (int i=0; i<pTJJammer->dim; i++) {
      for (int d=i; d<pTJJammer->dim; d++) {
	double dVal = glp_get_col_prim(lp,epsilonColTrack);
	//cout << "   e" << i << d << " = " << dVal << endl;
	if (i==d) {
	  gsl_matrix_set(updateLambdas,i,d,1.0 + dVal);
	}
	else {
	  gsl_matrix_set(updateLambdas,i,d,dVal);
	  gsl_matrix_set(updateLambdas,d,i,dVal);
	}
	epsilonColTrack++;
      }
    }
  }

  return;
}


void LPClass::ExportStrain_Gurobi()
{
  // first the lattice vecs; start by creating the matrix (I + epsilon), which is updateLambdas
  if (pTJJammer->strictJam==0) {
    epsilonColTrack = 0;
    double dVal = gVar[epsilonColTrack].get(GRB_DoubleAttr_X);
    epsilonColTrack++;
    for (int d=0;d<pTJJammer->dim;d++)
      gsl_matrix_set(updateLambdas,d,d,1.0 + dVal); //Only need to store a scalar, but the matrix helps later 
  }
  else if (pTJJammer->strictJam==1) {
    epsilonColTrack = 0;
    for (int i=0; i<pTJJammer->dim; i++) {
      for (int d=i; d<pTJJammer->dim; d++) {
	double dVal = gVar[epsilonColTrack].get(GRB_DoubleAttr_X);
	//cout << "   e" << i << d << " = " << dVal << endl;
	if (i==d) {
	  gsl_matrix_set(updateLambdas,i,d,1.0 + dVal);
	}
	else {
	  gsl_matrix_set(updateLambdas,i,d,dVal);
	  gsl_matrix_set(updateLambdas,d,i,dVal);
	}
	epsilonColTrack++;
      }
    }
  }

  return;
}

void LPClass::ApplyStrain()
{
  // now multiply: pTJJammer->lambdas = (I + epsilon)*pTJJammer->lambdas
  if (pTJJammer->strictJam==0) {
    gsl_matrix_memcpy(updateTemp,pTJJammer->lambdas);
    gsl_matrix_scale(updateTemp,gsl_matrix_get(updateLambdas,0,0)); //SCalar operation!
  }
  else if (pTJJammer->strictJam==1) {
    for (int l=0; l<pTJJammer->dim; l++) {
      for (int i=0; i<pTJJammer->dim; i++) {
	updateTempVal = 0.0;
	for (int d=0; d<pTJJammer->dim; d++) {
	  updateTempVal += gsl_matrix_get(updateLambdas,l,d)*gsl_matrix_get(pTJJammer->lambdas,d,i);
	}
	gsl_matrix_set(updateTemp,l,i,updateTempVal);   // now updateTemp is the new pTJJammer->lambdas matrix, (I+epsilon)*pTJJammer->lambdas
      }
    }
  }

  return;
}


void LPClass::ExportMovements()
//Also keeps track of which particle moved farthest (relative to the total packing displacement)
//if termCriterion is "maxDisp"
{

  gsl_matrix *globalTotDisp = gsl_matrix_alloc(pTJJammer->dim,1); //The total movement of the entire packing (in global coords)
  gsl_matrix **globalDisp = new gsl_matrix*[pTJJammer->N]; //The movements (in global coords) of each individual sphere
  for (unsigned int d = 0 ; d < (unsigned int) pTJJammer->dim ; d++) {gsl_matrix_set(globalTotDisp,d,0,0.0);}
  for (int i = 0 ; i < pTJJammer->N ; i++) {
    globalDisp[i] = gsl_matrix_alloc(pTJJammer->dim,1);
  }

  // now update the coordinates
  greatestMaxMoveVal = 0.0;
  for (int i=0; i<pTJJammer->N; i++) {//Loop through the spheres
    for (int d=0; d<pTJJammer->dim; d++) {//Loop through pTJJammer->dimensions
      if      (pTJJammer->solver == "GLPK"  ) {
	InternalDeltaLocalCoords[d] = glp_get_col_prim(lp,epsilonNum + 1 + pTJJammer->dim*i + d);
      }
      else if (pTJJammer->solver == "Gurobi") {
	InternalDeltaLocalCoords[d] = gVar[epsilonNum + pTJJammer->dim*i + d].get(GRB_DoubleAttr_X);//Note: gVar starts from 0, not 1.
      }
      pTJJammer->localCoords[pTJJammer->dim*i +d] += InternalDeltaLocalCoords[d];
      if (pTJJammer->localCoords[pTJJammer->dim*i +d] > 1.0) {//Apply periodic boundary movement
	pTJJammer->localCoords[pTJJammer->dim*i +d] -= 1.0;
      }
      else if (pTJJammer->localCoords[pTJJammer->dim*i +d] < 0.0) {
	pTJJammer->localCoords[pTJJammer->dim*i +d] += 1.0;
      }
    }

    // now calc the maximum movement - only if NNLs are used
    if (pTJJammer->useNNL || pTJJammer->termCriterion == "maxDisp") {
      maxMoveVal = 0.0;
      for (int l=0; l<pTJJammer->dim; l++) {
	oldLambdasNewLocals[l] = 0.0;
	for (int d=0; d<pTJJammer->dim; d++) {
	  oldLambdasNewLocals[l] += gsl_matrix_get(pTJJammer->lambdas,l,d)*pTJJammer->localCoords[pTJJammer->dim*i +d];
	}
      }
      for (int l=0; l<pTJJammer->dim; l++) {
	maxMoveVals[l] = 0.0;
	for (int d=0; d<pTJJammer->dim; d++) {
	  if (l==d) {
	    epsilonVal = gsl_matrix_get(updateLambdas,l,d) - 1.0;
	  }
	  else {
	    epsilonVal = gsl_matrix_get(updateLambdas,l,d);
	  }
	  maxMoveVals[l] += epsilonVal*oldLambdasNewLocals[d];
	}
	maxMoveVals[l] += InternalDeltaLocalCoords[l];
	if (pTJJammer->termCriterion == "maxDisp") {
	  gsl_matrix_set(globalDisp[i],l,0,maxMoveVals[l]);
	}
	maxMoveVal += maxMoveVals[l]*maxMoveVals[l];
      }
      if (maxMoveVal > greatestMaxMoveVal) {
	greatestMaxMoveVal = maxMoveVal;
      }
    }
  }//Sphere loop

  //If the termination criterion is "maxDisp", do those calculations now:
  if (pTJJammer->termCriterion == "maxDisp") {
    GetMaxDisp(globalTotDisp,globalDisp);
  }

  if (pTJJammer->useNNL) {
    greatestMaxMoveVal = sqrt(greatestMaxMoveVal);       // the scaling is necessary to match InternalNNLextraDist
    pTJJammer->maxSingleMove += greatestMaxMoveVal;
	if(Verbosity>7)
		cout << "  greatestMaxMoveVal = " << greatestMaxMoveVal << endl;
	if(Verbosity>7)
		cout << "  maxSingleMove      = " << pTJJammer->maxSingleMove << "\n";
  }

  // now update the pTJJammer->lambdas matrix officially
  gsl_matrix_memcpy(pTJJammer->lambdas,updateTemp);

  latVol = pTJJammer->getVol();
  if(Verbosity>7)
	  cout << "  lattice volume     = " << latVol << "\n";

  // now check for overlap
  if (pTJJammer->lpIters == pTJJammer->maxIters -1) {                 // if it's the last iteration, resize only so greatest overlap is at contact!
	  ifResize = pTJJammer->resizeIfOverlap(false);
  }
  else {
	  ifResize = pTJJammer->resizeIfOverlap(true);
  }

  if (ifResize) {
	  //cout << "overlap detected, resize necessary\n";
	  latVol = pTJJammer->getVol();
	  if(Verbosity>7)
		  cout << "  new lattice volume is " << latVol << "\n";
  }

  //Deallocate:
  gsl_matrix_free(globalTotDisp);
  for (int i = 0 ; i < pTJJammer->N ; i++) {
    gsl_matrix_free(globalDisp[i]);
  }
  delete [] globalDisp;

  return;
}


void LPClass::GetMaxDisp(gsl_matrix  *globalTotDisp,
			 gsl_matrix **globalDisp)
//For the "maxDisp" termination criterion
{
  pTJJammer->maxDisp = 0.0; //Reset

  //First, get the global movement:
  for (int i = 0 ; i < pTJJammer->N ; i++) {
    gsl_matrix_add(globalTotDisp,globalDisp[i]);
  }
  gsl_matrix_scale(globalTotDisp,1.0 / ((double) pTJJammer->N));
  for (int i = 0 ; i < pTJJammer->N ; i++) {
    gsl_matrix_sub(globalDisp[i],globalTotDisp); //Remove packing's motion
    //Find the resulting movement:
    double thisDisp = 0.0;
    for (unsigned int d = 0 ; d < (unsigned int) pTJJammer->dim ; d++) {
      double tempVal = gsl_matrix_get(globalDisp[i],d,0);
      thisDisp += tempVal*tempVal;
    }
    if (thisDisp > pTJJammer->maxDisp) {//Compare the displacements squared, since this is faster
      pTJJammer->maxDisp = thisDisp;
    }
  }

  pTJJammer->maxDisp = sqrt(pTJJammer->maxDisp) / pTJJammer->biggestRad; //Scale back to true displacement (L2-norm), non-pTJJammer->dimensionalize in terms of biggest radius

  cout << "  maxDisp = " << pTJJammer->maxDisp << endl;

  return;
}


void LPClass::ReviewResults()
//Check the values that are being output as the solution
//  Just in case anything interesting is going on...
{
	if(Verbosity>7)
		cout << "  Check optimization results" << endl;

	//Check the strain matrix:
	double trace = 0.0;
	for (int d = 0 ; d < pTJJammer->dim ; d++) {
		trace += gsl_matrix_get(updateLambdas,d,d);
	}
	trace -= pTJJammer->dim;
	if(Verbosity>7)
		cout << "    Tr(strain) = " << trace << endl;
	if (trace == 0) {
		if(Verbosity>7)
			cout << "    Strain trace constraint ACTIVE" << endl;
	}

	//Check individual strain elements:
	for (int i = 0 ; i < pTJJammer->dim ; i++) {
		//Diagonals:
		if (gsl_matrix_get(updateLambdas,i,i) - 1.0 == pTJJammer->InternalCompMax) {
			if(Verbosity>7)
				cout << "  strain ("<<i+1<<","<<i+1<<") hit InternalCompMax (limited expansion)" << endl;
		}
		if (gsl_matrix_get(updateLambdas,i,i) - 1.0 == -pTJJammer->InternalCompMax) {
			if(Verbosity>7)
				cout << "  strain ("<<i+1<<","<<i+1<<") hit InternalCompMax (limited compression)" << endl;
		}
		if (pTJJammer->strictJam == 1) {//Only if it matters...
			for (int j = i+1 ; j < pTJJammer->dim ; j++) {//Off-diagonals
				if (fabs(gsl_matrix_get(updateLambdas,i,j)) == pTJJammer->InternalShearMax) {
					if(Verbosity>7)
						cout << "  strain ("<<i+1<<","<<j+1<<") hit InternalShearMax (limited shear)" << endl;
				}
			}
		}
	}

  //Check sphere translations:
  int numHitTransMax = 0;
  for (int i=0; i<pTJJammer->N; i++) {//Loop through the spheres
    for (int d=0; d<pTJJammer->dim; d++) {//Loop through pTJJammer->dimensions
      double dVal=0.0;
      if      (pTJJammer->solver == "GLPK"  ) {
	dVal = glp_get_col_prim(lp,epsilonNum +1 +pTJJammer->dim*i +d);
      }
      else if (pTJJammer->solver == "Gurobi") {
	dVal = gVar[epsilonNum +pTJJammer->dim*i +d].get(GRB_DoubleAttr_X);//Note: gVar starts from 0, not 1.
      }
      if (fabs(dVal) == pTJJammer->InternalTransMax) {
	numHitTransMax++;
      }
	}
  }
  if (numHitTransMax > 0) {
	  if(Verbosity>7)
		  cout << "    "<<numHitTransMax<<" translations hit InternalTransMax" << endl;
  }

  UpdateLPFile(pTJJammer->lpIters,
	       trace,
	       (double) numHitTransMax / (double) pTJJammer->dim);
	       

  return;
}



void LPClass::ExportLP_Gurobi()
//For now, the LP export file has this format:
//
//numConstr
//(loop through the constrs)
//[sphere 1] [sphere 2] ["1" for active, "0" for not]
//
//
{
  cout<<"    Exporting LP\n";

  stringstream ss; ss<<pTJJammer->lpIters+1;
  string filename = pTJJammer->parentDirectory + pTJJammer->outputFilename + "_" + ss.str() + "LP.dat";
  try {
    ofstream outFile(filename.c_str() , ios::out);
    outFile<<gModel->get(GRB_IntAttr_NumConstrs)<<"\n";
    if (outFile.is_open()) {
      int varCount=1; //To go through the variables in the ia,ja,ar arrays
      for (int i = 0 ; i < (gModel->get(GRB_IntAttr_NumConstrs)) ; i++) {//Loop through the constraints
	//Output the spheres involved, and whether the constraint is active
	if (gExpr[i].size() == (unsigned int)(pTJJammer->dim)) {//Trace constraint
	  outFile<<"-1\t-1\t"<<gConstr[i].get(GRB_IntAttr_CBasis)<<"\n";
	}
	else if (gExpr[i].size() == (unsigned int)(epsilonNum + 2*pTJJammer->dim)) {//A constraint between a pair of spheres
	  const int sphere1 = (ja[varCount+epsilonNum  ] - epsilonNum - 1) / pTJJammer->dim;
	  const int sphere2 = (ja[varCount+epsilonNum+1] - epsilonNum - 1) / pTJJammer->dim;
	  //cout<<"Pair: "<<sphere1<<" , "<<sphere2<<"\n";
	  //cout<<"Check: "<<gExpr[i]<<"\n";
	  outFile<<sphere1<<"\t"<<sphere2<<"\t"<<gConstr[i].get(GRB_IntAttr_CBasis)<<"\n";
	}
	else {
	  cout<<"Constraint "<<i<<" has "<<gExpr[i].size()<<" terms\n";
	}
	varCount += gExpr[i].size();
      }
      outFile.close(); //and close
    }
  }
  catch (...) {
    cerr << "    Error trying to write " << filename << endl;
	throw;//has error, inform the caller
  }

  cout << "    Wrote to " << filename << endl;

  return;
}

  
void LPClass::UpdateLPFile(int    lpIters, double strainTrace, double freeSpheres)
//Save to a file all of the statistics
//  of each optimization iteration that
//  are interesting
{
	if(Verbosity>5)
	{
		if(Verbosity>7)
			cout << "  Write optimization stats to file" << endl;
		//Open the file (to make or to append):
		string filename = pTJJammer->parentDirectory + pTJJammer->outputFilename + "_Stats.dat";
		try {
			ofstream outFile;
			bool haveFile = false;
			if (lpIters > 0) {
				try {
					haveFile = boost::filesystem::exists(filename);
				}
				catch (...) 
				{
					std::cerr<<"Error using boost::filesystem!\n";
					throw;
				}
			}
			if (haveFile) {
				outFile.open(filename.c_str() , ios::app);
			}
			else {//Create the file and add the header:
				cout << "    Create a new file" << endl;
				outFile.open(filename.c_str() , ios::out);
				if (!outFile.is_open()) {
					cerr << "    Failed to open " << filename << endl;
					throw;
				}
				outFile << "#Iteration\t";
				outFile << "#phi\t";
				outFile << "#Tr(strain)\t";
				outFile << "#FreeSpheres\t";
				outFile << "#Timer\t";
				outFile << endl;
			}
			//Proceed to write this iteration's data:
			pTJJammer->timer.End();
			outFile.precision(pTJJammer->printPrecision);
			outFile << lpIters << "\t";
			outFile << pTJJammer->getSphereVol() / pTJJammer->getVol() << "\t";
			outFile << strainTrace << "\t";
			outFile << freeSpheres << "\t";
			outFile << pTJJammer->timer.timePassed << "\t";

			outFile << endl; //End the line
			outFile.close(); //and close
		}
		catch (...) {
			cerr << "    Error trying to write Optimization stats" << endl;
			throw;//has error, inform the caller
		}

		if(Verbosity>7)
			cout << "    Wrote to " << filename << endl;
	}
	return;
}

//////////////////////////////////////////////////////
//CLEANUP FOR LPCLASS
void LPClass::Cleanup()
//Deallocate arrays, etc
{
  // free memory
  if (strcmp(pTJJammer->solver.c_str(),"GLPK") == 0) {
    cout << "  Deallocate GLPK Object" << endl;
    glp_delete_prob(lp);
    //glp_free_env();
  }
  else if (strcmp(pTJJammer->solver.c_str(),"Gurobi") == 0) {
    delete gModel;
    delete gEnv;
    delete [] gVar;
    delete [] gConstr;
    delete [] gExpr;
  }

  delete [] ia;
  delete [] ja;
  delete [] ar;
  delete [] ineq;
  delete [] bounds;
  delete [] tempVec1;
  delete [] tempVec2;
  delete [] maxMoveVals;
  delete [] oldLambdasNewLocals;
  delete [] InternalDeltaLocalCoords;
  gsl_matrix_free(updateLambdas);
  gsl_matrix_free(updateTemp);

  return;
}



TJJammer::TJJammer()
{
	// random number generator initialization
	gsl_rng_env_setup();
	RNGtype = gsl_rng_mt19937;
	RNG = gsl_rng_alloc(RNGtype);
	time_t randomTime = time(nullptr);
	gsl_rng_set(RNG,randomTime);
	parentDirectory="./";
	outputFilename="_TJOutput.txt";
	tempFolder="";
	useNNL=true;
	autoSave=false;
	N=0;
	dim=0;
	maxInitRad=false;
	transMax=compMax=shearMax=0.01;
	strictJam=1;
	delta=0.45;
	NNLextraDist=1.0;
	maxIters=100000;
	termCriterion="latticeVol";
	phiTerm=1.0;
	termTol=2e-10;
	MCStartSweeps=0;
	MCMidSweeps=0;
	mc_dispLimit=0;
	randLatVecs=false;
	maxLatLength=1.0;
	minLatLength=0.7;
	minAngDeg=70;
	manualMinVol=0.6;
	radLatTol=0.75;
	topFail=10000;
	maxVolNum=10000;
	inputPhi=0.1;
	overBoxes=0;
	maxBoxNum=0;
	randOverlapMove=false;
	resizeSpace=1e-6;
	resizeTol=1e-12;
	printEvery=0;
	printThisIter_LP=false;
	printThisIter_pkg=false;
	printPrecision=15;
	useTimer=false;
	scaleUnitRadius=false;
	rescaleIC=1.0;
	printLPEvery=0;
	solver="Gurobi";
	localCoords=nullptr;
	globalCoords=nullptr;
	radii=nullptr;
	distTempG=nullptr;
	distTempL=nullptr;
	neighborListDelta=nullptr;
	neighborListOverlap=nullptr;
	lpIters=0;
	maxSingleMove=0.0;
	biggestRad=0.0;
	feasibleTol=1e-16;
	pivotTol=0.1;
	GLPK_PresolveVar=1;
	GLPK_BasisFactType=1;
	GRB_Method="Barrier";
	GRB_Presolve=0;
	grbThreads=1;
}

SpherePacking TJJammer::GetPacking()
{
	return result;
}
void TJJammer::SetPacking(const SpherePacking & p)
{
	this->dim=p.GetDimension();
	this->N=p.NumParticle();
	// initialize the variable coordinates - necessary regardless of whether an initial config (readInitFile = true) is provided
	lambdas = gsl_matrix_calloc(dim,dim);
	localCoords = new double[N*dim];
	globalCoords = new double[N*dim];
	radii = new double[N];
	distTempL = new double[dim];
	distTempG = new double[dim];

	for (int i=0; i<N; i++) {
		radii[i] = 0.0;
		for (int d=0; d<dim; d++) {
			localCoords[dim*i +d] = 0.0;
			globalCoords[dim*i +d]= 0.0;
		}
	}

	for (int i=0; i<dim; i++) {
		distTempL[i] = 0.0;
		distTempG[i] = 0.0;
	}

	// read in the radii configuration file if necessary
	{
		for (int i=0; i<N; i++) 
		{
			radii[i]=p.GetCharacteristics(i);
		}
	}
	// lattice vectors next
	for (int i=0; i<dim; i++) {
		for (int j=0; j<dim; j++) {
			double val=p.GetBasisVector(i).x[j];
			gsl_matrix_set(lambdas,j,i,val);
		}
	}

	// then coordinates and radii
	for (int i=0; i<N; i++) {
		for (int d=0; d<dim; d++) {
			globalCoords[dim*i +d]=p.GetCartesianCoordinates(i).x[d];
		}
	}

	// calculate inverse
	inverseLambdas = getInverse(lambdas);

	// calculate local coordinates
	for (int i=0; i<N; i++) {
		for (int j=0; j<dim; j++) {
			for (int d=0; d<dim; d++) {
				double val = gsl_matrix_get(inverseLambdas,j,d);
				localCoords[dim*i +j] = localCoords[dim*i +j] + val*globalCoords[dim*i +d];
			}
		}
	}

	for (int i=0; i<N; i++) {
		if (radii[i] > biggestRad) {
			biggestRad = radii[i];
		}
	}
}
void TJJammer::MC_Init()
{
  mc_dispLimit = radii[0]; //In local coords

  return;
}

double TJJammer::MC_Sweep(int numTries)
//Try N MC moves
//keep track of the success rate
{
  cout << "MC_Sweep()\n";
  inverseLambdas = getInverse(lambdas);
  int numSuccess = 0;
  double mc_maxDisp = 0.0;
  for (int i = 0 ; i < numTries ; i++) {
    //Pick a random sphere:
    int thisSphere = floor(N*gsl_rng_uniform_pos(RNG));
    //cout << "try sphere " << thisSphere << endl;
    double thisDisp = MC_Move(thisSphere);
    if (thisDisp > 0.0) {
      numSuccess++;
      if (thisDisp > mc_maxDisp) {
	mc_maxDisp = thisDisp;
      }
    }
  }

  double mc_accRate = (double) numSuccess / (double) numTries;
  //Print results:
  cout << "  MC Sweep Summary:\n";
  cout << "    mc_accRate   = " << numSuccess << " / " << numTries << " = " << mc_accRate << "\n";
  cout << "    mc_dispLimit = " << mc_dispLimit << endl;
  //Rescale maximum move:
  MC_AdjustDispLimit(mc_accRate);

  return(mc_maxDisp);
}


double TJJammer::MC_Move(int i)
//Try to move particle i
//Return resulting displacement (0 if failure)
{
  double *posG = new double[dim];
  double *oldPosL = new double[dim];
  double thisMoveDist = 0.0;
  //cout << "    Move:\n";
  for (int d = 0 ; d < dim ; d++) {
    //Store old position:
    oldPosL[d] =  localCoords[i*dim + d];
  }

  for (int j = 0 ; j < dim ; j++) {
    //Find global:
    posG[j] = 0.0;
    for (int d = 0 ; d < dim ; d++) {
      posG[j] += gsl_matrix_get(lambdas,j,d) * oldPosL[d];
    }
    //Move to a new position:
    double thisMove = mc_dispLimit * (2.0*gsl_rng_uniform_pos(RNG) - 1.0);
    //cout << "      " << posG[j] << " + " << thisMove << endl;
    posG[j] += thisMove;
    thisMoveDist += thisMove*thisMove;
  }
  thisMoveDist = sqrt(thisMoveDist);

  //Find new local position:
  for (int j = 0 ; j < dim ; j++) {
    int k = i*dim+j; //We'll use this a few times in this loop:
    localCoords[k] = 0.0;
    for (int d = 0 ; d < dim ; d++) {
      localCoords[k] += gsl_matrix_get(inverseLambdas,j,d) * posG[d];
    }
    //Put back in the FC:
    if (localCoords[k] > 1.0)
      localCoords[k] -= 1.0;
    else if (localCoords[k] < 0.0)
      localCoords[k] += 1.0;
  }

  //Check for overlaps (crude for now because this shouldn't take long compared to the LP)
  bool foundOverlap = false;
  if (overBoxes > 0) {
    int indexNum = int(pow(2.0*double(overBoxes) +1.0,dim));
    int selfIndexSkip = (indexNum -1)/2;              // this number is the self-image in the same unit cell
    for (int j = 0 ; j < N ; j++) {
      for (int k=0; k < indexNum; k++) {
	if (j==i && k == selfIndexSkip) {    // skip self image in same box
	  continue;
	}
	int quotient = k;
	for (int d=0; d<dim; d++) {
	  int remain = quotient % (2*overBoxes +1);
	  double shift = double(remain - overBoxes);        // this is the image that is shift lattice cells over
	  distTempL[d] = localCoords[dim*i + d] - localCoords[dim*j + d] + shift;
	  quotient = quotient/(2*overBoxes +1);
	}

	double dist = getGlobalLength(distTempL);

	if (dist < radii[i] + radii[j]) {             // note that this is for additive diameters!
	  foundOverlap = true;
	  //cout << "    Overlap with " << j << endl;
	  break;
	}
      }
    }
  }
  else { //No overBoxes -- use L/2 method
    for (int j = 0 ; j < N ; j++) {
      if (i == j) continue; //Don't check self-overlap with nearest image!
      for (int d=0; d<dim; d++) {
	distTempL[d] = localCoords[dim*i + d] - localCoords[dim*j + d];
	if (distTempL[d] > 0.5) {
	  distTempL[d] = distTempL[d] - 1.0;
	}
	else if (distTempL[d] < -0.5) {
	  distTempL[d] = distTempL[d] + 1.0;
	}
      }

      double dist = getGlobalLength(distTempL);

      if (dist < radii[i] + radii[j]) {             // note that this is for additive diameters!
	foundOverlap = true;
	//cout << "    Overlap with " << j << endl;
	break;
      }
    }
  }

  if (foundOverlap) {//Reset position:
    for (int d = 0 ; d < dim ; d++) {
      localCoords [i*dim+d] = oldPosL[d];
      thisMoveDist = 0.0;
    }
  }

  //Deallocate:
  delete [] posG;
  delete [] oldPosL;
  
  return(thisMoveDist);
}


void TJJammer::MC_AdjustDispLimit(double accRate)
{
  double lambda = 0.5; //Relaxation factor (to prevent large overshoot)
  mc_dispLimit *= (1.0 - lambda) + lambda / (2.0 * (1.0 - accRate));
  if (mc_dispLimit > 2.0*radii[0]) {//put a limit on it (kinda arbitrary)
    mc_dispLimit = 2.0*radii[0];
  }
  else if (mc_dispLimit < 1.0e-12) {
    mc_dispLimit = 1.0e-12;
  }

  return;
}

void TJJammer::RescalePacking()
{
  //Rescale lattice matrix:
  gsl_matrix_scale(lambdas,rescaleIC);
  inverseLambdas = getInverse(lambdas); //Re-compute
  //Rescale coordinates:
  int k=0;
  for (int i = 0 ; i < N ; i++) {
    for (int d = 0 ; d < dim ; d++) {
      globalCoords[k] *= rescaleIC;
      k++;
    }
  }
  //Local coordinates remain unchanged.
  return;
}
void TJJammer::DeleteOldNNL()
{
  for (int i=0; i<N; i++) {
    if (neighborListDelta[i] != 0) {
      delete neighborListDelta[i];
      neighborListDelta[i] = 0;
    }
    if (neighborListOverlap[i] != 0) {
      delete neighborListOverlap[i];
      neighborListOverlap[i] = 0;
    }
  }
}
int TJJammer::calcNNLs() 
{
  // variables
  NeighborNode *addressTempNNL = 0;
  double dist = 0.0;
  int totNeighbors = 0;    // total number of neighbor contacts (e.g., if 6 neighbors per particle, this is 3*N)

  // ensure deletion
  DeleteOldNNL();

  // calculate neighbor lists
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      // eliminate self-neighbor
      if (i == j) {
	continue;
      }

      // use the L/2 method to calculate distance
      for (int d=0; d<dim; d++) {
	distTempL[d] = localCoords[dim*i + d] - localCoords[dim*j + d];
	if (distTempL[d] > 0.5) {
	  distTempL[d] = distTempL[d] - 1.0;
	}
	else if (distTempL[d] < -0.5) {
	  distTempL[d] = distTempL[d] + 1.0;
	}
      }

      dist = getGlobalLength(distTempL);

      // create the two neighbor lists
      if (dist <= InternalDelta + radii[i] + radii[j] + InternalNNLextraDist) {
	// the first distance check is for the inclusion sphere
	totNeighbors = totNeighbors +1;
	// add a node at the beginning of the list
	addressTempNNL = neighborListDelta[i];
	neighborListDelta[i] = new NeighborNode;
	neighborListDelta[i]->next = addressTempNNL;
	neighborListDelta[i]->index = j;

	// now the neighborList for overlap calc
	if (dist <= radii[i] + radii[j] + InternalNNLextraDist) {
	  addressTempNNL = neighborListOverlap[i];
	  neighborListOverlap[i] = new NeighborNode;
	  neighborListOverlap[i]->next = addressTempNNL;
	  neighborListOverlap[i]->index = j;
	}
      }
    }
  }

  // returns
  return totNeighbors/2;
}
void TJJammer::MaximizeRadii()
//I know it's not fast, but it's not the real timesink in this program...
{
  cout<<"Maximize radii\n";

  double maxExpansion = 1.0; //Multiply the radii by this number.

  //Find the maximum allowed expansion
  if (useNNL) {
    NeighborNode *addressTempNNL = 0;
    for (int i=0; i<N; i++) {
      addressTempNNL = neighborListDelta[i];
      while (addressTempNNL) {
	int j = addressTempNNL->index;
	for (int d=0; d<dim; d++) {
	  distTempL[d] = localCoords[dim*i + d] - localCoords[dim*j + d];
	  if      (distTempL[d] >  0.5) distTempL[d] -= 1.0;
	  else if (distTempL[d] < -0.5) distTempL[d] += 1.0;
	}
	double thisExp = (getGlobalLength(distTempL) - resizeSpace) / (radii[i]+radii[j]);
	if (thisExp > maxExpansion) maxExpansion = thisExp;
      }
    }
  }
  else {//No NNL
    // box variables
    int indexNum = int(pow(2.0*double(overBoxes) +1.0,dim));
    int selfIndexSkip = (indexNum -1)/2;             // this number is the self-image in the same unit cell
    double shift = 0.0;
    int quotient = 0;
    int remain = 0;
    for (int i=0; i<N; i++) {
      for (int j=i; j<N; j++) {
	for (int k=0; k < indexNum; k++) {
	  if (j==i && k == selfIndexSkip) {    // skip self image in same box
	    continue;
	  }
	  quotient = k;
	  for (int d=0; d<dim; d++) {
	    remain = quotient % (2*overBoxes +1);
	    shift = double(remain - overBoxes);        // this is the image that is shift lattice cells over
	    distTempL[d] = localCoords[dim*i + d] - localCoords[dim*j + d] + shift;
	    quotient = quotient/(2*overBoxes +1);
	  }
	  double thisExp = (getGlobalLength(distTempL) - resizeSpace) / (radii[i]+radii[j]);
	  if (thisExp > maxExpansion) maxExpansion = thisExp;
	}
      }
    }
  }

  //Perform expansion:
  for (int i = 0 ; i < N ; i++)
    radii[i] *= maxExpansion;

  return;
}

bool TJJammer::CheckTermination(double lastLastVol, double latVol, double phi)
{
	if (termCriterion == "latticeVol") {
		if (fabs(lastLastVol - latVol) <= termTol) {
			if(Verbosity>5)
				cout << "Terminated due to failure to decrease volume over 2 runs\n";
			return(true); //Yes, terminate the LP iterations.
		}
	}
	else if (termCriterion == "maxDisp") {
		if (maxDisp < termTol) {
			if(Verbosity>5)
				cout << "Terminated due to particle movements being too small\n";
			return(true);
		}
	}
	if (phi >= phiTerm) {
		if(Verbosity>5)
			cout << "Terminated due to maximum density being met\n";
		return(true);
	}

	return(false); //No termination criteria were met
}

int TJJammer::printPacking(string thisOutputFilename, bool addIterNum)
{
	if(Verbosity>5)
		cout << "Printing packing..." << endl;
	// variables
	ofstream outFile;
	int success = 0;
	double globalVal = 0.0;
	//string thisOutputFilename = parentDirectory + outputFilename;
	stringstream intConvert;

	if (addIterNum) {
		intConvert << lpIters;
		thisOutputFilename.append("_");
		thisOutputFilename.append(intConvert.str());
	}

	// open file
	thisOutputFilename.append(".dat");
	if(Verbosity>5)
		cout << "  print location = " << thisOutputFilename << endl;
	outFile.open(thisOutputFilename.c_str(),ios::out);

	if (!outFile.is_open()) {
		cerr << "Can't create output file " << thisOutputFilename << endl;
		return success;
	}

	//Test for monodisperse -vs- polydisperse packing:
	double radVal = radii[0];
	bool monoPacking = true;
	for (int i = 1 ; i < N ; i++) {if (radii[i] != radVal) {monoPacking = false; break;}}

	// set the output precision
	outFile.precision(printPrecision);

	// print the packing (Alek's format)
	if (monoPacking) {
		if(Verbosity>5)
			cout << "Printing monodisperse packing format" << endl;
		outFile << dim << "\t HS\t mono" << endl;
	}
	else {
		outFile << dim << "\t HS\t poly" << endl;
	}
	outFile << N << "\t" << 1 << endl;
	outFile << N << endl;
	if (monoPacking) {outFile << 2.0 * radVal << endl;}
	//Lattice matrix:
	for (int i=0;i<dim;i++) {
		for (int j=0; j<dim; j++) {
			outFile << gsl_matrix_get(lambdas,j,i) << "\t";
		}
		outFile << endl;
	}
	for (int i=0; i<dim; i++) {
		outFile << "T\t";
	}
	outFile << endl;
	//Sphere loop:
	for (int i=0; i<N; i++) {
		for (int d=0; d<dim; d++) {
			globalVal = 0.0;
			for (int l=0; l<dim; l++) {
				globalVal = globalVal + gsl_matrix_get(lambdas,d,l)*localCoords[dim*i +l];    // calculate Lambda*localVec = global vec
			}
			outFile << globalVal << "\t";
		}
		if (!monoPacking) {outFile << radii[i];}
		outFile << endl;
	}

	outFile.close();

	// returns
	success = 1;
	return success;
}
void TJJammer::SaveProgress()
//How to check for existence of a directory
//How to make a directory
{
  cout << "Saving progress..." << endl;

  //Figure out if there's a TEMP Directory in the current parent directory
  bool haveFolder = false;
  try {
    haveFolder = boost::filesystem::is_directory(parentDirectory + tempFolder);
  }
  catch(...) {
    cout << "Directory doesn't exist" << endl;
    haveFolder = false;
  }
  //bool haveFolder = true;

  if (!haveFolder) {
    cout << "Making directory " << parentDirectory << tempFolder << endl;
    if (boost::filesystem::create_directory(parentDirectory + tempFolder) == 0) {
      cerr << "ERROR CREATING DIRECTORY FOR SAVE INFO" << endl;
    }
  }
  else {
    cout << "Temporary save data directory already exists" << endl;
  }

  //Write the file there:
  printPacking(parentDirectory + tempFolder + "SavePacking" , false); //appends ".dat" automatically

  //Write save log:
  cout << "Writing save log..." << endl;
  ofstream outFile;
  string saveLogName = parentDirectory + tempFolder + "SaveInfo.txt";
  outFile.open(saveLogName.c_str(), ios::out);
  if (!outFile.is_open()) {
    cerr << "Unable to open Save info file at " << saveLogName << endl;
  }
  else {
    outFile << "//Automatically-saved data\n";
    outFile << "#lastIteration\t" << lpIters << endl;
    outFile.close();
  }

  return;
}



void TJJammer::ClearSaveData()
//Clear out the temporary save information since we made it to the end of the program!
{
  cout << "Clearing temporary save data" << endl;

  if (tempFolder.empty()) {//No directory?  no clear!
    cout << "tempFolder is an empty string; don't try to delete.\n";
    return;
  }

  try {
    if (boost::filesystem::is_directory(parentDirectory + tempFolder)) {
      if (boost::filesystem::remove_all(parentDirectory + tempFolder) == 0) {
	cerr << "Error deleting temporary save data" << endl;
      }
    }
    else {
      cout << "Temporary save data doesn't exist; skipping\n";
    }
  }
  catch (...) {
    cout << "Temporary save data doesn't exist; skipping\n";
  }

  return;
}



void TJJammer::PrintTime(TimerClass timer)
{
  cout << "Printing Timer results..." << endl;
  //Get the file:
  ofstream outFile;
  string timerFile = "TIMEFILE.txt";

  outFile.open(timerFile.c_str(),ios::out | ios::app);

  if (outFile.is_open()) {
    outFile << timer.timePassed << endl;
    outFile.close();
    cout << "Success!" << endl;
  }
  else {
    cerr << "ERROR: could not write to " << timerFile << endl;
  }

  return;
}


void TJJammer::TJIterations()
{

  // variables
  int printSuccess = 0;
  int ifResize = 0;
  double latVol = 0.0;
  int totNeighbors = 0;    // this is the return if NNLs are used of calcNNLs
  double lastVol = 0.0;
  double lastLastVol = 0.0;

  //If the "rescale" function is asked for, do it here:
  if (rescaleIC > 1.0) RescalePacking();

  //Find the total volume of the spheres:
  double sphereVol = getSphereVol();

  // rescale all pertinent variables, e.g., InternalDelta, etc., by biggestRad. Don't rescale tolerances!!!
  // also note that resizeSpace is not rescaled because resizeSpace is in units of the contact distance b/t two overlapping spheres
  InternalDelta = delta*biggestRad;
  InternalNNLextraDist = NNLextraDist*biggestRad;
  InternalTransMax = transMax*biggestRad;
  InternalCompMax = compMax*biggestRad;
  InternalShearMax = shearMax*biggestRad;
  if(Verbosity>5)
  {
	  cout << "biggestRad = " << biggestRad << "\n";
	  cout << "InternalDelta = " << delta << "\n";
	  cout << "InternalNNLextraDist = " << InternalNNLextraDist << "\n";
	  cout << "InternalTransMax = " << InternalTransMax << "\n";
	  cout << "InternalCompMax = " << InternalCompMax << "\n";
	  cout << "InternalShearMax = " << InternalShearMax << "\n";
  }

  //Initializing the NNL, if required
  if (useNNL) {
    // create the array that holds the neighborNode first addresses
    neighborListDelta = new NeighborNode *[N];
    neighborListOverlap = new NeighborNode *[N];
    for (int i=0; i<N; i++) {
      neighborListDelta[i] = 0;
      neighborListOverlap[i] = 0;
    }
    totNeighbors = calcNNLs();
  if(Verbosity>5)
     cout << "initial NNL maxSize = " << totNeighbors << "\n";
  }

  //Resizing the initial configuration (if needed)
  ifResize = resizeIfOverlap(false);
  if (ifResize) {
    latVol = getVol();
  if(Verbosity>5)
     cout << "Initial resize was required. New lattice volume is " << latVol << "\n";
  }

  //Initial Monte Carlo Sequence (equilibrating the IC)
  MC_Init();
  if (MCStartSweeps > 1.0) 
  {
    for (int mcIndex = 0 ; mcIndex < (int) MCStartSweeps ; mcIndex++) {
  if(Verbosity>7)
       cout << mcIndex+1 << " / " << MCStartSweeps << "\n";
      double mc_maxDisp = MC_Sweep(N);
      //Re-calculate NNL if needed:
      maxSingleMove += mc_maxDisp;
      if (useNNL && (maxSingleMove > InternalNNLextraDist/2.0)) {
	maxSingleMove = 0.0;
  if(Verbosity>7)
 	cout << "Recalculating NNLs\n";
	totNeighbors = calcNNLs();
      }
    }
  }
  else if (MCStartSweeps > 0.0) {//Partial sweep
    int numTries = (int) (MCStartSweeps * N);
    MC_Sweep(numTries);
  }

  //If we wanted to maximize the initial density, rescale the radii here so that the first contact is made:
  if (maxInitRad) {
    MaximizeRadii();
  }


  ////////////////////////////////////////////////////////////
  //
  //  The SLP iterations
  //
  ////////////////////////////////////////////////////////////


  LPClass LP(*this);
  timer.Start();
  saveTimer.Start();
  while (lpIters < maxIters) 
  {
	  if(Verbosity>7)
		  cout << "------------------------------------------" << endl;
	  if(Verbosity>7)
		  cout << "LP Iteration " << lpIters+1 << endl;
	  if (lpIters > 0) {
		  lastLastVol = lastVol;


		  if (useNNL && (maxSingleMove > InternalNNLextraDist/2.0)) {
			  maxSingleMove = 0.0;
			  if(Verbosity>7)
				  cout << "Recalculating NNLs\n";
			  totNeighbors = calcNNLs();
		  }
	  }
	  lastVol = latVol;

	  //Set flags:
	  printThisIter_pkg = (printEvery   > 0) && ((lpIters+1) % printEvery   == 0);
	  printThisIter_LP  = (printLPEvery > 0) && ((lpIters+1) % printLPEvery == 0);

	  //Run the LP routines:
	  LP.LPRoutine(lpIters,totNeighbors);

	  //Do some MC sweeps between iterations, if asked
	  if (MCMidSweeps > 0.0) {
		  int numTries = (int) (MCMidSweeps * N);
		  double mc_maxDisp = MC_Sweep(numTries);
		  //Re-calculate NNL if needed:
		  maxSingleMove += mc_maxDisp;
		  if (useNNL && (maxSingleMove > InternalNNLextraDist/2.0)) {
			  maxSingleMove = 0.0;
			  cout << "Recalculating NNLs\n";
			  totNeighbors = calcNNLs();
		  }
	  }

	  latVol = getVol();
	  double phi = sphereVol / latVol;
	  if(Verbosity>6)
		  cout << "  phi                = " << phi << endl;
	  bool checkTerm = CheckTermination(lastLastVol,latVol,phi);
	  if (checkTerm) {
      break;
    }
    lpIters++;

    if (printThisIter_pkg) {
      printSuccess = printPacking(parentDirectory + outputFilename , true);
      if (printSuccess == 0) {
  if(Verbosity>7)
 	cerr << "Packing did not print for some reason at iters = " << lpIters << "\n";
      }
      else {
  if(Verbosity>7)
 	cout << "Printed packing at iters = " << lpIters << "\n";
      }
    }
    //Check if we should print some save data:
    saveTimer.End();
    if (saveTimer.timePassed > 60*5) {//Save every 5 minutes
      saveTimer.Start(); //Restart the save timer
      if (autoSave) {
	SaveProgress();
      }
    }
  }//LP loop

  if (useTimer == 1) {//Take note of when the simulation started
    timer.End();
    PrintTime(timer);
  }

  if (lpIters == maxIters) {
	  if(Verbosity>5)
		  cout << "Terminated due to maximum iterations exceeded\n";
  }
  if(Verbosity>5)
  {
	  cout << "iterations completed = " << lpIters << "\n";
	  cout << "lastLastVol          = " << lastLastVol << "\n";
	  cout << "lastVol              = " << lastVol << "\n";
	  cout << "latVol               = " << latVol << "\n";
	  cout << "phi               = " << sphereVol/latVol << "\n";
  }
  if (useTimer == 1) timer.Display();

  // print the last packing
  printSuccess = printPacking(parentDirectory + outputFilename , false);

  //If the last print was successful, then delete any temporary save info:
  if (printSuccess && autoSave) {
    ClearSaveData();
  }

  if(sphereVol/latVol>0.01 && sphereVol/latVol<1.01)
  {
	GeometryVector basis[ ::MaxDimension];
	for (int i=0;i<dim;i++) 
	{
		basis[i]=GeometryVector(dim);
		for (int j=0; j<dim; j++) 
			basis[i].x[j]= gsl_matrix_get(lambdas,j,i);
	}

	//ratio is the volume divided by the product of the length of all basis vectors
	double ratio= ::Volume(&basis[0], dim);
	for(int i=0; i<dim; i++)
		ratio/=std::sqrt(basis[i].Modulus2());

	if(ratio>0.01)
	{
		result=SpherePacking(dim, basis, radii[0]*2);
		for(int i=0; i<N; i++)
		{
			GeometryVector rel(dim);
			for(int j=0; j<dim; j++)
				rel.x[j]=localCoords[dim*i+j];
			result.Insert(radii[i], rel);
		}
	}
	else
	{
		//extremely small ratio, simulation box might be weired
		result=SpherePacking();
	}
  }
  else
  {
	  result=SpherePacking();
  }
  // free memory
  gsl_rng_free(RNG);
  gsl_matrix_free(lambdas);
  gsl_matrix_free(inverseLambdas);
  delete [] localCoords;
  delete [] globalCoords;
  delete [] radii;
  delete [] distTempL;
  delete [] distTempG;

  if (useNNL) {
    DeleteOldNNL();
    delete [] neighborListDelta;
    delete [] neighborListOverlap;
  }

  return;
}


bool IsRattler_LP(DimensionType dim, const std::vector<GeometryVector> & contacts)
//Use an LP to find out if sphere i is a rattler
//See Steve's 12/5/2012 notes for explanation
{
	//The return variable:
	bool thisIsRattler = false; //Until proven otherwise
	//Allocate:
	gsl_matrix** rij = new gsl_matrix*[contacts.size()]; //Maximum possible number of contacts
	for (int j = 0; j < contacts.size(); j++) {
		rij[j] = gsl_matrix_alloc(dim, 1);
	}

	//Gather the contacts (made by non-rattlers):
	int numContacts = 0;
	for (int j = 1; j <= contacts.size(); j++) 
	{//Look through its contacts
		//Log it into the vector list:
		//gsl_matrix_memcpy(rij[numContacts], posG[contactList[i][j]]); //get r_ij
		//gsl_matrix_sub(rij[numContacts], posG[i]);
		GeometryVector vec = contacts[j - 1];
		for (int k = 0; k < dim; k++)
			(*gsl_matrix_ptr(rij[numContacts], k, 0)) = vec.x[k];
		numContacts++;
	}//Done testing all of the contacts on the sphere

	//Test for dim+1 contacts at least:
	if (numContacts <= dim) {
		thisIsRattler = true;
	}
	else {//Continue the analysis:
		//Build the LP:
		glp_prob *lp;
		glp_smcp solverParams;
		//glp_bfcp otherParams;
		lp = glp_create_prob();
		//glp_set_prob_name(lp,"LP");
		glp_set_obj_dir(lp, GLP_MIN); //Minimization problem

		// specify control parameters: IMPORTANT!!!
		glp_init_smcp(&solverParams);
		//solverParams.presolve = GLPK_PresolveVar;
		//solverParams.tol_bnd = feasibleTol;
		//solverParams.tol_dj = feasibleTol;
		solverParams.msg_lev = GLP_MSG_ERR; //Only print error and warning messages to the terminal

		glp_add_cols(lp, dim); //Add variables
		for (int j = 0; j < dim; j++) {//Note: counting in the lp object starts at 1, not 0!
			glp_set_obj_coef(lp, j + 1, gsl_matrix_get(rij[0], j, 0));
			glp_set_col_bnds(lp, j + 1, GLP_DB, -1.0, 1.0);
		}
		glp_add_rows(lp, numContacts); //Add constraints
		for (int j = 1; j <= numContacts; j++) {
			glp_set_row_bnds(lp, j, GLP_UP, 0.0, 0.0); //First "0.0" is "lower bound" and gets ignored
		}
		//Make the constraints into the nice kind of matrix that glpk likes:
		int    *ia = new int[numContacts*dim + 1]; //GLPK and its silly start-at-one thing...
		int    *ja = new int[numContacts*dim + 1];
		double *ar = new double[numContacts*dim + 1];
		int index = 0;
		for (int j = 0; j < numContacts; j++) {
			for (int k = 0; k < dim; k++) {
				index++;
				ia[index] = j + 1; //Stupid +1's...
				ja[index] = k + 1;
				ar[index] = gsl_matrix_get(rij[j], k, 0);
			}
		}
		//Solve the LP:
		glp_load_matrix(lp, index, ia, ja, ar);
		int solverStatus = glp_simplex(lp, &solverParams);

		if (solverStatus == 0) {//No problems
			if (glp_get_obj_val(lp) < 0) {
				//cout << "Found rattler (obj val = " << glp_get_obj_val(lp) << ")" << endl;
				thisIsRattler = true;
			}
		}
		else {
			cout << "ERROR running GLPK!" << endl;
		}

		//Dealloacte the LP stuff:
		//cout << "Delete LP stuff..." << flush;
		glp_delete_prob(lp);
		//glp_free_env();
		delete[] ia;
		delete[] ja;
		delete[] ar;
		//cout << "done" << endl;
	}

	//Deallocate:
	for (int j = 0; j < contacts.size(); j++) {
		gsl_matrix_free(rij[j]);
	}
	delete[] rij;

	return(thisIsRattler);
}

bool IsRattler_Count(DimensionType dim, const std::vector<GeometryVector> & contacts)
{
	return contacts.size() <= dim;
}

#endif
