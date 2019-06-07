 
/*
   Example 5

   Interface:    Linear-Algebraic (IJ)
 
   Compile with: make ex5

   Sample run:   mpirun -np 4 ex5

   Description:  This example solves the 2-D Laplacian problem with zero boundary
                 conditions on an n xxx n grid.  The number of unknowns is N=n^2.
                 The standard 5-point stencil is used, and we solve for the
                 interior nodes only.

                 This example solves the same problem as Example 3.  Available
                 solvers are AMG, PCG, and PCG with AMG or Parasails
                 preconditioners.  */





  

#include "Element.h"
#include "solver.h"

using namespace std;

int hypre_FlexGMRESModifyPCAMGExample(HYPRE_Solver precond_data, int iterations, double rel_residual_norm)
{
/*--------------------------------------------------------------------------
   hypre_FlexGMRESModifyPCAMGExample -

    This is an example (not recommended)
   of how we can modify things about AMG that
   affect the solve phase based on how FlexGMRES is doing...For
   another preconditioner it may make sense to modify the tolerance..

 *--------------------------------------------------------------------------*/
   if (rel_residual_norm > .1)
      HYPRE_BoomerAMGSetNumSweeps(precond_data, 10);
   else
      HYPRE_BoomerAMGSetNumSweeps(precond_data, 1);
   return 0;
}















const double            hypreSolver::SCALE_FACTOR = 1.0e12;
bool                    hypreSolver::MPIINITIALISED = false;
bool                    hypreSolver::INITIALISED = false;
bool                    hypreSolver::USE_GRAVITY = false;
double                  hypreSolver::TOLERANCE = 1.0E-29;















/**
// This is the constructor for the solver class. Given the network it allocates memory
// for various C type arrays. These are set to the maximum length possible which is when
// all pores contain the fluid for which the flowrate is going to be obtained.
*/
hypreSolver::hypreSolver(const vector<Element*>& network, const vector<Element*>& inlet, const vector<Element*>& outlet,
               int outletIdx, int debugMode, string matFileName, bool matlabFormat) 
 : m_debugMode(debugMode), m_netelems(network), m_inPores(inlet), m_outPores(outlet), m_numPoresp1(outletIdx), 
   m_matrixFileName(matFileName),    m_poreiRows(outletIdx+1,-1) // m_probSize(outletIdx-1),

{

	m_rowedPores.reserve(network.size());
}






/**
 * This function computes the outlet flow of a
 * given fluid from the network given the inlet and outlet pressures. The solved pressure
 * field is not retained. The matrix is solved using a sparse matrix solver (BiCG) for
 * quick solution times. No assumptions are being made about the structure of the matrix.
*/
double hypreSolver::flowrate(double inletPrs, double outletPrs, const Fluid &fluid, double& flowError, double& elapsed , double a, bool b, bool c, bool resistivitySolve)
{
    clock_t  startSolve(clock());


	if (fluid.isOil())
		for (size_t i = 0; i < m_netelems.size(); ++i)
			m_netelems[i]->clearAllOSolverFlag();
	else
		for (size_t i = 0; i < m_netelems.size(); ++i)
			m_netelems[i]->clearAllWSolverFlag();

	bool outletConnectionExists(false);                       // For each pore at inlet boundary a path to the outlet is tried found.
	for (size_t bdrP = 0; bdrP < m_inPores.size(); ++bdrP)   // The pores existing within this path is marked in the process
	{
		if (m_inPores[bdrP]->connectedToOutlet(&fluid))
		{ 
			outletConnectionExists = true;
		}
	} 
 
	if (!outletConnectionExists) return 0.0;                 // No fluid connection exists between inlet and outlet















   int i;
   int myid(0);//, num_procs(1); 
   bool printError(false);



	int solver_id(0);
	int vis(0), print_system(0);


	HYPRE_IJMatrix AAA;
	HYPRE_ParCSRMatrix parcsr_A;
	HYPRE_IJVector bbb;
	HYPRE_ParVector par_b;
	HYPRE_IJVector xxx;
	HYPRE_ParVector par_x;

	HYPRE_Solver solver, precond;

	//int argc(0); char **argv;

	/* Initialize MPI */
	if(!MPIINITIALISED)   {  MPIINITIALISED=true;	}




   /* Preliminaries: want at least one processor per row */
   // if (n*n < num_procs) n = sqrt(num_procs) + 1;
   // N = n*n; /* global number of rows */

   int ilower(0), iupper;
   int local_size(0);//, extra;
   //- /* Each processor knows only of its own rows - the range is denoted by ilower   and upper.  Here we partition the rows. We account for the fact that  N may not divide evenly by the number of processors. */
   //- local_size = N/num_procs;
   //- extra = N - local_size*num_procs;
   //- ilower = local_size*myid;
   //- ilower += hypre_min(myid, extra);
   //- iupper = local_size*(myid+1);
   //- iupper += hypre_min(myid+1, extra);
   //- iupper = iupper - 1;
   //- /* How many rows do I have? */
   //- local_size = iupper - ilower + 1;
   vector<int> rowSizes(m_numPoresp1,0);
   	m_rowedPores.resize(0);
	for(int poreIdx = 1; poreIdx < m_numPoresp1; ++poreIdx)
    {
        if(m_netelems[poreIdx]->canBePassedToSolver(&fluid) && m_netelems[poreIdx]->isInsideSolverBox())
        {
				rowSizes[local_size]=m_netelems[poreIdx]->connectionNum()+1;
				m_poreiRows[poreIdx]=local_size;
				m_rowedPores.push_back(poreIdx);

				++local_size;
        }
    }
	iupper=local_size-1;
	



		/* Create the matrix.
			Note that this is a square matrix, so we indicate the row partition
			size twice (since number of rows = number of cols) */
		HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &AAA);

		HYPRE_IJMatrixSetRowSizes(AAA, &(rowSizes[0]));

		/* Choose a parallel csr format storage (see the User's Manual) */
		HYPRE_IJMatrixSetObjectType(AAA, HYPRE_PARCSR);

		/* Initialize before setting coefficients */
		HYPRE_IJMatrixInitialize(AAA);



		HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&xxx);
		HYPRE_IJVectorSetObjectType(xxx, HYPRE_PARCSR);
		HYPRE_IJVectorInitialize(xxx);

		/* Create the rhs and solution */
		HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&bbb);
		HYPRE_IJVectorSetObjectType(bbb, HYPRE_PARCSR);
		HYPRE_IJVectorInitialize(bbb);




		fillMatrixHypre( AAA, bbb, xxx, inletPrs, outletPrs, &fluid, false, resistivitySolve);

	/* Assemble after setting the coefficients */
	HYPRE_IJMatrixAssemble(AAA);



	/* Get the parcsr matrix object to use */
	HYPRE_IJMatrixGetObject(AAA, (void**) &parcsr_A);


	HYPRE_IJVectorAssemble(bbb);
	HYPRE_IJVectorGetObject(bbb, (void **) &par_b);

	HYPRE_IJVectorAssemble(xxx);
	HYPRE_IJVectorGetObject(xxx, (void **) &par_x);



	/*  Print out the system  - files names will be IJ.out.AAA.XXXXX and IJ.out.bbb.XXXXX, where XXXXX = processor id */
	if (print_system)
	{
		HYPRE_IJMatrixPrint(AAA, "IJ.out.AAA");
		HYPRE_IJVectorPrint(bbb, "IJ.out.bbb");
	}




	///. Choose a solver and solve the system
  {/* AMG */
	if (solver_id == 0)
	{
		int num_iterations;
		double final_res_norm;

		/* Create solver */
		HYPRE_BoomerAMGCreate(&solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_BoomerAMGSetPrintLevel(solver, 0);  /* print solve info + parameters */
		HYPRE_BoomerAMGSetCoarsenType(solver, 6); /* Falgout coarsening */
		HYPRE_BoomerAMGSetRelaxType(solver, 3);   /* G-S/Jacobi hybrid relaxation */
		HYPRE_BoomerAMGSetNumSweeps(solver, 1);   /* Sweeeps on each level */
		HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
		HYPRE_BoomerAMGSetTol(solver, 1e-12);      /* conv. tolerance */

		/* Now setup and solve! */
		HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
		HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

		/* Run info - needed logging turned on */
		HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
		HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		if (myid == 0 && printError)
		{
			printf("\n");
			printf("Iterations = %d\n", num_iterations);
			printf("Final Relative Residual Norm = %e\n", final_res_norm);
			printf("\n");
		}

		/* Destroy solver */
		HYPRE_BoomerAMGDestroy(solver);
	}
	/* PCG */
	else if (solver_id == 50)
	{
		int num_iterations;
		double final_res_norm;

		/* Create solver */
		HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
		HYPRE_PCGSetTol(solver, 1e-9); /* conv. tolerance */
		HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
		HYPRE_PCGSetPrintLevel(solver, 0); /* prints out the iteration info */
		HYPRE_PCGSetLogging(solver, 0); /* needed to get run info later */

		/* Now setup and solve! */
		HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
		HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

		/* Run info - needed logging turned on */
		HYPRE_PCGGetNumIterations(solver, &num_iterations);
		HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		if (myid == 0 && printError)
		{
			printf("\n");
			printf("Iterations = %d\n", num_iterations);
			printf("Final Relative Residual Norm = %e\n", final_res_norm);
			printf("\n");
		}

		/* Destroy solver */
		HYPRE_ParCSRPCGDestroy(solver);
	}
	/* PCG with AMG preconditioner */
	else if (solver_id == 1)
	{
		int num_iterations;
		double final_res_norm;

		/* Create solver */
		HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
		HYPRE_PCGSetTol(solver, 1e-9); /* conv. tolerance */
		HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
		HYPRE_PCGSetPrintLevel(solver, 0); /* print solve info */
		HYPRE_PCGSetLogging(solver, 0); /* needed to get run info later */

		/* Now set up the AMG preconditioner and specify any parameters */
		HYPRE_BoomerAMGCreate(&precond);
		HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print amg solution info */
		HYPRE_BoomerAMGSetCoarsenType(precond, 6);
		HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
		HYPRE_BoomerAMGSetNumSweeps(precond, 1);
		HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
		HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

		/* Set the PCG preconditioner */
		HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							 (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);

		/* Now setup and solve! */
		HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
		HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

		/* Run info - needed logging turned on */
		HYPRE_PCGGetNumIterations(solver, &num_iterations);
		HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		if (myid == 0 && printError)
		{
			printf("\n");
			printf("Iterations = %d\n", num_iterations);
			printf("Final Relative Residual Norm = %e\n", final_res_norm);
			printf("\n");
		}

		/* Destroy solver and preconditioner */
		HYPRE_ParCSRPCGDestroy(solver);
		HYPRE_BoomerAMGDestroy(precond);
	}
	/* PCG with Parasails Preconditioner */
	else if (solver_id == 8)
	{
		int    num_iterations;
		double final_res_norm;

		int      sai_max_levels = 1;
		double   sai_threshold = 0.1;
		double   sai_filter = 0.05;
		int      sai_sym = 1;

		/* Create solver */
		HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
		HYPRE_PCGSetTol(solver, 1e-9); /* conv. tolerance */
		HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
		HYPRE_PCGSetPrintLevel(solver, 0); /* print solve info */
		HYPRE_PCGSetLogging(solver, 0); /* needed to get run info later */

		/* Now set up the ParaSails preconditioner and specify any parameters */
		HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &precond);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_ParaSailsSetParams(precond, sai_threshold, sai_max_levels);
		HYPRE_ParaSailsSetFilter(precond, sai_filter);
		HYPRE_ParaSailsSetSym(precond, sai_sym);
		HYPRE_ParaSailsSetLogging(precond, 0);

		/* Set the PCG preconditioner */
		HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
					(HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);

		/* Now setup and solve! */
		HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
		HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);


		/* Run info - needed logging turned on */
		HYPRE_PCGGetNumIterations(solver, &num_iterations);
		HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		if (myid == 0 && printError)
		{
			printf("\n");
			printf("Iterations = %d\n", num_iterations);
			printf("Final Relative Residual Norm = %e\n", final_res_norm);
			printf("\n");
		}

		/* Destory solver and preconditioner */
		HYPRE_ParCSRPCGDestroy(solver);
		HYPRE_ParaSailsDestroy(precond);
	}
	/* Flexible GMRES with  AMG Preconditioner */
	else if (solver_id == 61)
	{
		int    num_iterations;
		double final_res_norm;
		int    restart = 30;
		int    modify = 1;


		/* Create solver */
		HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_FlexGMRESSetKDim(solver, restart);
		HYPRE_FlexGMRESSetMaxIter(solver, 1000); /* max iterations */
		HYPRE_FlexGMRESSetTol(solver, 1e-9); /* conv. tolerance */
		HYPRE_FlexGMRESSetPrintLevel(solver, 0); /* print solve info */
		HYPRE_FlexGMRESSetLogging(solver, 0); /* needed to get run info later */


		/* Now set up the AMG preconditioner and specify any parameters */
		HYPRE_BoomerAMGCreate(&precond);
		HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print amg solution info */
		HYPRE_BoomerAMGSetCoarsenType(precond, 6);
		HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
		HYPRE_BoomerAMGSetNumSweeps(precond, 1);
		HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
		HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

		/* Set the FlexGMRES preconditioner */
		HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							 (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);


		if (modify)
		/* this is an optional call  - if you don't call it, hypre_FlexGMRESModifyPCDefault
		is used - which does nothing.  Otherwise, you can define your own, similar to
		the one used here */
		HYPRE_FlexGMRESSetModifyPC( solver, (HYPRE_PtrToModifyPCFcn) hypre_FlexGMRESModifyPCAMGExample);


		/* Now setup and solve! */
		HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);
		HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);

		/* Run info - needed logging turned on */
		HYPRE_FlexGMRESGetNumIterations(solver, &num_iterations);
		HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
		if (myid == 0 && printError)
		{
			printf("\n");
			printf("Iterations = %d\n", num_iterations);
			printf("Final Relative Residual Norm = %e\n", final_res_norm);
			printf("\n");
		}

		/* Destory solver and preconditioner */
		HYPRE_ParCSRFlexGMRESDestroy(solver);
		HYPRE_BoomerAMGDestroy(precond);

	}
	else
	{
		if (myid ==0) printf("Invalid solver id specified.\n");
	}
  }

	/* Save the solution for GLVis visualization, see vis/glvis-ex5.sh */
	if (vis)
	{
		FILE *file;
		char filename[255];

		int nvalues = local_size;
		int *rows = new  int[nvalues];
		double *values = new  double[nvalues];

		for (i = 0; i < nvalues; i++)       rows[i] = ilower + i;

		/* get the local solution */
		HYPRE_IJVectorGetValues(xxx, nvalues, rows, values);

		sprintf(filename, "%s.%06d", "vis/ex5.sol", myid);
		if ((file = fopen(filename, "w")) == NULL)
		{
			printf("Error: can't open output file %s\n", filename);
			MPI_Finalize();
			exit(1);
		}

		/* save solution */
		for (i = 0; i < nvalues; i++)     fprintf(file, "%.14e\n", values[i]);

		fflush(file);
		fclose(file);

		free(rows);
		free(values);

		/* save global finite element mesh */
		//- if (myid == 0)  GLVis_PrintGlobalSquareMesh("vis/ex5.mesh", n-1);
	}

	int m_nMatrixRows=m_rowedPores.size();
	vector<double> porePressure(m_nMatrixRows,0.0);
	vector<int>  poreRows(m_nMatrixRows,0);
	for(int i = 0; i < m_nMatrixRows; ++i)    poreRows[i]=i;
	HYPRE_IJVectorGetValues(xxx, m_nMatrixRows, &poreRows[0], &porePressure[0]);

 
	for(int j = 0; j < m_nMatrixRows; ++j)                                           // Pass back the results from the
		((Pore*)(m_netelems[m_rowedPores[j]]))->setSolverResults(&fluid, resistivitySolve, porePressure[j]);   // solver. These values will// subsequently be used for calculating



    for(size_t inp = 0; inp < m_inPores.size(); ++inp)                              // pressure profiles on the lattice
    {
        Pore *inletPore = dynamic_cast< Pore* >(m_inPores[inp]);
        if(inletPore)
        {
            double localInletPrs(inletPrs);
            if(!resistivitySolve && USE_GRAVITY)
            {
                localInletPrs += inletPore->rhogh(fluid.density(), inletPore->node()->xPos(),
                    inletPore->node()->yPos(), inletPore->node()->zPos());
            }

            inletPore->setSolverResults(&fluid, resistivitySolve, localInletPrs);
        }
    }
 
 
    for(size_t op = 0; op < m_outPores.size(); ++op)
    {
        Pore *outletPore = dynamic_cast< Pore* >(m_outPores[op]);
        if(outletPore)
        {
            double localOutletPrs(outletPrs);
            if(!resistivitySolve && USE_GRAVITY)
            {
                localOutletPrs += outletPore->rhogh(fluid.density(), outletPore->node()->xPos(),
                    outletPore->node()->yPos(), outletPore->node()->zPos());
            }

            outletPore->setSolverResults(&fluid, resistivitySolve, localOutletPrs);
        }
	}








	int nErrors = 0;

	double flowRate = getFlowRate(xxx, &fluid, flowError, nErrors, resistivitySolve);

	elapsed += (double)(clock() - startSolve) / CLOCKS_PER_SEC;

	if (nErrors>0) cout<< "\n" << nErrors << " errors when solving for " << (fluid.isOil() ? "Oil" : "Water") << " flowrate           "; cout.flush();




	//. Clean up
	HYPRE_IJMatrixDestroy(AAA);
	HYPRE_IJVectorDestroy(bbb);
	HYPRE_IJVectorDestroy(xxx);

	return flowRate;

}
 


/**
* Only the pressures in the pores connected to the outlet are actually needed to compute
* the flowrate leaving the network. The conductances and index pointers to the solution
* vector are stored in a vector that is filled during matrix filling.
*/
double hypreSolver::getFlowRate(HYPRE_IJVector xxx, const Fluid* fluid, double& maxError, int& nErrors, bool resistSolve) const
{
    double flowOut(0.0), flowIn(0.0);

	vector<double> outletPrs(m_networkOutlets.size(),0.0);
	vector<int>  outletRows(m_networkOutlets.size(),0);
    for(size_t i = 0; i < m_networkOutlets.size(); ++i) 
		outletRows[i]=m_networkOutlets[i].first();
	HYPRE_IJVectorGetValues(xxx, m_networkOutlets.size(), &outletRows[0], &outletPrs[0]);
    for(size_t i = 0; i < m_networkOutlets.size(); ++i)
    {
        double porePrs = outletPrs[i];
        double outletPrs = m_networkOutlets[i].third();
        double conductance = m_networkOutlets[i].second() / SCALE_FACTOR;
        if (porePrs < -1.4e24 && porePrs > -1.6e24)
        {
				continue;
		}
		if (porePrs<outletPrs-2)
		{
			++nErrors;
		}
        flowOut += conductance * (porePrs - outletPrs /*- rhogh*/);
    }




	vector<double> inletPrs(m_networkInlets.size(),0.0);
	vector<int>  inletRows(m_networkInlets.size(),0);
	for(size_t i = 0; i < m_networkInlets.size(); ++i)   	inletRows[i]=m_networkInlets[i].first();
	HYPRE_IJVectorGetValues(xxx, m_networkInlets.size(), &inletRows[0], &inletPrs[0]);
    for(size_t j = 0; j < m_networkInlets.size(); ++j)
    {
        double porePrs = inletPrs[j];
        double inletPrs = m_networkInlets[j].third();
        double conductance = m_networkInlets[j].second() / SCALE_FACTOR;
        if (porePrs < -1.4e24 && porePrs > -1.6e24)
        {
				continue;
		}
		if (porePrs>inletPrs+2 && m_debugMode)
		{
			++nErrors;
		}
        flowIn += conductance * (inletPrs - porePrs /*+ rhogh*/);
    }

    double flowError = fabs(flowOut - flowIn) / flowOut;

       if ((flowError>1.0 || flowError<-1.0)  && m_debugMode)
        {
				cout<<"   q_i-q_o = "<<flowIn<<" "<<flowOut<<"   ";cout.flush();
		}
	maxError=max(maxError,flowError);

	return (flowOut + flowIn) / 2.0;
}










/**
// Both the matrix and rhs vectors are filled. The matrix is filled taking advantage of it's
// sparse property. Hence it is represented in terms a three arrays; value, row pointer and
// column index. The matrix size and number of non-zero elements are also determined. The rhs
// vector is stored in a C type array rather than the MV format since the required size is
// not known during filling. It'll have to be converted before solving the system.
*/
void hypreSolver::fillMatrixHypre(HYPRE_IJMatrix AAA, HYPRE_IJVector bbb, HYPRE_IJVector xxx, double inletPrs, double outletPrs, const Fluid *fluid, bool writeVelocity, bool resistivitySolve)
{
	int nnz;
	double values[1000];
	int cols[1000];

	int row(0), poreIdx, conn;
	m_networkOutlets.clear();                               // Ensure we don't have left over soliutions from before
	m_networkInlets.clear();
	for(poreIdx = 1; poreIdx < m_numPoresp1; ++poreIdx)
	{
		Pore *currPore = dynamic_cast<Pore*>(m_netelems[poreIdx]);
		if(currPore == NULL)
		{
			cerr << "\nError: Did not retrieve correct pore in solver **** \n" << endl;
			exit(-1);
		}


		if(currPore->canBePassedToSolver(fluid) && currPore->isInsideSolverBox())
		{

			register double conductanceSum(0.0);
			register double rhs(0.0);
			FourSome<int, double, double, double> inletPoint(-1, 0, 0, 0);
			FourSome<int, double, double, double> outletPoint(-1, 0, 0, 0);




			nnz=1;
			values[0]=1.0e-64;
			cols[0]=row;
			for(conn = 0; conn < currPore->connectionNum(); ++conn)
			{
				double conductance(0.0);
				double deltaGrav(0.0);

				const Element* throat = currPore->connection(conn);
				Element *nextPore = currPore->getConnectionProp(conn, conductance, deltaGrav, fluid, resistivitySolve);

				USE_GRAVITY = false;
				ensure(conductance >= 0.0);
				if(!USE_GRAVITY) deltaGrav = 0.0;
				conductance *= SCALE_FACTOR;            // Some of the preconditioners will assume a value to be zero if very small. Let's just scale the matrix a bit...


				if(conductance > 1.0e-165)
				{
					int nextPoreIndex(nextPore->node()->index());

					if(writeVelocity)
					{
						pair<const Element*, double> veloEntry(throat, conductance / SCALE_FACTOR);
						m_throatConductances.push_back(veloEntry);
					}
					conductanceSum += conductance;
					values[0]+=conductance;
					if(USE_GRAVITY) rhs += conductance*deltaGrav;
 
					if(nextPore->isOnInletSlvrBdr())                            // Inlet///////////////////////////////////////////////
					{
						double localInletPrs(inletPrs);
						if(!resistivitySolve && USE_GRAVITY && nextPore->isEntryOrExitRes())
							localInletPrs += currPore->rhogh(fluid->density(), currPore->node()->xPos(), currPore->node()->yPos(), currPore->node()->zPos());
						else if(!resistivitySolve && USE_GRAVITY && !nextPore->isEntryOrExitRes()) 
							localInletPrs += nextPore->rhogh(fluid->density(), nextPore->node()->xPos(), nextPore->node()->yPos(), nextPore->node()->zPos());


						inletPoint.first(row); inletPoint.second(inletPoint.second()+conductance);  inletPoint.third(localInletPrs);  inletPoint.fourth(deltaGrav); 
						rhs += conductance * localInletPrs;
					}
					else if(nextPore->isOnOutletSlvrBdr())                       // Outlet////////////////////////////////////////////////
					{ 
						double localOutletPrs(outletPrs);
						if(!resistivitySolve && USE_GRAVITY && nextPore->isEntryOrExitRes())
							localOutletPrs += currPore->rhogh(fluid->density(), currPore->node()->xPos(), currPore->node()->yPos(), currPore->node()->zPos());
						else if(!resistivitySolve && USE_GRAVITY && !nextPore->isEntryOrExitRes())
							localOutletPrs += nextPore->rhogh(fluid->density(), nextPore->node()->xPos(), nextPore->node()->yPos(), nextPore->node()->zPos());


						outletPoint.first(row); outletPoint.second(outletPoint.second()+conductance);  outletPoint.third(localOutletPrs);  outletPoint.fourth(deltaGrav); 

						rhs += conductance * localOutletPrs;
					}
					else                                                        // Interior point
					{
						values[nnz]=-conductance;
						cols[nnz++]=m_poreiRows[nextPoreIndex];
					}
				}
			}



			if (conductanceSum>1.0e-165)
			{
			
				
				HYPRE_IJMatrixAddToValues(AAA, 1, &nnz, &row, cols, values);

				HYPRE_IJVectorSetValues(bbb, 1, &row, &rhs);
					
				if(inletPoint.first()>= 0)
				{
					m_networkInlets.push_back(inletPoint);
					HYPRE_IJVectorSetValues(xxx, 1, &row, &inletPrs);
				}
				if(outletPoint.first()>= 0)
				{
					m_networkOutlets.push_back(outletPoint);
					HYPRE_IJVectorSetValues(xxx, 1, &row, &outletPrs);
				}
				else
					HYPRE_IJVectorSetValues(xxx, 1, &row, &inletPrs);


				++row;

			}
			else
			{
				 currPore->setSolverResults(fluid, resistivitySolve, -1.5e24);   // solver. These values will
			}
		}
		else
		{
			 currPore->setSolverResults(fluid, resistivitySolve, -1.5e24);   // solver. These values will
		}
	}

}












