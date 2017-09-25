#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBasicSort.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "spectral.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

// ***********************************************************
int DG_Prob::Eigenvectors(const double Dt,
                          const Epetra_Map & Map)
{
  printf("Entrou em Eigenvectors\n");
 
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  //MPI::COMM_WORLD.Barrier();
  Comm.Barrier();
  Teuchos::RCP<Epetra_FECrsMatrix> M = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, Map,0));//&NNz[0]);
  Teuchos::RCP<Epetra_FEVector> RHS = Teuchos::rcp(new Epetra_FEVector(Map,1));
  
  DG_MatrizVetor_Epetra(Dt,M,RHS);

  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp(new Epetra_CrsMatrix(Copy, Map,0
                                                                       /* &NNz[0]*/) );
  Epetra_Export Exporter(Map,Map);
  A->PutScalar(0.0);
  A->Export(*(M.ptr()),Exporter,Add);
  A->FillComplete();

  using std::cout;
 // int nx = 5;
  bool boolret;
  int MyPID = Comm.MyPID();
  
  bool verbose = true;
  bool debug = false;
  std::string which("LR");
  
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM,LM,SR,LR,SI,or LI).");
  
  typedef double ScalarType;
  typedef Teuchos::ScalarTraits<ScalarType>          SCT;
  typedef SCT::magnitudeType               MagnitudeType;
  typedef Epetra_MultiVector                          MV;
  typedef Epetra_Operator                             OP;
  typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
  typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;
  
 
 // double rho = 2*nx+1;
  
  // Compute coefficients for discrete convection-diffution operator
 // const double one = 1.0;
  
 // int NumEntries, info;
  
  //************************************
  // Start the block Arnoldi iteration
  //***********************************
  //
  //  Variables used for the Block Krylov Schur Method
  //    
  int nev = 10;
  int blockSize = 1;
  int numBlocks = 20;
  int maxRestarts = 500;
  //int stepSize = 5;
  double tol = 1e-8;
  
  // Create a sort manager to pass into the block Krylov-Schur solver manager
  // -->  Make sure the reference-counted pointer is of type Anasazi::SortManager<>
  // -->  The block Krylov-Schur solver manager uses Anasazi::BasicSort<> by default,
  //      so you can also pass in the parameter "Which", instead of a sort manager.
  Teuchos::RCP<Anasazi::SortManager<MagnitudeType> > MySort =     
    Teuchos::rcp( new Anasazi::BasicSort<MagnitudeType>( which ) );
  
  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }
  if (debug) {
    verbosity += Anasazi::Debug;
  }
  //
  // Create parameter list to pass into solver manager
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Sort Manager", MySort );
  //MyPL.set( "Which", which );  
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  //MyPL.set( "Step Size", stepSize );
  MyPL.set( "Convergence Tolerance", tol );
  
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
  ivec->Random();
  
  // Create the eigenproblem.
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec) );
  
  // Inform the eigenproblem that the operator A is symmetric
  //MyProblem->setHermitian(rho==0.0); 
  
  // Set the number of eigenvalues requested
  MyProblem->setNEV( nev );
  
  // Inform the eigenproblem that you are finishing passing it information
  boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }
  
  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);
  
  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && MyPID==0 && verbose) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
  }
  
  // Get the Ritz values from the eigensolver
  std::vector<Anasazi::Value<double> > ritzValues = MySolverMgr.getRitzValues();
  
  // Output computed eigenvalues and their direct residuals
  if (verbose && MyPID==0) {
    int numritz = (int)ritzValues.size();
    cout.setf(std::ios_base::right, std::ios_base::adjustfield);
    cout<<endl<< "Computed Ritz Values"<< endl;
    if (MyProblem->isHermitian()) {
      cout<< std::setw(16) << "Real Part"
	  << endl;
      cout<<"-----------------------------------------------------------"<<endl;
      for (int i=0; i<numritz; i++) {
        cout<< std::setw(16) << ritzValues[i].realpart 
	    << endl;
      }  
      cout<<"-----------------------------------------------------------"<<endl;
    } 
    else {
      cout<< std::setw(16) << "Real Part"
	  << std::setw(16) << "Imag Part"
	  << endl;
      cout<<"-----------------------------------------------------------"<<endl;
      for (int i=0; i<numritz; i++) {
        cout<< std::setw(16) << ritzValues[i].realpart 
	    << std::setw(16) << ritzValues[i].imagpart 
	    << endl;
      }  
      cout<<"-----------------------------------------------------------"<<endl;
    }  
  }
  
  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ScalarType,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<ScalarType> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  std::vector<int> index = sol.index;
  int numev = sol.numVecs;
  
  if (numev > 0) {
    // Compute residuals.
    Teuchos::LAPACK<int,double> lapack;
    std::vector<double> normA(numev);
    
    if (MyProblem->isHermitian()) {
      // Get storage
      Epetra_MultiVector Aevecs(Map,numev);
      Teuchos::SerialDenseMatrix<int,double> B(numev,numev);
      B.putScalar(0.0); 
      for (int i=0; i<numev; i++) {B(i,i) = evals[i].realpart;}
      
      // Compute A*evecs
      OPT::Apply( *A, *evecs, Aevecs );
      
      // Compute A*evecs - lambda*evecs and its norm
      MVT::MvTimesMatAddMv( -1.0, *evecs, B, 1.0, Aevecs );
      MVT::MvNorm( Aevecs, normA );
      
      // Scale the norms by the eigenvalue
      for (int i=0; i<numev; i++) {
        normA[i] /= Teuchos::ScalarTraits<double>::magnitude( evals[i].realpart );
      }
    } else {
      // The problem is non-Hermitian.
      int i=0;
      std::vector<int> curind(1);
      std::vector<double> resnorm(1), tempnrm(1);
      Teuchos::RCP<MV> tempAevec;
      Teuchos::RCP<const MV> evecr, eveci;
      Epetra_MultiVector Aevec(Map,numev);
      
      // Compute A*evecs
      OPT::Apply( *A, *evecs, Aevec );
      
      Teuchos::SerialDenseMatrix<int,double> Breal(1,1), Bimag(1,1);
      while (i<numev) {
        if (index[i]==0) {
          // Get a view of the current eigenvector (evecr)
          curind[0] = i;
          evecr = MVT::CloneView( *evecs, curind );
	  
          // Get a copy of A*evecr
          tempAevec = MVT::CloneCopy( Aevec, curind );
	  
          // Compute A*evecr - lambda*evecr
          Breal(0,0) = evals[i].realpart;
          MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );
	  
          // Compute the norm of the residual and increment counter
          MVT::MvNorm( *tempAevec, resnorm );
          normA[i] = resnorm[0]/Teuchos::ScalarTraits<MagnitudeType>::magnitude( evals[i].realpart );
          i++;
        } else {
          // Get a view of the real part of the eigenvector (evecr)
          curind[0] = i;
          evecr = MVT::CloneView( *evecs, curind );
	  
          // Get a copy of A*evecr
          tempAevec = MVT::CloneCopy( Aevec, curind );
	  
          // Get a view of the imaginary part of the eigenvector (eveci)
          curind[0] = i+1;
          eveci = MVT::CloneView( *evecs, curind );
	  
          // Set the eigenvalue into Breal and Bimag
          Breal(0,0) = evals[i].realpart;
          Bimag(0,0) = evals[i].imagpart;
	  
          // Compute A*evecr - evecr*lambdar + eveci*lambdai
          MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );
          MVT::MvTimesMatAddMv( 1.0, *eveci, Bimag, 1.0, *tempAevec );
          MVT::MvNorm( *tempAevec, tempnrm );
	  
          // Get a copy of A*eveci
          tempAevec = MVT::CloneCopy( Aevec, curind );
	  
          // Compute A*eveci - eveci*lambdar - evecr*lambdai
          MVT::MvTimesMatAddMv( -1.0, *evecr, Bimag, 1.0, *tempAevec );
          MVT::MvTimesMatAddMv( -1.0, *eveci, Breal, 1.0, *tempAevec );
          MVT::MvNorm( *tempAevec, resnorm );
	  
          // Compute the norms and scale by magnitude of eigenvalue
          normA[i] = lapack.LAPY2( tempnrm[i], resnorm[i] ) /
            lapack.LAPY2( evals[i].realpart, evals[i].imagpart );
          normA[i+1] = normA[i];
	  
          i=i+2;
        }
      }
    }
    
    // Output computed eigenvalues and their direct residuals
    if (verbose && MyPID==0) {
      cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      cout<<endl<< "Actual Residuals"<<endl;
      if (MyProblem->isHermitian()) {
        cout<< std::setw(16) << "Real Part"
	    << std::setw(20) << "Direct Residual"<< endl;
        cout<<"-----------------------------------------------------------"<<endl;
        for (int i=0; i<numev; i++) {
          cout<< std::setw(16) << evals[i].realpart 
	      << std::setw(20) << normA[i] << endl;
        }  
        cout<<"-----------------------------------------------------------"<<endl;
      } 
      else {
        cout<< std::setw(16) << "Real Part"
	    << std::setw(16) << "Imag Part"
	    << std::setw(20) << "Direct Residual"<< endl;
        cout<<"-----------------------------------------------------------"<<endl;
        for (int i=0; i<numev; i++) {
          cout<< std::setw(16) << evals[i].realpart 
	      << std::setw(16) << evals[i].imagpart 
	      << std::setw(20) << normA[i] << endl;
        }  
        cout<<"-----------------------------------------------------------"<<endl;
      }  
    }
  }
  
#ifdef EPETRA_MPI
  MPI_Finalize();
#endif  
  return 0;
}
