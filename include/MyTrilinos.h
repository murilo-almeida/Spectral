#ifndef _MyTrilinos
#define _MyTrilinos
//#include "Epetra_ConfigDefs.h"
#include "Amesos_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_Time.h"
//#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Teuchos_ParameterList.hpp"
//#include "Galeri_Maps.h"
//#include "Galeri_CrsMatrices.h"
#include "AztecOO.h"
//#include "AztecOO_Version.h"

//using namespace Teuchos;
//using namespace Galeri;

//#include "Trilinos_Util.h"
//#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
#include "CrsMatrixTranspose.h"
#include "Teuchos_RCP.hpp"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

#endif
