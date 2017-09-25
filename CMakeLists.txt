cmake_minimum_required(VERSION 3.9)
project(Spectral_test)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_INSTALL_PREFIX "/Users/murilo/MeusProjetos/Spectral")
MESSAGE(" Cmake Instal Prefix = ${CMAKE_INSTALL_PREFIX}") 


SET(Trilinos_PATH "/opt/local/Trilinos")

find_package(mpi)
find_package(boost)

# Your "do-configure" script that invokes CMake should set
# TRILINOS_PATH to the path to your Trilinos install.
# You do _not_ need to edit this line.
FIND_PACKAGE(Trilinos CONFIGS "TrilinosConfig.cmake" HINTS "/opt/local/Trilinos/lib/cmake/Trilinos")

# If FIND_PACKAGE successfully found your Trilinos install, it will
# set the Boolean flag Trilinos_FOUND.  The following IF statement
# fails with a FATAL_ERROR if Trilinos was not found.  If it _was_
# found, it prints out the values of some Trilinos configuration
# details.  You may find them useful for building your application
# that uses Trilinos.
IF(Trilinos_FOUND)
   MESSAGE("\nFound Trilinos!  Here are the details: ")
   MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
   MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
#   MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
#   MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
   MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
#   MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
#   MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
#   MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
   MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
   MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
   MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
#   MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
   MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
#   MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
#   MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS =
#     ${Trilinos_Fortran_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
   MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
   MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
   MESSAGE("End of Trilinos details\n")
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos!")
ENDIF()


SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
MESSAGE("Setting CXX compiler = ${CMAKE_CXX_COMPILER}")
set(CMAKE_CXX_FLAGS  ${Trilinos_CXX_COMPILER_FLAGS})

MESSAGE("Setting C compiler.")
SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})



###################################################################

set(MY_INCLUDE ${Boost_INCLUDE_DIR}
		 ${MPI_CXX_INCLUDE_PATH}
		 "/opt/local/Trilinos/include"
		 "/opt/local/include/newmat"
		 "${CMAKE_INSTALL_PREFIX}/include")

set(LIB_DIRS "/opt/local/lib" "${CMAKE_INSTALL_PREFIX}/lib" "/opt/local/Trilinos/lib")
set(MY_RPATH "${CMAKE_INSTALL_PREFIX}/lib" "/opt/local/Trilinos/lib")
set(LIBS "newmat" ${Trilinos_LIBRARIES})
message("MY_RPATH = ${MY_RPATH}")

add_subdirectory(src/Stdel)
add_subdirectory(src/DG)
