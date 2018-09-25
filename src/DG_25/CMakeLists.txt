cmake_minimum_required(VERSION 3.9)

set(INSTALL_DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

set(CMAKE_INSTALL_RPATH ${MY_RPATH})

#Bring the headers, such as Stdel.h into the project

include_directories(${MY_INCLUDE})
link_directories(${LIB_DIRS})

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "*.cpp")

set(myexec "dg")
#Generate the shared library from the sources
add_executable(${myexec} ${SOURCES})
target_link_libraries(${myexec} Spectral ${LIBS})

get_target_property(MY_R ${myexec} BUILD_RPATH)
message("BUILD_RPATH ${MY_R}")

get_target_property(MY_R ${myexec} INSTALL_RPATH)
message("INSTALL_RPATH ${MY_R}")

install(TARGETS ${myexec} DESTINATION ${INSTALL_DESTINATION})


