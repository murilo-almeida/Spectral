cmake_minimum_required(VERSION 3.9)

set(INSTALL_DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

set(CMAKE_INSTALL_RPATH ${MY_RPATH})

#Bring the headers, such as Stdel.h into the project

include_directories(${MY_INCLUDE})
link_directories(${LIB_DIRS})

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "*.cpp")

#Generate the shared library from the sources
add_executable(dg ${SOURCES})
target_link_libraries(dg Spectral ${LIBS})

get_target_property(MY_R dg BUILD_RPATH)
message("BUILD_RPATH ${MY_R}")

get_target_property(MY_R dg INSTALL_RPATH)
message("INSTALL_RPATH ${MY_R}")

install(TARGETS dg DESTINATION ${INSTALL_DESTINATION})


