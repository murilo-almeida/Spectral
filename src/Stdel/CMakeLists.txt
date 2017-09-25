cmake_minimum_required(VERSION 3.9)

set(INSTALL_DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

#Bring the headers, such as Stdel.h into the project

include_directories(${MY_INCLUDE})
link_directories(${LIB_DIRS})

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "*.cpp")

#Generate the shared library from the sources
add_library(Spectral SHARED ${SOURCES})
target_link_libraries(Spectral ${LIBS})

install(TARGETS Spectral DESTINATION ${INSTALL_DESTINATION})


