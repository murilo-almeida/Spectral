cmake_minimum_required(VERSION 3.9)

set(INSTALL_DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

set(CMAKE_INSTALL_RPATH ${MY_RPATH})

#Bring the headers, such as Stdel.h into the project

include_directories(${MY_INCLUDE})
link_directories(${LIB_DIRS})

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "*.cpp")

set(exec_file "ed2d")
#Generate the executable
add_executable(${exec_file} ${SOURCES})
target_link_libraries(${exec_file}Spectral ${LIBS})

get_target_property(MY_R ${exec_file} BUILD_RPATH)
message("BUILD_RPATH ${MY_R}")

get_target_property(MY_R ${exec_file} INSTALL_RPATH)
message("INSTALL_RPATH ${MY_R}")

install(TARGETS ${exec_file} DESTINATION ${INSTALL_DESTINATION})


