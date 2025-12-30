include_directories(${NR_SOURCE_DIR}/src)

message (STATUS "Building with Fortran compiler ${CMAKE_Fortran_COMPILER}.")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
   set (CMAKE_Fortran_FLAGS_RELEASE -O3 -fno-range-check -g)
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
   set (CMAKE_Fortran_FLAGS_RELEASE -O3 -g)
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")

message (STATUS "Release compiler flags = ${CMAKE_Fortran_FLAGS_RELEASE}")
message (STATUS "Debug compiler flags = ${CMAKE_Fortran_FLAGS_DEBUG}")

set (FFLAGS ${CMAKE_Fortran_FLAGS_RELEASE})

file (GLOB src *.f)
add_library (nr ${src})
target_compile_options (nr PRIVATE ${FFLAGS})

set (DIRECTORY_PERMISSIONS "OWNER_WRITE" "OWNER_READ" "OWNER_EXECUTE" "GROUP_READ" "GROUP_EXECUTE" "WORLD_READ" "WORLD_EXECUTE")
set (INSTALL_PATH ${NR_SOURCE_DIR}/lib)
install (DIRECTORY DESTINATION ${INSTALL_PATH} DIRECTORY_PERMISSIONS ${DIRECTORY_PERMISSIONS})
install (TARGETS nr DESTINATION ${INSTALL_PATH})

