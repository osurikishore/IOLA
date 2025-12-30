include_directories(${SONDELIB_SOURCE_DIR}/src)

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

set (CMAKE_Fortran_MODULE_DIRECTORY ${SONDELIB_SOURCE_DIR}/include)
set (FFLAGS ${CMAKE_Fortran_FLAGS_RELEASE})

file (GLOB ftn_objs *.f *.F90)
set (FTN_OBJS ${ftn_objs})
set (LLIBS "m")
add_library (sondelib ${FTN_OBJS})

set (DIRECTORY_PERMISSIONS "OWNER_WRITE" "OWNER_READ" "OWNER_EXECUTE" "GROUP_READ" "GROUP_EXECUTE" "WORLD_READ" "WORLD_EXECUTE")
set (INSTALL_PATH ${SONDELIB_SOURCE_DIR}/lib)
install (DIRECTORY DESTINATION ${INSTALL_PATH} DIRECTORY_PERMISSIONS ${DIRECTORY_PERMISSIONS})
install (TARGETS sondelib DESTINATION ${INSTALL_PATH})

