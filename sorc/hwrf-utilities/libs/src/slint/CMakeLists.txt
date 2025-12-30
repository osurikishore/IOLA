include_directories(${NWP_WMS_SOURCE_DIR}/src/slint)

message (STATUS "Building with Fortran compiler ${CMAKE_Fortran_COMPILER}.")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
   set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fbounds-check -g")
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
   set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -g")
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")

message (STATUS "Release compiler flags = ${CMAKE_Fortran_FLAGS_RELEASE}")
message (STATUS "Debug compiler flags = ${CMAKE_Fortran_FLAGS_DEBUG}")

set (LLIBS "m")

file (GLOB ftn_objs *.F90)
set (FTN_OBJS ${ftn_objs})

set (FFLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
add_library (slint ${FTN_OBJS})

install (TARGETS slint DESTINATION ${INSTALL_PATH})

