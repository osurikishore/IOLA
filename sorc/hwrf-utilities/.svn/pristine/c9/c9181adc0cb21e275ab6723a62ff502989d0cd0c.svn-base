include_directories(${TEMPDROP_SONDE_SOURCE_DIR}/src)

set (EXECUTABLE "tempdrop_sonde.exe")

message (STATUS "Building with Fortran compiler ${CMAKE_Fortran_COMPILER}.")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
   set (CMAKE_Fortran_FLAGS_RELEASE -O3 -g)
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
   set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -g -traceback")
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")

message (STATUS "Release compiler flags = ${CMAKE_Fortran_FLAGS_RELEASE}")
message (STATUS "Debug compiler flags = ${CMAKE_Fortran_FLAGS_DEBUG}")

set (FFLAGS ${CMAKE_Fortran_FLAGS_RELEASE})

set (DDEFS "NETCDF" "SONDE_LIB_PATH" "SLATEC_LIB_PATH" "NR_LIB_PATH")
foreach (item ${DDEFS})
   list (APPEND DEFINITIONS "-D${item}")
endforeach (item ${DDEFS})
add_definitions (${DEFINITIONS})

if (EXISTS ${NETCDF})
   message (STATUS "Building with NETCDF in path = ${NETCDF}")
else (EXISTS ${NETCDF})
   message (FATAL_ERROR "The required NETCDF path was either not specified by \
   	    the user and/or does not exist. Aborting!!!")
endif (EXISTS ${NETCDF})
include_directories (${NETCDF}/include)
find_library (NETCDF_LIB_PATH netcdf HINTS ${NETCDF})

if (EXISTS ${SONDE_LIB_PATH})
   message (STATUS "Building with SONDELIB in path = ${SONDE_LIB_PATH}")
else (EXISTS ${SONDE_LIB_PATH})
   message (FATAL_ERROR "The required SONDELIB path was either not specified by \
   	    the user and/or does not exist. Aborting!!!")
endif (EXISTS ${SONDE_LIB_PATH})
if (EXISTS ${SLATEC_LIB_PATH})
   message (STATUS "Building with SLATEC in path = ${SLATEC_LIB_PATH}")
else (EXISTS ${SLATEC_LIB_PATH})
   message (FATAL_ERROR "The required SLATEC path was either not specified by \
   	    the user and/or does not exist. Aborting!!!")
endif (EXISTS ${SLATEC_LIB_PATH})
include_directories(${NR_INC_PATH})
if (EXISTS ${NR_LIB_PATH})
   message (STATUS "Building with NR in path = ${NR_LIB_PATH}")
else (EXISTS ${NR_LIB_PATH})
   message (FATAL_ERROR "The required NR path was either not specified by \
   	    the user and/or does not exist. Aborting!!!")
endif (EXISTS ${NR_LIB_PATH})

file (GLOB ftn_objs *.F90)
set (FTN_OBJS ${ftn_objs})
add_executable(${EXECUTABLE} ${FTN_OBJS})

set (LLIBS "m" ${NETCDF}/lib/libnetcdf.a ${SONDE_LIB_PATH}/libsondelib.a ${SLATEC_LIB_PATH}/libslatec.a ${NR_LIB_PATH}/libnr.a)
target_link_libraries (${EXECUTABLE} ${LLIBS})
set_target_properties (${EXECUTABLE} PROPERTIES FFLAGS "${FFLAGS}")

set (INSTALL_PATH ${TEMPDROP_SONDE_SOURCE_DIR}/bin)
set (DIRECTORY_PERMISSIONS "OWNER_WRITE" "OWNER_READ" "OWNER_EXECUTE" "GROUP_READ" "GROUP_EXECUTE" "WORLD_READ" "WORLD_EXECUTE")
install (DIRECTORY DESTINATION ${INSTALL_PATH} DIRECTORY_PERMISSIONS ${DIRECTORY_PERMISSIONS})
install (TARGETS ${EXECUTABLE} DESTINATION ${INSTALL_PATH})
set (DEMO_PATH ${TEMPDROP_SONDE_SOURCE_DIR}/demo)
file (GLOB DEMO_FILES "${DEMO_PATH}/*")
foreach (item ${DEMO_FILES})
   list (APPEND FILES_TO_DEPLOY "${item}") 
endforeach (item ${DEMO_FILES})
install (FILES ${FILES_TO_DEPLOY} DESTINATION ${INSTALL_PATH})

