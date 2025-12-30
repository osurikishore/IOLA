include_directories(${DA_UPDATE_ANALYSIS_SOURCE_DIR}/src)

set (EXECUTABLE "da_update_analysis.exe")

message (STATUS "Building with Fortran compiler ${CMAKE_Fortran_COMPILER}.")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
   set (CMAKE_Fortran_FLAGS_RELEASE -O3 -Wextra -std=f95 -fcheck=all -fbacktrace -Wall)
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
   set (CMAKE_Fortran_FLAGS_RELEASE -O3)
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")

message (STATUS "Release compiler flags = ${CMAKE_Fortran_FLAGS_RELEASE}")
message (STATUS "Debug compiler flags = ${CMAKE_Fortran_FLAGS_DEBUG}")

set (FFLAGS ${CMAKE_Fortran_FLAGS_RELEASE})

set (DDEFS "NETCDF" "SLINT_LIB_PATH" "SLINT_INC_PATH" "FFTW_LIB_PATH" "FFTW_INC_PATH")
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

if (EXISTS ${SLINT_LIB_PATH})
   message (STATUS "Building with SLINT_LIB in path = ${SLINT_LIB_PATH}")
else (EXISTS ${SLINT_LIB_PATH})
   message (FATAL_ERROR "The required SLINT_LIB path was either not specified by \
   	    the user and/or does not exist. Aborting!!!")
endif (EXISTS ${SLINT_LIB_PATH})
if (EXISTS ${SLINT_INC_PATH})
   message (STATUS "Building with SLINT_INC in path = ${SLINT_INC_PATH}")
else (EXISTS ${SLINT_INC_PATH})
   message (FATAL_ERROR "The required SLINT_INC path was either not specified by \
   	    the user and/or does not exist. Aborting!!!")
endif (EXISTS ${SLINT_INC_PATH})
include_directories(${SLINT_INC_PATH})

if (EXISTS ${FFTW_LIB_PATH})
   message (STATUS "Building with FFTW in path = ${FFTW_LIB_PATH}")
else (EXISTS ${FFTW_LIB_PATH})
   message (FATAL_ERROR "The required FFTW path was either not specified by \
   	    the user and/or does not exist. Aborting!!!")
endif (EXISTS ${FFTW_LIB_PATH})
if (EXISTS ${FFTW_INC_PATH})
   message (STATUS "Building with FFTW_INC in path = ${FFTW_INC_PATH}")
else (EXISTS ${FFTW_INC_PATH})
   message (FATAL_ERROR "The required FFTW_INC path was either not specified by \
   	    the user and/or does not exist. Aborting!!!")
endif (EXISTS ${FFTW_INC_PATH})
include_directories(${FFTW_INC_PATH})

file (GLOB ftn_objs *.F90)
set (FTN_OBJS ${ftn_objs})
add_executable(${EXECUTABLE} ${FTN_OBJS})

set (LLIBS "m" ${NETCDF}/lib/libnetcdff.dylib ${SLINT_LIB_PATH}/libslint.a ${FFTW_LIB_PATH}/libfftw3.a)
target_link_libraries (${EXECUTABLE} ${LLIBS})
set_target_properties (${EXECUTABLE} PROPERTIES FFLAGS "${FFLAGS}")

set (INSTALL_PATH ${DA_UPDATE_ANALYSIS_SOURCE_DIR}/bin)
set (DIRECTORY_PERMISSIONS "OWNER_WRITE" "OWNER_READ" "OWNER_EXECUTE" "GROUP_READ" "GROUP_EXECUTE" "WORLD_READ" "WORLD_EXECUTE")
install (DIRECTORY DESTINATION ${INSTALL_PATH} DIRECTORY_PERMISSIONS ${DIRECTORY_PERMISSIONS})
install (TARGETS ${EXECUTABLE} DESTINATION ${INSTALL_PATH})
set (DEMO_PATH ${DA_UPDATE_ANALYSIS_SOURCE_DIR}/demo)
file (GLOB DEMO_FILES "${DEMO_PATH}/*")
foreach (item ${DEMO_FILES})
   list (APPEND FILES_TO_DEPLOY "${item}") 
endforeach (item ${DEMO_FILES})
install (FILES ${FILES_TO_DEPLOY} DESTINATION ${INSTALL_PATH})

