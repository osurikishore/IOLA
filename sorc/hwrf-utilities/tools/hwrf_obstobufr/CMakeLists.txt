include_directories(${OBS_TO_BUFR_SOURCE_DIR}/src)

set (EXECUTABLE "obs_to_bufr.exe")

message (STATUS "Building with Fortran compiler ${CMAKE_Fortran_COMPILER}.")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
   set (CMAKE_Fortran_FLAGS_RELEASE -O3 -g)
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
   set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -g -traceback -C")
   set (CMAKE_Fortran_FLAGS_DOUBLE_REAL "-real-size 64")
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")

message (STATUS "Release compiler flags = ${CMAKE_Fortran_FLAGS_RELEASE}")
message (STATUS "Debug compiler flags = ${CMAKE_Fortran_FLAGS_DEBUG}")

set (FFLAGS ${CMAKE_Fortran_FLAGS_RELEASE} ${CMAKE_Fortran_FLAGS_DOUBLE_REAL})

set (DDEFS "BUFR_LIB_PATH" "SLATEC_LIB_PATH" "NR_LIB_PATH")
foreach (item ${DDEFS})
   list (APPEND DEFINITIONS "-D${item}")
endforeach (item ${DDEFS})
add_definitions (${DEFINITIONS})

file (GLOB ftn_objs *.F90)
set (FTN_OBJS ${ftn_objs})
add_executable(${EXECUTABLE} ${FTN_OBJS})

set (LLIBS "m" ${BUFR_LIB_PATH}/libbufr_i4r8.a ${SLATEC_LIB_PATH}/libslatec.a ${NR_LIB_PATH}/libnr.a) 
target_link_libraries (${EXECUTABLE} ${LLIBS})
set_target_properties (${EXECUTABLE} PROPERTIES FFLAGS "${FFLAGS}")

set (INSTALL_PATH ${OBS_TO_BUFR_SOURCE_DIR}/bin)
set (DIRECTORY_PERMISSIONS "OWNER_WRITE" "OWNER_READ" "OWNER_EXECUTE" "GROUP_READ" "GROUP_EXECUTE" "WORLD_READ" "WORLD_EXECUTE")
install (DIRECTORY DESTINATION ${INSTALL_PATH} DIRECTORY_PERMISSIONS ${DIRECTORY_PERMISSIONS})
install (TARGETS ${EXECUTABLE} DESTINATION ${INSTALL_PATH})
set (DEMO_PATH ${OBS_TO_BUFR_SOURCE_DIR}/demo)
file (GLOB DEMO_FILES "${DEMO_PATH}/*")
foreach (item ${DEMO_FILES})
   list (APPEND FILES_TO_DEPLOY "${item}") 
endforeach (item ${DEMO_FILES})
install (FILES ${FILES_TO_DEPLOY} DESTINATION ${INSTALL_PATH})

