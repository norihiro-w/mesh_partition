# Source files
IF(USE_METIS_SOURCE)
	include_directories(${METIS_PATH}/libmetis)
	include_directories(${METIS_PATH}/programs)
	
	set(METIS_SOURCES
	  ${METIS_PATH}/programs/mpmetis.c
	  ${METIS_PATH}/programs/cmdline_mpmetis.c 
	  ${METIS_PATH}/programs/io.c 
	  ${METIS_PATH}/programs/stat.c
	  metis_main.h
	)
	ADD_LIBRARY(METIS_lib STATIC ${METIS_SOURCES} )

	ADD_EXECUTABLE(partmesh-egs partition.cpp)

	target_link_libraries(partmesh-egs
	   Math
	   Mesh		
	   METIS_lib
       metis
	)
ELSE() 
        ADD_EXECUTABLE(partmesh-egs partition.cpp)
	target_link_libraries(partmesh-egs
		Math
		Mesh
	)	
ENDIF()



