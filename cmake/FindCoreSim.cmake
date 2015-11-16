# CORE_INCLUDE_DIRS
# CORE_LIBRARIES

set(CORE_DIR $ENV{CORE_DIR})

find_package(SimModSuite REQUIRED)

find_path(PCU_INCLUDE_DIR PCU.h
          HINTS ${CORE_DIR}
          PATH_SUFFIXES include)

find_path(APF_INCLUDE_DIR apf.h
          HINTS ${CORE_DIR}
          PATH_SUFFIXES include)

find_path(APF_SIM_INCLUDE_DIR apfSIM.h
          HINTS ${CORE_DIR}
          PATH_SUFFIXES include)

find_path(GMI_INCLUDE_DIR gmi.h
          HINTS ${CORE_DIR}
          PATH_SUFFIXES include)

find_path(GMI_SIM_INCLUDE_DIR gmi_sim.h
          HINTS ${CORE_DIR}
          PATH_SUFFIXES include)

find_path(SPR_INCLUDE_DIR spr.h
          HINTS ${CORE_DIR}
          PATH_SUFFIXES include)

find_path(MA_INCLUDE_DIR ma.h
          HINTS ${CORE_DIR}
          PATH_SUFFIXES include)

find_path(MTH_INCLUDE_DIR mth.h
          HINTS ${CORE_DIR}
          PATH_SUFFIXES include)

find_path(MDS_INCLUDE_DIR apfMDS.h
          HINTS ${CORE_DIR}
          PATH_SUFFIXES include)


set(CORE_INCLUDE_DIRS ${PCU_INCLUDE_DIR} 
			 ${APF_INCLUDE_DIR} 
			 ${APF_SIM_INCLUDE_DIR} 
			 ${GMI_INCLUDE_DIR} 
			 ${GMI_SIM_INCLUDE_DIR}
			 ${SPR_INCLUDE_DIR}
			 ${MA_INCLUDE_DIR}
                         ${MTH_INCLUDE_DIR}
			 ${MDS_INCLUDE_DIR}
	                 ${SIMMODSUITE_INCLUDE_DIRS})

find_library(PCU_LIB pcu
             HINTS ${CORE_DIR}
             PATH_SUFFIXES lib)

find_library(APF_LIB apf
             HINTS ${CORE_DIR}
             PATH_SUFFIXES lib)

find_library(APF_SIM_LIB apf_sim
             HINTS ${CORE_DIR}
             PATH_SUFFIXES lib)

find_library(GMI_LIB gmi
             HINTS ${CORE_DIR}
             PATH_SUFFIXES lib)

find_library(GMI_SIM_LIB gmi_sim
             HINTS ${CORE_DIR}
             PATH_SUFFIXES lib)

find_library(SPR_LIB spr
             HINTS ${CORE_DIR}
             PATH_SUFFIXES lib)

find_library(MA_LIB ma
             HINTS ${CORE_DIR}
             PATH_SUFFIXES lib)

find_library(MTH_LIB mth
             HINTS ${CORE_DIR}
             PATH_SUFFIXES lib)

find_library(MDS_LIB mds
             HINTS ${CORE_DIR}
             PATH_SUFFIXES lib)

set(CORE_LIBRARIES ${MA_LIB}
                   ${MDS_LIB}
                   ${APF_SIM_LIB} 
                   ${APF_LIB} 
                   ${GMI_SIM_LIB} 
                   ${GMI_LIB} 
                   ${SPR_LIB} 
                   ${PCU_LIB}
                   ${MTH_LIB}
                   ${SIMMODSUITE_LIBS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CORE DEFAULT_MSG CORE_LIBRARIES CORE_INCLUDE_DIRS)

mark_as_advanced(CORE_LIBRARIES CORE_INCLUDE_DIRS)

