cmake_minimum_required(VERSION 2.8)

SET(CTEST_DO_SUBMIT ON)
set(CTEST_DO_MEMCHECK ON)
set(CTEST_DO_COVERAGE OFF) 

SET(CTEST_TEST_TYPE Nightly)
set(CTEST_PROJECT_NAME "Biotissue")
set(CTEST_NIGHTLY_START_TIME "01:00:00 EST")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=Biotissue")
set(CTEST_DROP_SITE_CDASH TRUE)

set(CTEST_SITE Jenga.scorec.rpi.edu)

set(CTEST_DASHBOARD_ROOT "/fasttmp/mersoj/develop/test/biotissue/cdash" )
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
#set(CTEST_BUILD_CONFIGURATION RelWithDebInfo)
set(CTEST_BUILD_CONFIGURATION Debug)
set(CTEST_BUILD_NAME  "linux-gcc-${CTEST_BUILD_CONFIGURATION}")
set(CTEST_BUILD_FLAGS -j4)

set(CTEST_SOURCE_NAME src)
set(CTEST_BINARY_NAME build)

find_program(CTEST_COMMAND NAMES ctest)
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "/fasttmp/mersoj/develop/install/openmpi/1.10.7/share/openmpi/openmpi-valgrind.supp")
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full")

set(SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
set(BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

# we set these variables because ctest complains if they don't get set
# we reset them for each branch later
set(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
set(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")


find_program(CTEST_GIT_COMMAND NAMES git)
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

set(COMMON_OPTIONS 
  "-DCMAKE_C_COMPILER:FILEPATH=mpicc"
  "-DCMAKE_CXX_COMPILER:FILEPATH=mpicxx"
  "-DBUILD_TESTS=ON"
  "-DLOGRUN=TRUE"
  "-DSCOREC_DIR=$ENV{DEVROOT}/install/core/Debug/lib/cmake/SCOREC"
  "-DCMAKE_CXX_FLAGS_DEBUG:STRING=-g -O0 -fprofile-arcs -ftest-coverage"
  "-DCMAKE_C_FLAGS_DEBUG:STRING=-g -O0 -fprofile-arcs -ftest-coverage"
  "-DCMAKE_EXE_LINKER_FLAGS_DEBUG:STRING=-g -O0 -fprofile-arcs -ftest-coverage"
  "-Dlas_DIR=$ENV{DEVROOT}/install/las/Debug/lib/cmake"
  "-Dlas_core_DIR=$ENV{DEVROOT}/install/las/Debug/lib/cmake"
  "-DWRITE_MICRO_ITER=0"
  "-DWRITE_MICRO_STEP=0"
  "-DVERBOSITY=0"
  )

SET(DEVELOP_CONFIGURE_OPTIONS "${COMMON_OPTIONS}"
  "-DCMAKE_PREFIX_PATH=$ENV{DEVROOT}/test/amsi/cdash/build/develop/amsi"
)
SET(MASTER_CONFIGURE_OPTIONS "${COMMON_OPTIONS}"
  "-DCMAKE_PREFIX_PATH=$ENV{DEVROOT}/test/amsi/cdash/build/master/amsi"
)

function(git_exec CMD ACTION)
  string(REPLACE " " ";" CMD2 "${CMD}")
  message("Running \"git ${CMD}\"")
  execute_process(COMMAND "${CTEST_GIT_COMMAND}" ${CMD2}
    WORKING_DIRECTORY "${CTEST_SOURCE_DIRECTORY}"
    RESULT_VARIABLE RETVAR)
  if(RETVAR)
    message(FATAL_ERROR "${ACTION} failed (code ${RETVAR})!")
  else()
    message("${ACTION} succeeded")
  endif()
endfunction(git_exec)


function(checkout_branch BRANCH_NAME)
  git_exec("checkout ${BRANCH_NAME}"
           "Checking out branch ${BRANCH_NAME}")
endfunction(checkout_branch)



function(run_tests BINARY_DIRECTORY SOURCE_DIRECTORY BRANCH_NAME CONFIGURE_OPTIONS)
  set(BRANCH_NAME "${BRANCH_NAME}")
  set(BUILD_DIRECTORY "${BINARY_DIRECTORY}/${BRANCH_NAME}")
  # These two variables need to be set for ctest to not complain!
  set(CTEST_BINARY_DIRECTORY "${BUILD_DIRECTORY}")
  set(CTEST_SOURCE_DIRECTORY "${SOURCE_DIRECTORY}")
  # set update to checkout and pull the requested branch
  checkout_branch(${BRANCH_NAME})
  ctest_empty_binary_directory("${BUILD_DIRECTORY}")
  # set subproject name we are operating on
  set_property(GLOBAL PROPERTY SubProject "${BRANCH_NAME}")
  set_property(GLOBAL PROPERTY Label "${BRANCH_NAME}")
  ctest_start("${CTEST_TEST_TYPE}" "${SOURCE_DIRECTORY}" "${BUILD_DIRECTORY}")
  ctest_update(SOURCE "${SOURCE_DIRECTORY}")
  ctest_configure(BUILD "${BUILD_DIRECTORY}" SOURCE "${SOURCE_DIRECTORY}" OPTIONS "${CONFIGURE_OPTIONS}")
  ctest_build(BUILD "${BUILD_DIRECTORY}" CONFIGURATION "${CTEST_BUILD_CONFIGURATION}")
  ctest_test(BUILD "${BUILD_DIRECTORY}")
  if(CTEST_DO_COVERAGE)
    ctest_coverage(BUILD "${BUILD_DIRECTORY}")
  endif(CTEST_DO_COVERAGE)
  if(CTEST_DO_MEMCHECK)
    ctest_memcheck(BUILD "${BUILD_DIRECTORY}")
  endif(CTEST_DO_MEMCHECK)
  if(CTEST_DO_SUBMIT)
    ctest_submit(PARTS Start Update Configure Build Test Coverage MemCheck
        RETRY_COUNT 4
        RETRY_DELAY 30
        RETURN_VALUE SUBMIT_ERROR)
    if(SUBMIT_ERROR)
      message(WARNING "Could not submit ${BRANCH_NAME} results to CDash (code ${SUBMIT_ERROR})!")
    else()
      message("Submitted ${BRANCH_NAME} results to CDash")
    endif()
  endif()
endfunction(run_tests)

if(CTEST_DO_SUBMIT)
  ctest_submit(FILES "${CTEST_SCRIPT_DIRECTORY}/Project.xml"
      RETURN_VALUE HAD_ERROR)
  if(HAD_ERROR)
    message(WARNING "Cannot submit Biotissue Project.xml!")
  endif()
endif()

run_tests( "${BINARY_DIRECTORY}" "${SOURCE_DIRECTORY}" develop "${DEVELOP_CONFIGURE_OPTIONS}")
run_tests( "${BINARY_DIRECTORY}" "${SOURCE_DIRECTORY}" master "${MASTER_CONFIGURE_OPTIONS}")
