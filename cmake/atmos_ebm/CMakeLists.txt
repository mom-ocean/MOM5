set(ATMOS_DIR "${SRC_DIR}/atmos_ebm")
set(ATMOS_PARAM_DIR "${SRC_DIR}/atmos_param")
set(ATMOS_SHARED_DIR "${SRC_DIR}/atmos_shared")
set(ATMOS_COUPLED_DIR "${SRC_DIR}/atmos_coupled")
set(ATMOS_SPECTRAL_DIR "${SRC_DIR}/atmos_spectral")
set(FMS_DIR "${SRC_DIR}/shared")

add_compile_definitions(OVERLOAD_C8 OVERLOAD_C4 OVERLOAD_R4)

# Gather source files
file(GLOB_RECURSE ATMOS_SOURCE CONFIGURE_DEPENDS
  ${ATMOS_DIR}/*.[fF]90  ${ATMOS_PARAM_DIR}/*.[fF]90 
  ${ATMOS_SHARED_DIR}/*.[fF]90 ${ATMOS_SHARED_DIR}/*.[fF]90 
  ${ATMOS_SPECTRAL_DIR}/tools/*.[fF]90 )

# message(STATUS "ATMOS_SOURCE ${ATMOS_SOURCE}")

add_library(atmos_ebm ${ATMOS_SOURCE})
target_include_directories(atmos_ebm PRIVATE
  $<TARGET_FILE_DIR:FMS> ${FMS_DIR}/include/ ${ATMOS_PARAM_DIR}/cosp )

target_link_libraries(atmos_ebm FMS)
