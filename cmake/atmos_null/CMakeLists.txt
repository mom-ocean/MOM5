set(ATMOS_DIR "${SRC_DIR}/atmos_null")
set(ATMOS_PARAM_DIR "${SRC_DIR}/atmos_param")
set(ATMOS_SHARED_DIR "${SRC_DIR}/atmos_shared")
set(FMS_DIR "${SRC_DIR}/shared")

# common files for both builds
file(GLOB_RECURSE ATMOS_SOURCE LIST_DIRECTORIES false FOLLOW_SYMLINKS
  ${ATMOS_DIR}/*.[fF]90 
  ${ATMOS_SHARED_DIR}/diag_integral/*.[fF]90 
  ${ATMOS_PARAM_DIR}/diag_integral/*.[fF]90 
  ${ATMOS_PARAM_DIR}/monin_obukhov/*.[fF]90 )

add_library(atmos_null ${ATMOS_SOURCE})
target_include_directories(atmos_null PRIVATE
  $<TARGET_FILE_DIR:FMS> 
  ${FMS_DIR}/include/ 
  ${FMS_DIR}/include/ 
  ${ATMOS_PARAM_DIR}/cosp)

target_link_libraries(atmos_null FMS)
