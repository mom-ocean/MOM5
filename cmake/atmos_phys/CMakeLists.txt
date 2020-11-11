set(ATMOS_PARAM_DIR "${SRC_DIR}/atmos_param")
set(ATMOS_SHARED_DIR "${SRC_DIR}/atmos_shared")
set(FMS_DIR "${SRC_DIR}/shared")

# common files for both builds
file(GLOB_RECURSE ATMOS_SOURCE LIST_DIRECTORIES false FOLLOW_SYMLINKS
  ${ATMOS_PARAM_DIR}/*.[fF]90 ${ATMOS_SHARED_DIR}/*.[fF]90 )

add_library(atmos_phys ${ATMOS_SOURCE})
target_include_directories(atmos_phys PRIVATE
  $<TARGET_FILE_DIR:FMS> ${FMS_DIR}/include/ ${ATMOS_PARAM_DIR}/cosp)

target_link_libraries(atmos_phys FMS)
