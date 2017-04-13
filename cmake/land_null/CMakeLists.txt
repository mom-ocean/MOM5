set(LAND_DIR "${SRC_DIR}/land_null")
set(LAND_PARAM_DIR "${SRC_DIR}/land_param")
set(FMS_DIR "${SRC_DIR}/shared")

# common files for both builds
file(GLOB_RECURSE LAND_SOURCE LIST_DIRECTORIES false FOLLOW_SYMLINKS
  ${LAND_DIR}/*.[fF]90 ${LAND_PARAM_DIR}/*.[fF]90 )

add_library(land_null ${LAND_SOURCE})

target_include_directories(land_null PRIVATE $<TARGET_FILE_DIR:FMS>)

#  $<TARGET_FILE_DIR:FMS> ${FMS_DIR}/include/)

target_link_libraries(land_null FMS)
