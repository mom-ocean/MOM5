set(SIS_DIR "${SRC_DIR}/ice_sis")
set(ICE_PARAM_DIR "${SRC_DIR}/ice_param")

# common files for both builds
file(GLOB_RECURSE SIS_SOURCE LIST_DIRECTORIES false FOLLOW_SYMLINKS
  ${SIS_DIR}/*.[fF]90 ${ICE_PARAM_DIR}/*.[fF]90 )

add_library(ice_sis ${SIS_SOURCE} ${OCEANSHARED_SOURCE})
target_include_directories(ice_sis PRIVATE
  $<TARGET_FILE_DIR:FMS> ${FMS_DIR}/include/ )

target_link_libraries(ice_sis FMS)
