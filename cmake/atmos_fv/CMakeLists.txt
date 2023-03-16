set(ATMOS_DIR "${SRC_DIR}/atmos_fv_dynamics")
set(ATMOS_COUPLED "${SRC_DIR}/atmos_coupled")
set(FMS_DIR "${SRC_DIR}/shared")

# common files for both builds
file(GLOB_RECURSE ATMOS_SOURCE LIST_DIRECTORIES false FOLLOW_SYMLINKS
  ${ATMOS_DIR}/tools/*.[fF]90 
  ${ATMOS_DIR}/model/*.[fF]90 
  ${ATMOS_DIR}/driver/coupled/*.[fF]90 
  ${ATMOS_COUPLED}/*.[fF]90 )

add_library(atmos_fv ${ATMOS_SOURCE})
target_include_directories(atmos_fv PRIVATE
  $<TARGET_FILE_DIR:FMS> 
  ${FMS_DIR}/include/ 
  ${FMS_DIR}/include/ 
  ${ATMOS_DIR}/model/ 
  $<TARGET_FILE_DIR:atmos_phys>)

target_link_libraries(atmos_fv atmos_phys FMS)
