set(MOM5_DIR "${SRC_DIR}/mom5")
set(FMS_DIR "${SRC_DIR}/shared")

# common files for both builds
file(GLOB_RECURSE MOM5_SOURCE LIST_DIRECTORIES false FOLLOW_SYMLINKS
  ${MOM5_DIR}/ocean_blobs/*.[fF]90 
  ${MOM5_DIR}/ocean_bgc/*.[fF]90 
  ${MOM5_DIR}/ocean_core/*.[fF]90 
  ${MOM5_DIR}/ocean_diag/*.[fF]90 
  ${MOM5_DIR}/ocean_param/*.[fF]90 
  ${MOM5_DIR}/ocean_tracers/*.[fF]90 
  ${MOM5_DIR}/ocean_wave/*.[fF]90 
  ${SRC_DIR}/ocean_shared/generic_tracers/*.[fF]90
)

set(MOM5_DIR "${SRC_DIR}/mom5")

add_library(ocean ${MOM5_SOURCE} ${OCEANSHARED_SOURCE})
target_include_directories(ocean PRIVATE
  $<TARGET_FILE_DIR:FMS> 
  $<TARGET_FILE_DIR:version>
  ${MOM5_DIR}/ocean_param/gotm-4.0/include 
  ${MOM5_DIR}/ocean_core 
  ${FMS_DIR}/include/ 
)

target_link_libraries(ocean FMS version)
