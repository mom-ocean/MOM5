set(MOM5_DIR "${SRC_DIR}/mom5")
set(FMS_DIR "${SRC_DIR}/shared")

# common files for both builds
file(GLOB_RECURSE MOM5_SOURCE LIST_DIRECTORIES false FOLLOW_SYMLINKS
  ${MOM5_DIR}/*.[fF]90 ${MOM5_DIR}/*.c)

file(GLOB_RECURSE OCEANSHARED_SOURCE LIST_DIRECTORIES false FOLLOW_SYMLINKS
  ${SRC_DIR}/ocean_shared/*.[fF]90)

# MOM5 ocean only
file(GLOB_RECURSE MOM5_SOLO_SOURCE LIST_DIRECTORIES false
  ${MOM5_DIR}/drivers/*.[fF]90)

# Remove these from previous list
list(REMOVE_ITEM MOM5_SOURCE ${MOM5_SOLO_SOURCE})

add_library(ocean ${MOM5_SOURCE} ${OCEANSHARED_SOURCE})
target_include_directories(ocean PRIVATE
  $<TARGET_FILE_DIR:FMS> ${MOM5_DIR}/ocean_core ${MOM5_DIR}/ocean_param/gotm-4.0/include ${FMS_DIR}/include/ ${SRC_DIR}/ocean_shared)

target_link_libraries(ocean FMS)
