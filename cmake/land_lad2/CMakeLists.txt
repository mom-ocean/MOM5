set(LAND_DIR "${SRC_DIR}/land_lad2")
set(FMS_DIR "${SRC_DIR}/shared")

add_compile_definitions(USE_LOG_DIAG_FIELD_INFO)

# common files for both builds
file(GLOB_RECURSE LAND_SOURCE LIST_DIRECTORIES false FOLLOW_SYMLINKS
  ${LAND_DIR}/*.[fF]90 )

add_library(land_lad2 ${LAND_SOURCE})

target_include_directories(land_lad2 PRIVATE 
  $<TARGET_FILE_DIR:FMS> ${FMS_DIR}/include/)

target_link_libraries(land_lad2 FMS)
