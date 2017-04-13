set(MOM5_DIR "${SRC_DIR}/mom5")
set(FMS_DIR "${SRC_DIR}/shared")

# MOM5 ocean only
file(GLOB_RECURSE MOM5_SOLO_SOURCE LIST_DIRECTORIES false
  ${MOM5_DIR}/drivers/*.[fF]90)

add_executable(MOM5_solo ${MOM5_SOLO_SOURCE} ${MOM5_SOURCE} ${OCEANSHARED_SOURCE})

target_include_directories(MOM5_solo PRIVATE
$<TARGET_FILE_DIR:FMS> $<TARGET_FILE_DIR:ocean> ${MOM5_DIR}/ocean_core ${FMS_DIR}/include/ ${SRC_DIR}/ocean_shared)

target_link_libraries(MOM5_solo FMS ocean ${NETCDF_FLIBS})
# target_link_libraries(MOM5 FMS ${NETCDF_FLIBS})
