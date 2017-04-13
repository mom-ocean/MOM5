set(MOM5_DIR "${SRC_DIR}/mom5")
set(FMS_DIR "${SRC_DIR}/shared")

# MOM5 ocean only
file(GLOB_RECURSE MOM5_SOLO_SOURCE LIST_DIRECTORIES false
  ${MOM5_DIR}/drivers/*.[fF]90)

add_executable(MOM5_SIS ${MOM5_SOLO_SOURCE} ${MOM5_SOURCE} ${OCEANSHARED_SOURCE})

target_include_directories(MOM5_SIS PRIVATE
$<TARGET_FILE_DIR:FMS> $<TARGET_FILE_DIR:ocean> $<TARGET_FILE_DIR:ice> ${MOM5_DIR}/ocean_core ${FMS_DIR}/include/ ${SRC_DIR}/ocean_shared)

target_link_libraries(MOM5_SIS FMS ocean ice ${NETCDF_FLIBS})
# target_link_libraries(MOM5 FMS ${NETCDF_FLIBS})
