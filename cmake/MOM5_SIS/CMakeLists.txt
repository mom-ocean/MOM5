# Coupled model driver
file(GLOB_RECURSE MOM5_SIS_SOURCE LIST_DIRECTORIES false
  ${SRC_DIR}/coupler/*.[fF]90)

add_executable(MOM5_SIS ${MOM5_SIS_SOURCE})

target_include_directories(MOM5_SIS PRIVATE
$<TARGET_FILE_DIR:FMS>
$<TARGET_FILE_DIR:ocean> 
$<TARGET_FILE_DIR:ice_sis>
$<TARGET_FILE_DIR:atmos_null> 
$<TARGET_FILE_DIR:land_null> 
${MOM5_DIR}/ocean_core 
${FMS_DIR}/include/ 
${SRC_DIR}/ocean_shared)

target_link_libraries(MOM5_SIS FMS ocean ice_sis atmos_null land_null ${NETCDF_FLIBS})
# target_link_libraries(MOM5 FMS ${NETCDF_FLIBS})
