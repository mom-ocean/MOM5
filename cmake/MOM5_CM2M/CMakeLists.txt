# Coupled model driver
file(GLOB_RECURSE MOM5_CM2M_SOURCE LIST_DIRECTORIES false
  ${SRC_DIR}/coupler/*.[fF]90)

add_executable(MOM5_CM2M ${MOM5_CM2M_SOURCE})

target_include_directories(MOM5_CM2M PRIVATE
  $<TARGET_FILE_DIR:FMS>
  $<TARGET_FILE_DIR:ocean> 
  $<TARGET_FILE_DIR:ice_sis>
  $<TARGET_FILE_DIR:atmos_fv> 
  $<TARGET_FILE_DIR:atmos_phys> 
  $<TARGET_FILE_DIR:land_lad> 
  ${MOM5_DIR}/ocean_core 
  ${FMS_DIR}/include/ 
  ${SRC_DIR}/ocean_shared)

target_link_libraries(MOM5_CM2M FMS ocean ice_sis atmos_phys atmos_fv 
                                     land_lad ${NETCDF_FLIBS})
