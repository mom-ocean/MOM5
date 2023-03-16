# Coupled model driver
file(GLOB_RECURSE MOM5_ICCM_SOURCE LIST_DIRECTORIES false
  ${SRC_DIR}/coupler/*.[fF]90)

add_executable(MOM5_ICCM ${MOM5_ICCM_SOURCE})

target_include_directories(MOM5_ICCM PRIVATE
  $<TARGET_FILE_DIR:FMS>
  $<TARGET_FILE_DIR:ocean> 
  $<TARGET_FILE_DIR:ice_sis>
  $<TARGET_FILE_DIR:atmos_phys> 
  $<TARGET_FILE_DIR:atmos_fv> 
  $<TARGET_FILE_DIR:atmos_bg> 
  $<TARGET_FILE_DIR:land_lad> 
  $<TARGET_FILE_DIR:version> 
  ${SRC_DIR}/ocean_shared
)

target_link_libraries(MOM5_ICCM 
  FMS 
  ocean 
  ice_sis 
  atmos_phys 
  atmos_fv 
  atmos_bg
  land_lad 
  version
  ${NETCDF_FLIBS}
)
