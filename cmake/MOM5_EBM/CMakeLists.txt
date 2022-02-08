# Coupled model driver
file(GLOB_RECURSE MOM5_EBM_SOURCE LIST_DIRECTORIES false
  ${SRC_DIR}/coupler/*.[fF]90)

add_compile_definitions(OVERLOAD_C8 OVERLOAD_C4 OVERLOAD_R4)

add_executable(MOM5_EBM ${MOM5_EBM_SOURCE})

target_include_directories(MOM5_EBM PRIVATE
  $<TARGET_FILE_DIR:FMS>
  $<TARGET_FILE_DIR:ocean> 
  $<TARGET_FILE_DIR:ice_sis>
  $<TARGET_FILE_DIR:atmos_ebm> 
  $<TARGET_FILE_DIR:land_lad> 
  $<TARGET_FILE_DIR:version> 
  ${MOM5_DIR}/ocean_core 
  ${FMS_DIR}/include/ 
  ${SRC_DIR}/ocean_shared
)

target_link_libraries(MOM5_EBM 
  FMS 
  ocean 
  ice_sis 
  atmos_ebm 
  land_lad 
  version
)
