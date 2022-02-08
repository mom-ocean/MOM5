set(MOM5_DIR "${SRC_DIR}/mom5")

# MOM5 ocean only
set(MOM5_SOLO_SOURCE ${MOM5_DIR}/drivers/ocean_solo.F90)

add_executable(MOM5_solo ${MOM5_SOLO_SOURCE})

target_include_directories(MOM5_solo PRIVATE
  $<TARGET_FILE_DIR:FMS> 
  $<TARGET_FILE_DIR:ocean> 
  $<TARGET_FILE_DIR:version>)

target_link_libraries(MOM5_solo 
  ocean 
  version 
  FMS 
  ${NETCDF_FLIBS}
)
