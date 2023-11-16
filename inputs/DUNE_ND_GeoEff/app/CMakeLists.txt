find_package(ROOT)

if ( ${ROOT_FOUND} )

   set(EDEPSIM_ROOT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/EDepSimForMakeProject_20200424.root)
   set(EDEPSIM_ROOT_MAKE_PROJECT_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/makeEDepSimProject.C)

  # Use example EDepSim file to produce libraries and header files
  add_custom_command(OUTPUT ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libEDepSimEvents.so
  COMMAND root -b -q ${EDEPSIM_ROOT_FILE} ${EDEPSIM_ROOT_MAKE_PROJECT_SCRIPT}
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/EDepSimEvents/EDepSimEvents.so ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libEDepSimEvents.so
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/EDepSimEvents/*.h ${CMAKE_SOURCE_DIR}/include/
  DEPENDS ${EDEPSIM_ROOT_FILE} ${EDEPSIM_ROOT_MAKE_PROJECT_SCRIPT}
  )

  add_custom_target(compile_EDepSimEvents DEPENDS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libEDepSimEvents.so)

  add_library(EDepSimEvents SHARED IMPORTED)
  set_property(TARGET EDepSimEvents PROPERTY IMPORTED_LOCATION  ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libEDepSimEvents.so)
  add_dependencies(EDepSimEvents compile_EDepSimEvents)


  include_directories(${ROOT_INCLUDE_DIRS})
  add_executable(runGeoEffEDepSimBeamDir runGeoEffEDepSimBeamDir.cpp)
  add_executable(runGeoEffFDEvtSim runGeoEffFDEvtSim.cpp)
  add_executable(runGeoEffFDEvtSim_hadronplot runGeoEffFDEvtSim_hadronplot.cpp)

  target_link_libraries(runGeoEffEDepSimBeamDir geoEff EDepSimEvents ${ROOT_LIBRARIES})
  target_link_libraries(runGeoEffFDEvtSim geoEff ${ROOT_LIBRARIES})
  target_link_libraries(runGeoEffFDEvtSim_hadronplot geoEff ${ROOT_LIBRARIES})

else ()

  Message("Couldn't find either ROOT. Skipping app. Library will still compile."
)
endif()
