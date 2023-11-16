find_package(Eigen3)

add_library(geoEff SHARED geoEff.cpp)
target_link_libraries(geoEff Eigen3::Eigen)

pybind11_add_module(pyGeoEff pyGeoEff.cpp)
target_link_libraries(pyGeoEff PRIVATE geoEff)
