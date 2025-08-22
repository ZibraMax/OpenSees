#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_truss OBJECT)

#target_sources(OPS_Element_truss
#  PRIVATE

target_sources(OPS_Element
  PRIVATE
    OriHinge.cpp
  PUBLIC
    OriHinge.h
)




target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

