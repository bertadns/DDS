cmake_minimum_required(VERSION 3.14)

project(depinning LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

#set(CMAKE_AUTOUIC ON)
#set(CMAKE_AUTOMOC ON)
#set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set( CMAKE_CXX_FLAGS "-O3")
set( CMAKE_CXX_FLAGS "-Wall")

#set(CMAKE_CXX_FLAGS "-lboost_program_options")

#find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core)
#find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Core)

#find_package(Boost COMPONENTS program_options REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})
#target_link_libraries(btarget ${Boost_LIBRARIES} )

add_executable(depinning
  main.cpp
  general.cpp
  general.h
  screw.cpp
  screw.h
  parser.cpp
  parser.h
)


#target_link_libraries(depinning Qt${QT_VERSION_MAJOR}::Core)

install(TARGETS depinning
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
