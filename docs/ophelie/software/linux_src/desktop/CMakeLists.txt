cmake_minimum_required(VERSION 3.21)
project(ophelie_desktop LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Ubuntu 24.04 system VTK 9.1 is built against Qt5 (VTK_QT_VERSION=5).
# VTK::mpi requires MPI::MPI_C — find MPI with C component before VTK.
find_package(MPI REQUIRED COMPONENTS C CXX)
find_package(Qt5 REQUIRED COMPONENTS Widgets OpenGL)
find_package(VTK REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  FiltersSources
  IOGeometry
  IOXML
  InteractionStyle
  RenderingCore
  RenderingOpenGL2
  GUISupportQt
)

add_executable(ophelie_desktop
  src/main.cpp
  src/MainWindow.cpp
  src/MainWindow.h
  src/ProjectModel.cpp
  src/ProjectModel.h
  src/SolverProcess.cpp
  src/SolverProcess.h
  src/VtkViewWidget.cpp
  src/VtkViewWidget.h
)

target_include_directories(ophelie_desktop PRIVATE src)
target_link_libraries(ophelie_desktop PRIVATE
  Qt5::Widgets
  Qt5::OpenGL
  MPI::MPI_C
  ${VTK_LIBRARIES}
)

vtk_module_autoinit(
  TARGETS ophelie_desktop
  MODULES ${VTK_LIBRARIES}
)

install(TARGETS ophelie_desktop RUNTIME DESTINATION bin)
