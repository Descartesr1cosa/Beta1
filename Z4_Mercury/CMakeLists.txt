cmake_minimum_required(VERSION 3.5)
project(Mercury)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile_commands.json" FORCE)

# Find MPI
find_package(MPI REQUIRED)

# Set Compile Options
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# include files directory
include_directories(${MPI_INCLUDE_PATH})
include_directories(../include)
include_directories(../Z4_Mercury/include)
include_directories(.)

# source files directiory
# Basic Files
aux_source_directory(../src/0_basic  DIRSRCS0)
aux_source_directory(../src/1_grid  DIRSRCS1)
aux_source_directory(../src/2_topology  DIRSRCS2)
aux_source_directory(../src/3_field  DIRSRCS3)

# Mercury Solver Files
aux_source_directory(./src  SRC)
# aux_source_directory(./src/1_boundary  SRC0)
# aux_source_directory(./src/2_operators  SRC1)
# aux_source_directory(./src/3_io  SRC2)
# aux_source_directory(./src/4_solver  SRC3)
# aux_source_directory(./src/physic  SRC4)
# aux_source_directory(./src/tools  SRC5)

# build the executable program
add_executable(Mercury main.cpp ${DIRSRCS0} ${DIRSRCS1} ${DIRSRCS2} ${DIRSRCS3} ${SRC})
#  ${SRC0} ${SRC1} ${SRC2} ${SRC3} ${SRC4} ${SRC5} ${SRC6})

# add the mpi linkings
target_link_libraries(Mercury PRIVATE ${MPI_LIBRARIES} stdc++fs)
