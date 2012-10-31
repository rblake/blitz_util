

find_library(blitz_lib blitz)
find_path(blitz_include blitz)
include_directories(${blitz_include})

option(BLITZ_RANGE_CHECKING "Turn on range checking in blitz.  Slows down code dramatically.  Turn off for production runs." ON)
if(BLITZ_RANGE_CHECKING)
  add_definitions(-DBZ_DEBUG)
endif(BLITZ_RANGE_CHECKING)

##################################################################

add_library(blitz_util
  math.cpp
)

target_link_libraries(blitz_util ${blitz_lib} lapack)

