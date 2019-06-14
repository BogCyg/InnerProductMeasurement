cmake_minimum_required( VERSION 3.6.2 )

# For new project it is sufficient to change only its name in the following line
set( PROJECT_NAME InnerProd )

project( ${PROJECT_NAME} )

set( CMAKE_BUILD_TYPE Debug )
#set( CMAKE_BUILD_TYPE Release )


#[[ADD_DEFINITIONS(-g++ -O2 -fsigned-char -freg-struct-return -Wall -W -Wshadow -Wstrict-prototypes -Wpointer-arith -Wcast-qual -Winline -Werror)]]

if( WIN32 )
	set( CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS /W3 /GR /EHsc /std:c++17 /D_UNICODE /DUNICODE" )
	set( CMAKE_CXX_FLAGS_DEBUG "/MDd /Zi /Ob0 /Od /RTC1 /std:c++17 /D_UNICODE /DUNICODE" )
	message( "Win settings chosen..." )
elseif( UNIX )
	set( CMAKE_CXX_FLAGS "-std=c++17 -Wall" )
	set( CMAKE_CXX_FLAGS_DEBUG "-g -std=c++17 -Wall" )
	message( "Linux settings chosen..." )
elseif( MAC )
	set( CMAKE_CXX_FLAGS "-std=c++17 -Wall" )
	set( CMAKE_CXX_FLAGS_DEBUG "-g -std=c++17 -Wall" )
	message( "Mac settings chosen..." )
endif()


# Inform CMake where the header files are
include_directories( include )


# Automatically add all *.cpp files to the project
file ( GLOB SOURCES "./src/*.cpp" "./src/908/CPP/Src/*.cpp" )
#file ( GLOB_RECURSE SOURCES "./src/*.cpp" )
add_executable( ${PROJECT_NAME} ${SOURCES} ../ttmath/ttmathuint_x86_64_msvc.obj )

# You can either disable asm with the penalty of worse run-time performance (#define TTMATH_NOASM 1)
# or assembly ("C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\x86_amd64\ml64.exe" /c ttmathuint_x86_64_msvc.asm)
# and add the following object to the project. 
#add_executable( ${PROJECT_NAME} ../ttmath/ttmathuint_x86_64_msvc.obj )	this does not work, the above works

# Set the default project 
set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME} )


message( "CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}" )
