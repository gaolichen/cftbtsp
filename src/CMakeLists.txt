# check whether c++11 is supported.
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)

function(enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()

# example usage
#enable_cxx_compiler_flag_if_supported("-Wall")
enable_cxx_compiler_flag_if_supported("-Werror")
#enable_cxx_compiler_flag_if_supported("-Wextra")
#enable_cxx_compiler_flag_if_supported("-pedantic")

SET(BTSP_LIB_SRCS 
    BootstrapRunner.cpp
    common.cpp
	ConformalBlock.cpp 
    CorrelatorSet.cpp
    CrossingEquations.cpp
    CftData.cpp
    gegenbauer_polynomial.cpp
    RandomAlgebricEquations.cpp
    Scalar4ptFunction.cpp
)

IF(WIN32)
	SET(BTSP_LIB_SRCS 
		${BTSP_LIB_SRCS}
        BootstrapRunner.h
        common.h
		ConformalBlock.h 
        CorrelatorSet.h
        CftData.h
        gegenbauer_polynomial.hpp
        Scalar4ptFunction.h
	)
	
	source_group("Header Files" FILES *.h)
ENDIF(WIN32)

SET(BTSP_MAIN_SRCS
	main.cpp)

enable_testing()

#I like to keep test files in a separate source directory called test
file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)

find_package(Boost 1.64 REQUIRED COMPONENTS
             unit_test_framework)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_library (btsplib SHARED ${BTSP_LIB_SRCS})
    add_executable (btsp ${BTSP_MAIN_SRCS})
    target_link_libraries(btsplib ${Boost_LIBRARIES})
    target_link_libraries(btsp ${Boost_LIBRARIES} btsplib)
    # the test code depends on c++ 11. Support of c++0x is not enough.
    # if the compiler does not support c++11, do not build test code.
    if(COMPILER_SUPPORTS_CXX11)
        add_executable(btsptest ${TEST_SRCS})
        target_link_libraries(btsptest ${Boost_LIBRARIES} btsplib)
        add_test(NAME btsptest COMMAND btsptest)
    endif()
endif()

