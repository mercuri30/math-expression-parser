cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(math_expression_parser)

include_directories(include)

add_executable(parser_test
    src/expression.cpp
    src/unit.cpp
)

target_link_libraries(parser_test cppunit)
