add_library(methods method.cpp)

target_include_directories(methods PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_link_libraries(methods PUBLIC math_ops)
