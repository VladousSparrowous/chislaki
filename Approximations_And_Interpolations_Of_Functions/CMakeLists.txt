add_library(polynom polynom.cpp)

target_include_directories(polynom PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_link_libraries(polynom PUBLIC math_ops)
