add_subdirectory(src)

add_executable(main main.cpp)

target_link_libraries(main PRIVATE methods math_ops)

target_include_directories(main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)