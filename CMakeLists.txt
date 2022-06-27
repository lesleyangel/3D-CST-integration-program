cmake_minimum_required(VERSION 3.0.0)              # CMAKE的最低版本
project(CST)                                       # 项目名

file(GLOB SRC_FILES                                # 将需要编译的文件弄到SRC_FILES里
    "${PROJECT_SOURCE_DIR}/*.h"
    "${PROJECT_SOURCE_DIR}/*.cpp")
add_executable(${CMAKE_PROJECT_NAME} ${SRC_FILES}) # 可执行文件包含的文件

target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/libdll/blas_win32_MTd.lib")      # 添加lib
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/libdll/lapack_win32_MTd.lib")    # 添加lib
target_compile_features(${CMAKE_PROJECT_NAME} PRIVATE cxx_std_11)   # 启用C++11

add_custom_command(
    TARGET ${CMAKE_PROJECT_NAME}
    POST_BUILD                                                     # 构架完成之后
    COMMAND ${CMAKE_COMMAND} -E copy_directory                     # 拷贝文件夹 
            "${PROJECT_SOURCE_DIR}/libdll/debug"                   # 把bebug文件夹下的内容
            "$<TARGET_FILE_DIR:${CMAKE_PROJECT_NAME}>"             # 拷贝到项目生成路径下
)
# # set(CMAKE_BUILD_TYPE DEBUG)
# # set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -g -ggdb")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# aux_source_directory(. DIR_SRCS)
# add_executable(TestProj ${DIR_SRCS})


# link_directories("${PROJECT_SOURCE_DIR}/libdll")
# link_libraries(main "${PROJECT_SOURCE_DIR}/libdll/blas_win32_MTd")
# link_libraries(main "${PROJECT_SOURCE_DIR}/libdll/lapack_win32_MTd")

# if(WINDOWS)
#     # Copy DLL files to output directory
#     if(CMAKE_CL_64)
#         add_custom_command(TARGET BarcodeReader POST_BUILD 
#         COMMAND ${CMAKE_COMMAND} -E copy_if_different
#         "${PROJECT_SOURCE_DIR}}/libdll/blas_win32_MTd.dll"              
#         $<TARGET_FILE_DIR:BarcodeReader>)
#         add_custom_command(TARGET BarcodeReader POST_BUILD 
#         COMMAND ${CMAKE_COMMAND} -E copy_if_different
#         "${PROJECT_SOURCE_DIR}}/libdll/lapack_win32_MTd.dll"              
#         $<TARGET_FILE_DIR:BarcodeReader>)
#     else()
#         add_custom_command(TARGET BarcodeReader POST_BUILD 
#         COMMAND ${CMAKE_COMMAND} -E copy_if_different
#         "${PROJECT_SOURCE_DIR}}/libdll/blas_win32_MTd.dll"              
#         $<TARGET_FILE_DIR:BarcodeReader>)
#         add_custom_command(TARGET BarcodeReader POST_BUILD 
#         COMMAND ${CMAKE_COMMAND} -E copy_if_different
#         "${PROJECT_SOURCE_DIR}}/libdll/lapack_win32_MTd.dll"              
#         $<TARGET_FILE_DIR:BarcodeReader>)
#     endif()
# endif()