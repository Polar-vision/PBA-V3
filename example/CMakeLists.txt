cmake_minimum_required(VERSION 3.10)
project(example)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:CONSOLE")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(3RDPARTY_DIR "E:/zuo/projects/backup/3rdparty")
# set(BA_LIBS "E:/zuo/projects/backup/ba/build/Release")
set(3RDPARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../3rdparty" CACHE PATH "Path to 3rdparty dependencies")
set(BA_LIBS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ba/build/Release" CACHE PATH "Path to BA build directory")
if(CMAKE_CL_64)
    # 确保 Vcpkg 的路径优先于系统路径
    set(VCPKG_INSTALLED_DIR "E:/zuo/vcpkg/installed") 
    set(VCPKG_TARGET_TRIPLET "x64-windows")
    list(INSERT CMAKE_PREFIX_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
    list(INSERT CMAKE_MODULE_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share")
    # 查找 Abseil
    # set(absl_DIR "E:/zuo/vcpkg/installed/x64-windows/share/absl")
    # find_package(absl REQUIRED 20240116) # 20250127
    # message("-- Found abseil version ${absl_VERSION}: ${absl_DIR}")
    # find_package(absl CONFIG REQUIRED)
    # find_package(abseil CONFIG REQUIRED 
    #     COMPONENTS container
    #     # VERSION 20250127  # 匹配错误提示的its_20250127命名空间
    #     PATHS "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share/abseil"
    #     NO_DEFAULT_PATH
    # )
    # 查找 Ceres
    set(cudss_DIR "E:/zuo/Program Files/NVIDIA cuDSS/v0.5/lib/12/cmake/cudss")
    set(Ceres_DIR "E:/zuo/projects/ceres-solver/install/lib/cmake/Ceres")
    find_package(Ceres REQUIRED)
else()
    set(Cholmod_LIBS_DIR "${3RDPARTY_DIR}/cholmod/x86/lib")
    message("-- Wind32 maga1 ${Cholmod_LIBS_DIR}")
endif()


# 源文件
set(SOURCES
    example.cpp
)

# 创建可执行程序
add_executable(example ${SOURCES})
set_target_properties(example PROPERTIES LINK_FLAGS "/SAFESEH:NO")
set(BA_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ba/include")

# 头文件路径
target_include_directories(example PRIVATE
    ${BA_INCLUDE_DIR}
    ${3RDPARTY_DIR}
)

target_link_libraries(example PRIVATE "${BA_LIBS_DIR}/ba.lib")
if(CMAKE_CL_64) # CMAKE_CL_64
    # target_link_libraries(example PRIVATE
    #     "${CERES_LIBS}/ceres.lib"
    # )
    target_include_directories(example PRIVATE
        # ${CERES_INCLUDE_DIR}
        ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include  # 添加 Vcpkg 包含路径
    )

    target_link_libraries(example PRIVATE ${CERES_LIBRARIES})
else()
    target_link_libraries(example PRIVATE
        "${Cholmod_LIBS_DIR}/libamd.lib"
        "${Cholmod_LIBS_DIR}/libcamd.lib"
        "${Cholmod_LIBS_DIR}/libccolamd.lib"
        "${Cholmod_LIBS_DIR}/libcholmod.lib"
        "${Cholmod_LIBS_DIR}/libcolamd.lib"
        "${Cholmod_LIBS_DIR}/libmetis_CHOLMOD.lib"
        "${Cholmod_LIBS_DIR}/libgoto2.lib"
    )
    message("-- Wind32 maga2 ${Cholmod_LIBS_DIR}")
endif()
# message("-- Wind32 XX ${Cholmod_LIBS}")

add_custom_command(TARGET example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${BA_LIBS_DIR}/ba.dll"
    $<TARGET_FILE_DIR:example>)

if(CMAKE_CL_64)
    # add_custom_command(TARGET example POST_BUILD
    #     COMMAND ${CMAKE_COMMAND} -E copy_if_different
    #     "E:/zuo/projects/ceres-solver/install/bin/ceres.dll"
    #     $<TARGET_FILE_DIR:example>
    # )
    add_custom_command(TARGET example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_RUNTIME_DLLS:example> $<TARGET_FILE_DIR:example>
        COMMAND_EXPAND_LISTS
    )
    # 复制Ceres及其所有依赖DLL
    # file(GLOB CERES_DLLS "${CERES_DIR}/bin/*.dll")
    # foreach(dll ${CERES_DLLS})
    #     add_custom_command(TARGET example POST_BUILD
    #         COMMAND ${CMAKE_COMMAND} -E copy_if_different
    #         ${dll}
    #         $<TARGET_FILE_DIR:example>)
    # endforeach()
else()
    add_custom_command(TARGET example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${BA_LIBS_DIR}/libgoto2.dll"
        $<TARGET_FILE_DIR:example>)
    message("-- Wind32 maga3 ${Cholmod_LIBS_DIR}")
endif()