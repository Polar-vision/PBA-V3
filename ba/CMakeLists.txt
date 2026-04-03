cmake_minimum_required(VERSION 3.10)
project(ba)

if (MSVC)
    message(STATUS "Using MSVC compiler")
    string(FIND "${CMAKE_CXX_FLAGS_RELEASE}" "/O2" O2_FOUND)
    if(O2_FOUND EQUAL -1)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
    endif()
    message(STATUS "CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")

    string(FIND "${CMAKE_C_FLAGS_RELEASE}" "/O2" O2_FOUND)
    if(O2_FOUND EQUAL -1)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2")
    endif()
    message(STATUS "CMAKE_C_FLAGS_RELEASE = ${CMAKE_C_FLAGS_RELEASE}")
    set(3RDPARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../3rdparty")
    
    if(CMAKE_CL_64)
        message(STATUS "MSVC 64-bit (x64)")
        # 确保 Vcpkg 的路径优先于系统路径
        set(VCPKG_INSTALLED_DIR "E:/zuo/vcpkg/installed") 
        set(VCPKG_TARGET_TRIPLET "x64-windows")
        list(INSERT CMAKE_PREFIX_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
        list(INSERT CMAKE_MODULE_PATH 0 "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share")
        # 查找 Abseil
        # set(absl_DIR "E:/zuo/vcpkg/installed/x64-windows/share/absl")
        # set(Eigen3_DIR "E:/zuo/vcpkg/installed/x64-windows/share/eigen3")
        # find_package(absl CONFIG REQUIRED)
        # find_package(absl REQUIRED 20240116) # 20250127
        # message("-- Found abseil version ${absl_VERSION}: ${absl_DIR}")
        
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
        # find_package(Ceres REQUIRED)
        # set(CERES_LIBS "E:/zuo/projects/ceres-solver/install/lib")
        # # CERES 设置
        # set(CERES_INCLUDE_DIR "E:/zuo/projects/ceres-solver/install/include")
        # # 查找Ceres和依赖项

        # find_package(absl CONFIG REQUIRED)
        # message(STATUS "Ceres include dirs: ${CERES_INCLUDE_DIR}")
        # message(STATUS "Ceres libraries: ${CERES_LIBS}")
    else()
        message(STATUS "MSVC 32-bit (Win32)")
        # Eigen 设置
        set(EIGEN3_INCLUDE_DIR "${3RDPARTY_DIR}/eigen3")
        set(SPECTRA_INCLUDE_DIR "${3RDPARTY_DIR}/spectra/include")
        set(Cholmod_LIBS_DIR "${3RDPARTY_DIR}/cholmod/x86/lib")
    endif()
endif()


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 源文件
set(SOURCES
    "src/BAExporter.cpp"
    # "src/BALPBAImp.cpp"
    # "src/BALSBAImp.cpp"
    # "src/dllmain.cpp"
    "src/IBA.cpp"
    # "src/IBALBA.cpp"
    "src/PBAImp.cpp"
    "src/SBAImp.cpp"
    # "src/stdafx.cpp"
)

# 头文件路径
set(BA_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
message(STATUS "BA_INCLUDE_DIR路径: ${BA_INCLUDE_DIR}")
# ----------------------------
# 1. 生成静态库（.lib）
# ----------------------------
add_library(ba STATIC ${SOURCES})
target_compile_definitions(ba PRIVATE BA_STATIC)  # 静态库不需要导出符号
target_include_directories(ba PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${SPECTRA_INCLUDE_DIR}
    ${BA_INCLUDE_DIR}
    ${3RDPARTY_DIR}
)
if(CMAKE_CL_64)
    target_include_directories(ba PRIVATE
        # ${CERES_INCLUDE_DIR}
        ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include  # 添加 Vcpkg 包含路径
    )
    # target_link_libraries(ba PRIVATE
    #     "${CERES_LIBS}/ceres.lib"
    # )

    target_link_libraries(ba PRIVATE ${CERES_LIBRARIES})
    # target_link_libraries(ba PRIVATE 
    #     Ceres::ceres 
    #     SuiteSparse::cholmod
    #     SuiteSparse::cxsparse
    #     absl::container  # 显式链接相同版本Abseil
    # )
    # 设置预编译头
    # target_precompile_headers(ba PRIVATE ${BA_INCLUDE_DIR}/stdafx.h)
else()
    target_link_libraries(ba PRIVATE
        "${Cholmod_LIBS_DIR}/libamd.lib"
        "${Cholmod_LIBS_DIR}/libcamd.lib"
        "${Cholmod_LIBS_DIR}/libccolamd.lib"
        "${Cholmod_LIBS_DIR}/libcholmod.lib"
        "${Cholmod_LIBS_DIR}/libcolamd.lib"
        "${Cholmod_LIBS_DIR}/libmetis_CHOLMOD.lib"
        "${Cholmod_LIBS_DIR}/libgoto2.lib"
    )
endif()

# ----------------------------
# 2. 生成动态库（.dll + .lib 导入库）
# ----------------------------
add_library(ba_shared SHARED ${SOURCES})
set_target_properties(ba_shared PROPERTIES LINK_FLAGS "/SAFESEH:NO")
target_compile_definitions(ba_shared PRIVATE BA_EXPORTS)  # 动态库需要导出符号
target_include_directories(ba_shared PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${SPECTRA_INCLUDE_DIR}
    ${BA_INCLUDE_DIR}
    ${3RDPARTY_DIR}
)
if(CMAKE_CL_64)
    target_include_directories(ba_shared PRIVATE
        # ${CERES_INCLUDE_DIR}
        ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include  # 添加 Vcpkg 包含路径
    )
    # target_link_libraries(ba_shared PRIVATE
    #     "${CERES_LIBS}/ceres.lib"
    # )
    target_link_libraries(ba_shared PRIVATE ${CERES_LIBRARIES})
    # target_link_libraries(ba_shared PRIVATE 
    #     Ceres::ceres 
    #     SuiteSparse::cholmod
    #     SuiteSparse::cxsparse
    #     absl::container  # 显式链接相同版本Abseil
    # )
    # target_precompile_headers(ba_shared PRIVATE ${BA_INCLUDE_DIR}/stdafx.h)
else()
    target_link_libraries(ba_shared PRIVATE
        "${Cholmod_LIBS_DIR}/libamd.lib"
        "${Cholmod_LIBS_DIR}/libcamd.lib"
        "${Cholmod_LIBS_DIR}/libccolamd.lib"
        "${Cholmod_LIBS_DIR}/libcholmod.lib"
        "${Cholmod_LIBS_DIR}/libcolamd.lib"
        "${Cholmod_LIBS_DIR}/libmetis_CHOLMOD.lib"
        "${Cholmod_LIBS_DIR}/libgoto2.lib"
    )
endif()



# 动态库的导出文件名（可选）
set_target_properties(ba_shared PROPERTIES OUTPUT_NAME "ba")

# 复制动态库依赖的 DLL（如 libgoto2.dll）
if(CMAKE_CL_64)
    # add_custom_command(TARGET ba_shared POST_BUILD
    #     COMMAND ${CMAKE_COMMAND} -E copy_if_different
    #     "E:/zuo/projects/ceres-solver/install/bin/ceres.dll"
    #     $<TARGET_FILE_DIR:ba_shared>
    # )
    
    add_custom_command(TARGET ba_shared POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_RUNTIME_DLLS:ba_shared> $<TARGET_FILE_DIR:ba_shared>
        COMMAND_EXPAND_LISTS
    )
else()
    add_custom_command(TARGET ba_shared POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${3RDPARTY_DIR}/cholmod/x86/dll/libgoto2.dll"
        $<TARGET_FILE_DIR:ba_shared>
    )
endif()
