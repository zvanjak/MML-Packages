# =============================================================================
# MML_DocsApp - Documentation Demo Application
# =============================================================================

set(docs_core_target MML_Core)
set(docs_packages_target MML_Packages)

# Configure this directory directly to consume prebuilt libraries; the three
# paths below point to package include/, the parent of MML Core's mml/, and build-Release/.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    cmake_minimum_required(VERSION 3.14)
    project(MML_Docs_Consumer LANGUAGES CXX)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    set(MML_INSTALLED_PREFIX "" CACHE PATH "Installed MML-Packages release prefix")
    if(MML_INSTALLED_PREFIX)
        find_package(MMLPackages REQUIRED CONFIG PATHS "${MML_INSTALLED_PREFIX}" NO_DEFAULT_PATH)
        set(docs_core_target MML::Core)
        set(docs_packages_target MML::Packages)
    else()
        set(MML_PACKAGE_INCLUDE_DIR "" CACHE PATH "Public MML-Packages include directory")
        set(MML_CORE_INCLUDE_DIR "" CACHE PATH "Directory containing the MML Core mml headers")
        set(MML_PREBUILT_DIR "" CACHE PATH "Existing MML-Packages build directory")
        foreach(required_path MML_PACKAGE_INCLUDE_DIR MML_CORE_INCLUDE_DIR MML_PREBUILT_DIR)
            if(NOT IS_DIRECTORY "${${required_path}}")
                message(FATAL_ERROR "${required_path} must point to an existing directory")
            endif()
        endforeach()

        add_library(MML_Core INTERFACE)
        target_include_directories(MML_Core INTERFACE
            "${MML_CORE_INCLUDE_DIR}"
            "${MML_CORE_INCLUDE_DIR}/mml"
            "${MML_CORE_INCLUDE_DIR}/mml/single_header"
        )
        add_library(MML_Packages INTERFACE)
        target_compile_definitions(MML_Packages INTERFACE MML_STATIC)
        target_include_directories(MML_Packages INTERFACE
            "${MML_PACKAGE_INCLUDE_DIR}"
            "${MML_PACKAGE_INCLUDE_DIR}/optimization"
            "${MML_PACKAGE_INCLUDE_DIR}/optimization/interfaces"
            "${MML_PACKAGE_INCLUDE_DIR}/optimization/core"
            "${MML_PACKAGE_INCLUDE_DIR}/optimization/algorithms"
            "${MML_PACKAGE_INCLUDE_DIR}/optimization/algorithms/LP"
            "${MML_PACKAGE_INCLUDE_DIR}/fourier"
            "${MML_PACKAGE_INCLUDE_DIR}/pde"
            "${MML_PACKAGE_INCLUDE_DIR}/statistics"
            "${MML_PACKAGE_INCLUDE_DIR}/symbolic"
        )
        foreach(package Fourier Optimization PDE Statistics Symbolic)
            string(TOLOWER "${package}" package_dir)
            set(package_library "${MML_PREBUILT_DIR}/src/${package_dir}/Release/${CMAKE_STATIC_LIBRARY_PREFIX}MML_${package}${CMAKE_STATIC_LIBRARY_SUFFIX}")
            if(NOT EXISTS "${package_library}")
                message(FATAL_ERROR "Missing prebuilt library: ${package_library}")
            endif()
            target_link_libraries(MML_Packages INTERFACE "${package_library}")
        endforeach()
    endif()
endif()

# Collect sources from subdirectories
file(GLOB MAIN_SOURCES "*.cpp")
file(GLOB OPTIMIZATION_SOURCES "optimization/*.cpp")
file(GLOB PDE_SOURCES "pde/*.cpp")
file(GLOB SYMBOLIC_SOURCES "symbolic/*.cpp")
file(GLOB MML_EXT_SOURCES "mml_ext/*.cpp")
set(MML_EXT_DEMO_SOURCES ${MML_EXT_SOURCES})
set(MML_EXT_MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/mml_ext/docs_demo_mml_ext_main.cpp")
list(REMOVE_ITEM MML_EXT_DEMO_SOURCES "${MML_EXT_MAIN_SOURCE}")

add_executable(MML_DocsApp
    ${MAIN_SOURCES}
    ${OPTIMIZATION_SOURCES}
    ${PDE_SOURCES}
    ${SYMBOLIC_SOURCES}
    ${MML_EXT_DEMO_SOURCES}
)

target_link_libraries(MML_DocsApp PRIVATE
    ${docs_core_target}
    ${docs_packages_target}
)

if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    target_include_directories(MML_DocsApp PRIVATE
        ${CMAKE_SOURCE_DIR}/test_data
        ${CMAKE_SOURCE_DIR}/test_beds
    )
endif()

# Organize in IDE (Visual Studio folders)
source_group("Main" FILES ${MAIN_SOURCES})
source_group("Optimization" FILES ${OPTIMIZATION_SOURCES})
source_group("PDE" FILES ${PDE_SOURCES})
source_group("Symbolic" FILES ${SYMBOLIC_SOURCES})
source_group("MML_Ext" FILES ${MML_EXT_DEMO_SOURCES})

add_executable(MML_Ext_DocsApp
    ${MML_EXT_DEMO_SOURCES}
    ${MML_EXT_MAIN_SOURCE}
)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    target_link_libraries(MML_Ext_DocsApp PRIVATE ${docs_core_target} ${docs_packages_target})
else()
    target_link_libraries(MML_Ext_DocsApp PRIVATE MML_Ext)
endif()

if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    target_include_directories(MML_Ext_DocsApp PRIVATE
        ${CMAKE_SOURCE_DIR}/test_data
        ${CMAKE_SOURCE_DIR}/test_beds
    )
endif()

source_group("MML_Ext" FILES ${MML_EXT_DEMO_SOURCES} ${MML_EXT_MAIN_SOURCE})
