# =============================================================================
# MML_Ext - Miscellaneous extension algorithms
# =============================================================================

add_library(MML_Ext INTERFACE)

target_include_directories(MML_Ext INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(MML_Ext INTERFACE MML_Core MML_PackageCommon)

add_library(MML::Ext ALIAS MML_Ext)
