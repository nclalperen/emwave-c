find_package(Check QUIET)

if (NOT Check_FOUND)
  add_library(check_stub STATIC
    ${CMAKE_SOURCE_DIR}/third_party/check_stub/check_stub.c)
  target_include_directories(check_stub PUBLIC
    ${CMAKE_SOURCE_DIR}/third_party/check_stub)
  add_library(Check::check ALIAS check_stub)
  set(Check_FOUND TRUE)
endif()

function(add_emwave_check_test target)
  cmake_parse_arguments(ARG "" "" "SOURCES" ${ARGN})
  add_executable(${target} ${ARG_SOURCES})
  target_link_libraries(${target} PRIVATE emwave_core Check::check)
  target_compile_definitions(${target} PRIVATE _CRT_SECURE_NO_WARNINGS)
  add_test(NAME ${target} COMMAND ${target})
endfunction()

add_emwave_check_test(config_loader_tests SOURCES test_config_loader.c)
target_compile_definitions(config_loader_tests PRIVATE CONFIGS_DIR="${CMAKE_SOURCE_DIR}/configs")

add_emwave_check_test(fdtd_core_tests SOURCES test_fdtd_core.c)

add_emwave_check_test(analysis_tests SOURCES test_analysis_alloc.c)
