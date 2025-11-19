find_package(Check QUIET)

if (NOT Check_FOUND)
  add_library(check_stub STATIC
    ${CMAKE_SOURCE_DIR}/third_party/check_stub/check_stub.c)
  target_include_directories(check_stub PUBLIC
    ${CMAKE_SOURCE_DIR}/third_party/check_stub)
  add_library(Check::check ALIAS check_stub)
  set(Check_FOUND TRUE)
endif()

set(EMWAVE_TEST_TMP_DIR "${CMAKE_BINARY_DIR}/test-tmp")
file(MAKE_DIRECTORY "${EMWAVE_TEST_TMP_DIR}")

function(add_emwave_check_test target)
  cmake_parse_arguments(ARG "" "" "SOURCES" ${ARGN})
  add_executable(${target} ${ARG_SOURCES})
  target_compile_features(${target} PRIVATE c_std_11)
  target_link_libraries(${target} PRIVATE emwave_core Check::check)
  if (MSVC)
    target_compile_definitions(${target} PRIVATE _CRT_SECURE_NO_WARNINGS _CRT_NONSTDC_NO_DEPRECATE)
  endif()
  if (WIN32)
    target_compile_definitions(${target} PRIVATE WIN32_LEAN_AND_MEAN NOMINMAX)
  endif()
  add_test(NAME ${target} COMMAND ${target})
endfunction()

add_emwave_check_test(config_loader_tests SOURCES test_config_loader.c)
target_compile_definitions(config_loader_tests PRIVATE
  CONFIGS_DIR="${CMAKE_SOURCE_DIR}/configs"
  TEST_TMP_DIR="${EMWAVE_TEST_TMP_DIR}")

add_emwave_check_test(fdtd_core_tests SOURCES test_fdtd_core.c)

add_emwave_check_test(analysis_tests SOURCES test_analysis_alloc.c)

add_emwave_check_test(config_runtime_tests SOURCES test_config_runtime.c)

add_emwave_check_test(material_library_tests SOURCES test_material_library.c)
