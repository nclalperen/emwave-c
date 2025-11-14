add_executable(config_loader_tests test_config_loader.c)
target_link_libraries(config_loader_tests PRIVATE emwave_core)
target_compile_definitions(config_loader_tests PRIVATE CONFIGS_DIR="${CMAKE_SOURCE_DIR}/configs")

add_executable(fdtd_core_tests test_fdtd_core.c)
target_link_libraries(fdtd_core_tests PRIVATE emwave_core)

add_test(NAME config_loader_tests COMMAND config_loader_tests)
add_test(NAME fdtd_core_tests COMMAND fdtd_core_tests)
