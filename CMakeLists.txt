cmake_minimum_required(VERSION 3.20)
project(EmWaveC C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Use MSYS2 config packages
find_package(SDL2 CONFIG REQUIRED)
find_package(SDL2_ttf CONFIG REQUIRED)

# Executable
add_executable(emwave src/main.c)

# Local includes if present
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include")
  target_include_directories(emwave PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
endif()

# Keep normal int main(...) (no SDL2main shim)
target_compile_definitions(emwave PRIVATE SDL_MAIN_HANDLED)

# Link ONLY SDL2 + SDL2_ttf (no SDL2main)
target_link_libraries(emwave PRIVATE SDL2::SDL2 SDL2_ttf::SDL2_ttf)

# Warnings
if (MSVC)
  target_compile_options(emwave PRIVATE /W4 /permissive-)
else()
  target_compile_options(emwave PRIVATE -Wall -Wextra -Wno-unused-parameter)
endif()

# Install to make a portable dist/
install(TARGETS emwave RUNTIME DESTINATION .)
if (EXISTS "${CMAKE_SOURCE_DIR}/assets")
  install(DIRECTORY assets/ DESTINATION assets)
endif()

# Copy a system font so TTF_OpenFont("DejaVuSans.ttf", ...) works on Windows
if (WIN32)
  add_custom_command(TARGET emwave POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "C:/Windows/Fonts/arial.ttf"
      "$<TARGET_FILE_DIR:emwave>/DejaVuSans.ttf")
endif()
