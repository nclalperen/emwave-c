# SPDX-License-Identifier: MIT

find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
  pkg_check_modules(PC_CHECK QUIET check)
endif()

find_path(Check_INCLUDE_DIR
  NAMES check.h
  HINTS ${PC_CHECK_INCLUDE_DIRS})

find_library(Check_LIBRARY
  NAMES check check_pic libcheck
  HINTS ${PC_CHECK_LIBRARY_DIRS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Check DEFAULT_MSG Check_LIBRARY Check_INCLUDE_DIR)

if (Check_FOUND AND NOT TARGET Check::check)
  add_library(Check::check UNKNOWN IMPORTED)
  set_target_properties(Check::check PROPERTIES
    IMPORTED_LOCATION "${Check_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${Check_INCLUDE_DIR}")
  if (PC_CHECK_LIBRARIES)
    set_target_properties(Check::check PROPERTIES
      INTERFACE_LINK_LIBRARIES "${PC_CHECK_LIBRARIES}")
  endif()
endif()
