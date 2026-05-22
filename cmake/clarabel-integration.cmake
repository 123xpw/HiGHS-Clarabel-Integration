# cmake/clarabel-integration.cmake
# Sets up Clarabel.cpp (Rust-based solver) integration for HiGHS.
# Included by the top-level CMakeLists.txt when USE_CLARABEL=ON.
# Provides targets:
#   libclarabel_c_static  -- INTERFACE lib linking to libclarabel_c.a
#   libclarabel_c_shared  -- INTERFACE lib linking to libclarabel_c.so
#   Eigen3::Eigen         -- INTERFACE lib with Eigen header-only include path

message(STATUS "Clarabel integration: enabled")

# --- Eigen3 (header-only, vendored submodule) ---
set(CLARABEL_EIGEN_DIR "${PROJECT_SOURCE_DIR}/extern/eigen")
if(NOT EXISTS "${CLARABEL_EIGEN_DIR}/Eigen/Core")
  message(FATAL_ERROR
    "Eigen submodule not found at ${CLARABEL_EIGEN_DIR}.\n"
    "Run:  git submodule update --init extern/eigen")
endif()

if(NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CLARABEL_EIGEN_DIR}")
  set(Eigen3_VERSION "3.4" CACHE INTERNAL "Eigen version from submodule")
  message(STATUS "Clarabel: Eigen3::Eigen target created (${CLARABEL_EIGEN_DIR})")
else()
  message(STATUS "Clarabel: Eigen3::Eigen target already exists")
endif()

# --- Clarabel.cpp Rust wrapper ---
set(CLARABEL_CPP_DIR "${PROJECT_SOURCE_DIR}/extern/Clarabel.cpp")
if(NOT EXISTS "${CLARABEL_CPP_DIR}/rust_wrapper/Cargo.toml")
  message(FATAL_ERROR
    "Clarabel.cpp submodule not found at ${CLARABEL_CPP_DIR}.\n"
    "Run:  git submodule update --init --recursive extern/Clarabel.cpp")
endif()

# Variables consumed by rust_wrapper/CMakeLists.txt
set(CLARABEL_PROJECT_VERSION "0.10.0")
set(CLARABEL_ROOT_DIR "${CLARABEL_CPP_DIR}")

# Make cargo visible: prefer user's installed cargo, fall back to ~/.cargo/bin
if(NOT DEFINED ENV{CARGO_HOME})
  set(CARGO_HOME_BIN "$ENV{HOME}/.cargo/bin")
  if(EXISTS "${CARGO_HOME_BIN}/cargo")
    set(ENV{PATH} "${CARGO_HOME_BIN}:$ENV{PATH}")
    message(STATUS "Clarabel: added ${CARGO_HOME_BIN} to PATH for cargo")
  endif()
endif()
find_program(CARGO_EXECUTABLE cargo HINTS "$ENV{HOME}/.cargo/bin" REQUIRED)
message(STATUS "Clarabel: cargo = ${CARGO_EXECUTABLE}")

# Add only the rust_wrapper; this creates:
#   libclarabel_c        (custom target: runs `cargo build [--release]`)
#   libclarabel_c_static (INTERFACE: links libclarabel_c.a)
#   libclarabel_c_shared (INTERFACE: links libclarabel_c.so/.dylib)
add_subdirectory(
  "${CLARABEL_CPP_DIR}/rust_wrapper"
  "${CMAKE_CURRENT_BINARY_DIR}/clarabel_rust_wrapper"
)

# Cache paths for use in other CMakeLists files
set(CLARABEL_CPP_INCLUDE_DIR "${CLARABEL_CPP_DIR}/include"
    CACHE INTERNAL "Clarabel C++ header directory")
set(CLARABEL_EIGEN_INCLUDE_DIR "${CLARABEL_EIGEN_DIR}"
    CACHE INTERNAL "Eigen header directory")

# Rust static libs on Linux need these system libraries at link time
if(UNIX AND NOT APPLE)
  set(CLARABEL_SYSTEM_LIBS dl pthread m CACHE INTERNAL "System libs for Rust static linking")
elseif(APPLE)
  set(CLARABEL_SYSTEM_LIBS "-framework Security" CACHE INTERNAL "")
else()
  set(CLARABEL_SYSTEM_LIBS "" CACHE INTERNAL "")
endif()

message(STATUS "Clarabel: C++ headers  -> ${CLARABEL_CPP_INCLUDE_DIR}")
message(STATUS "Clarabel: Eigen headers -> ${CLARABEL_EIGEN_INCLUDE_DIR}")
message(STATUS "Clarabel: system libs   -> ${CLARABEL_SYSTEM_LIBS}")
