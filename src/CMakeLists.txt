if(NOT USE_CLARABEL)
  return()
endif()

# Determine which HiGHS target to link against (FAST_BUILD vs legacy build)
if(TARGET highs)
  set(HIGHS_LINK_TARGET highs)
elseif(TARGET libhighs)
  set(HIGHS_LINK_TARGET libhighs)
else()
  message(FATAL_ERROR "verify_clarabel: no HiGHS library target found (highs / libhighs)")
endif()

# ── M2 + M3: data-mapping layer and solution write-back ─────────────────────
add_library(highs_clarabel_interface STATIC
  lp_data/HighsClarabelInterface.cpp
  lp_data/HighsClarabelSolution.cpp
)
target_compile_features(highs_clarabel_interface PRIVATE cxx_std_14)
target_link_libraries(highs_clarabel_interface PUBLIC
  ${HIGHS_LINK_TARGET}
  libclarabel_c_static
  Eigen3::Eigen
  ${CLARABEL_SYSTEM_LIBS}
)
target_include_directories(highs_clarabel_interface PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}          # so #include "HighsClarabelInterface.h" resolves
  ${CLARABEL_CPP_INCLUDE_DIR}
  ${CLARABEL_EIGEN_INCLUDE_DIR}
)

# ── M1 verification executable ───────────────────────────────────────────────
add_executable(verify_clarabel verify_clarabel.cpp)
target_compile_features(verify_clarabel PRIVATE cxx_std_14)
target_link_libraries(verify_clarabel PRIVATE
  highs_clarabel_interface
)
target_include_directories(verify_clarabel PRIVATE
  ${CLARABEL_CPP_INCLUDE_DIR}
  ${CLARABEL_EIGEN_INCLUDE_DIR}
)

message(STATUS "verify_clarabel: executable will link against ${HIGHS_LINK_TARGET} + Clarabel")

# ── M2 unit test ─────────────────────────────────────────────────────────────
add_executable(test_m2_mapping test_m2_mapping.cpp)
target_compile_features(test_m2_mapping PRIVATE cxx_std_14)
target_link_libraries(test_m2_mapping PRIVATE
  highs_clarabel_interface
)

# ── M3 unit test ─────────────────────────────────────────────────────────────
add_executable(test_m3_solution test_m3_solution.cpp)
target_compile_features(test_m3_solution PRIVATE cxx_std_14)
target_link_libraries(test_m3_solution PRIVATE
  highs_clarabel_interface
)

# ── M4 integration test ───────────────────────────────────────────────────────
add_executable(test_m4_integration test_m4_integration.cpp)
target_compile_features(test_m4_integration PRIVATE cxx_std_14)
target_link_libraries(test_m4_integration PRIVATE
  ${HIGHS_LINK_TARGET}
)

# ── M_MIP integration test ────────────────────────────────────────────────────
add_executable(test_m_mip_integration test_m_mip_integration.cpp)
target_compile_features(test_m_mip_integration PRIVATE cxx_std_14)
target_link_libraries(test_m_mip_integration PRIVATE
  ${HIGHS_LINK_TARGET}
)

