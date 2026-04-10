add_executable(test_geometry unit/test_geometry.c)
target_link_libraries(test_geometry PRIVATE dcreator_core)
add_test(NAME geometry COMMAND test_geometry)

add_executable(test_kernel unit/test_kernel.c)
target_link_libraries(test_kernel PRIVATE dcreator_core)
add_test(NAME kernel COMMAND test_kernel)

# Regression test: run the tool on the legacy DCREATOR_Example edge
# dislocation input and compare against the known-good legacy output
# (captured in tests/regression/reference_legacy.fcc).
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/regression/input.fcc"
   AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/regression/reference_legacy.fcc")
    # Step 1: run dcreator
    add_test(NAME regression_legacy_edge_run
        COMMAND ${CMAKE_COMMAND} -E chdir
            ${CMAKE_CURRENT_SOURCE_DIR}/regression
            $<TARGET_FILE:dcreator> legacy_edge.dissparam)
    # Step 2: compare against the reference via compare.py
    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_FOUND)
        add_test(NAME regression_legacy_edge_compare
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare.py
                ${CMAKE_CURRENT_SOURCE_DIR}/regression/out.fcc
                ${CMAKE_CURRENT_SOURCE_DIR}/regression/reference_legacy.fcc
                1e-10)
        set_tests_properties(regression_legacy_edge_compare
            PROPERTIES DEPENDS regression_legacy_edge_run)
    endif()
endif()
