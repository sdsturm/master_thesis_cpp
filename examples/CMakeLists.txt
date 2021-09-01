add_executable(example_01 example_01.cpp)
target_link_libraries(example_01 mthesis-lib)

add_executable(test_01 test_01.cpp)
target_link_libraries(test_01 mthesis-lib)

add_executable(gpof_synth_signal gpof_synth_signal.cpp)
target_link_libraries(gpof_synth_signal mthesis-lib)

add_executable(tlgf_hs tlgf_hs.cpp)
target_link_libraries(tlgf_hs mthesis-lib)
