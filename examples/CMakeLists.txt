add_executable(example_gpof_synth_signal example_gpof_synth_signal.cpp)
target_link_libraries(example_gpof_synth_signal
    mthesis-lib
    ${GNUPLOT_IOSTREAM_LINK_LIBS})

add_executable(example_tlgf_hs example_tlgf_hs.cpp)
target_link_libraries(example_tlgf_hs
    mthesis-lib
    ${GNUPLOT_IOSTREAM_LINK_LIBS})

add_executable(example_generic_gf example_generic_gf.cpp)
target_link_libraries(example_generic_gf mthesis-lib)

add_executable(example_dcim example_dcim.cpp)
target_link_libraries(example_dcim mthesis-lib)

add_executable(example_fmm example_fmm.cpp)
target_link_libraries(example_fmm
    mthesis-lib
    Boost::timer)

add_executable(example_legendre_p_recurrence example_legendre_p_recurrence.cpp)
target_link_libraries(example_legendre_p_recurrence
    mthesis-lib
    ${GNUPLOT_IOSTREAM_LINK_LIBS})








add_executable(test_01 test_01.cpp)
target_link_libraries(test_01 mthesis-lib)
