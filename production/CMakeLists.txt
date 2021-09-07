add_executable(thesis_numerical_integration_sommerfeld_identity
    thesis_numerical_integration_sommerfeld_identity.cpp)
target_link_libraries(thesis_numerical_integration_sommerfeld_identity
    mthesis-lib
    Boost::timer)

add_executable(thesis_complex_source_beam thesis_complex_source_beam.cpp)
target_link_libraries(thesis_complex_source_beam mthesis-lib)

add_executable(thesis_tlgf thesis_tlgf.cpp)
target_link_libraries(thesis_tlgf mthesis-lib)
