add_executable(production_num_int_sommerfeld_identity
    production_num_int_sommerfeld_identity.cpp)
target_link_libraries(production_num_int_sommerfeld_identity
    mthesis-lib
    Boost::timer)

add_executable(thesis_complex_source_beam thesis_complex_source_beam.cpp)
target_link_libraries(thesis_complex_source_beam mthesis-lib)

add_executable(thesis_tlgf thesis_tlgf.cpp)
target_link_libraries(thesis_tlgf mthesis-lib)
