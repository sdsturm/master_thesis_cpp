add_executable(production_num_int_sommerfeld_identity
    production_num_int_sommerfeld_identity.cpp)
target_link_libraries(production_num_int_sommerfeld_identity
    mthesis-lib
    Boost::timer)
