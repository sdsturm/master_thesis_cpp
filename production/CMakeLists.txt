add_executable(numerical_integration_sommerfeld_identity
    numerical_integration_sommerfeld_identity.cpp)
target_link_libraries(numerical_integration_sommerfeld_identity
    mthesis-lib
    Boost::timer)
