enable_testing()
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

project(mthesis-lib-tests)

add_executable(${PROJECT_NAME}
    mthesis-lib-tests.cpp
    test_legendre_p_recurrence.cpp
    )

target_link_libraries(${PROJECT_NAME}
    mthesis-lib
    Boost::unit_test_framework
    )

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
