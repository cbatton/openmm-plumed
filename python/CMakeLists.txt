set(WRAP_FILE PlumedPluginWrapper.cpp)
set(MODULE_NAME openmmplumed)

# Find python
find_package(Python COMPONENTS Interpreter Development REQUIRED)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import mpi4py; print(mpi4py.get_include())"
    OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Execute SWIG to generate source code for the Python module.

add_custom_command(
    OUTPUT "${WRAP_FILE}"
    COMMAND "${SWIG_EXECUTABLE}"
        -python -c++
        -o "${WRAP_FILE}"
        -I"${OPENMM_DIR}/include"
        -I"${MPI_CXX_INCLUDE_PATH}"
        -I"${MPI4PY_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/plumedplugin.i"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/plumedplugin.i"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# Compile the Python module.

add_custom_target(PythonInstall DEPENDS "${WRAP_FILE}")
set(OPENMMPLUMED_HEADER_DIR "${CMAKE_SOURCE_DIR}/openmmapi/include")
set(OPENMMPLUMED_LIBRARY_DIR "${CMAKE_BINARY_DIR}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
add_custom_command(TARGET PythonInstall
    COMMAND "${PIP_EXECUTABLE}" install .
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
