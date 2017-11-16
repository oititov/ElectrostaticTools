find_package(PythonLibs)
find_package(PythonInterp)
  
  if( (PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND) OR ( (DEFINED ET_PYTHON_INCLUDE_PATH) AND (DEFINED ET_PYTHON_LIBRARY) ) )
    message("Python libraries found. Python bindings will be built.")
    
    set(ET_PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_DIRS} CACHE PATH "Python include path.")
    set(ET_PYTHON_LIBRARY ${PYTHON_LIBRARIES} CACHE PATH "Python library.")
    
    include_directories(${ET_PYTHON_INCLUDE_PATH})
    
    SET(CMAKE_SWIG_FLAGS "")
    SET_SOURCE_FILES_PROPERTIES(SOURCE ../elec_tools.i PROPERTIES CPLUSPLUS ON)
    SWIG_ADD_MODULE(molecule_python python ../elec_tools.i )
    SWIG_LINK_LIBRARIES(molecule_python ${ET_PYTHON_LIBRARY} ${LIBS})
    set_target_properties(_molecule_python PROPERTIES PREFIX "_" OUTPUT_NAME "ElectrostaticTools")

    if (NOT DEFINED ET_PYTHON_LIBDIR)
      if (PYTHON_VERSION_MAJOR LESS 3)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
      endif (PYTHON_VERSION_MAJOR LESS 3)
      if (PYTHON_VERSION_MAJOR GREATER 2)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print (get_python_lib())" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
      endif (PYTHON_VERSION_MAJOR GREATER 2) 
    endif (NOT DEFINED ET_PYTHON_LIBDIR)
    
    set(ET_PYTHON_LIBDIR ${PYTHON_SITE_PACKAGES} CACHE PATH "Python bindings destination.")
    
    install(TARGETS ${SWIG_MODULE_molecule_python_REAL_NAME} DESTINATION ${ET_PYTHON_LIBDIR})
    install(FILES ${CMAKE_BINARY_DIR}/src/bindings/Python/ElectrostaticTools.py DESTINATION ${ET_PYTHON_LIBDIR})
  endif( (PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND) OR ( (DEFINED ET_PYTHON_INCLUDE_PATH) AND (DEFINED ET_PYTHON_LIBRARY) ) )
