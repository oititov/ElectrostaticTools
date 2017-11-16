find_package(SWIG)

option(ET_PYTHON_BINDINGS "Build Python bindings." ON )
option(ET_PERL_BINDINGS "Build Perl bindings." ON )

if(SWIG_FOUND)
  message("Swig found. Will attempt to build script bindings.")

  if (CMAKE_MAJOR_VERSION GREATER 2)
  cmake_policy(PUSH)
  cmake_policy(SET CMP0054 OLD)
  endif (CMAKE_MAJOR_VERSION GREATER 2)
  
  include(${SWIG_USE_FILE})
    
  if (ET_PERL_BINDINGS)
    add_subdirectory(Perl)
  endif (ET_PERL_BINDINGS)

  if (ET_PYTHON_BINDINGS)
    add_subdirectory(Python)
  endif (ET_PYTHON_BINDINGS)
  
  if (CMAKE_MAJOR_VERSION GREATER 2)
  cmake_policy(POP)
  endif (CMAKE_MAJOR_VERSION GREATER 2)
  
endif(SWIG_FOUND)
