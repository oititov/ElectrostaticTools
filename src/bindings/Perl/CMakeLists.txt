find_package(PerlLibs)
  if(PERLLIBS_FOUND OR ((DEFINED ET_PERL_INCLUDE_PATH) AND (DEFINED ET_PERL_LIBRARY)))
    message("Perl libraries found. Perl bindings will be built.")
    
    set(ET_PERL_INCLUDE_PATH ${PERL_INCLUDE_PATH} CACHE PATH "Perl include path.")
    set(ET_PERL_LIBRARY ${PERL_LIBRARY} CACHE PATH "Perl library.")
        
    include_directories(${ET_PERL_INCLUDE_PATH})
    
    SET(CMAKE_SWIG_FLAGS "")
    SET_SOURCE_FILES_PROPERTIES(SOURCE ../elec_tools.i PROPERTIES CPLUSPLUS ON)
    SWIG_ADD_MODULE(molecule_perl perl ../elec_tools.i )
    SWIG_LINK_LIBRARIES(molecule_perl ${ET_PERL_LIBRARY} ${LIBS})
    set_target_properties(molecule_perl PROPERTIES PREFIX "" OUTPUT_NAME "ElectrostaticTools")

    set(ET_PERL_LIBDIR ${PERL_SITELIB} CACHE PATH "Perl bindings destination.")
    
    install(TARGETS ${SWIG_MODULE_molecule_perl_REAL_NAME} DESTINATION ${ET_PERL_LIBDIR})
    install(FILES ${CMAKE_BINARY_DIR}/src/bindings/Perl/ElectrostaticTools.pm DESTINATION ${ET_PERL_LIBDIR})
  endif(PERLLIBS_FOUND OR ((DEFINED ET_PERL_INCLUDE_PATH) AND (DEFINED ET_PERL_LIBRARY)))
