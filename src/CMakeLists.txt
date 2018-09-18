cmake_minimum_required(VERSION 3.5)

project(HypoTestInv LANGUAGES CXX)

# Using ROOT's CONFIG mode search
# The Makefile method also included (non-explicitly, though): Foam Html MathMore Minuit XMLParser
find_package(ROOT CONFIG REQUIRED COMPONENTS RooStats RooFit RooFitCore HistFactory)


if(TARGET ROOT::Minuit2)
    message(STATUS "Minuit2 included in ROOT")
else()
    message(STATUS "Minuit2 missing from ROOT, using standalone edition")
    find_package(Minuit2 CONFIG REQUIRED)
endif()


# Set flags and things automagically (not ideal, but simplest for now)
include("${ROOT_USE_FILE}")
message(STATUS "Root use_file Path is: " ${ROOT_USE_FILE})


# This function adds a LinkDef and ROOT dictionary for a library
function(ADD_LIBRARY_LINK_DEF NAME)
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/RootDict/${NAME}Dict_LinkDef.h"
        "#ifdef __CINT__\n"
        "#pragma link C++ nestedclasses;\n"
        "#pragma link C++ nestedtypedefs;\n"
        "#pragma link C++ defined_in \"${CMAKE_CURRENT_SOURCE_DIR}/${NAME}.h\";\n"
        "#endif\n")

    root_generate_dictionary(
        G__${NAME}
        ${NAME}.h
      LINKDEF
        "${CMAKE_CURRENT_BINARY_DIR}/RootDict/${NAME}Dict_LinkDef.h"
        )

    add_library(${NAME} G__${NAME}.cxx ${NAME}.cxx)
    target_link_libraries(${NAME} INTERFACE ${ROOT_LIBRARIES})

    # Link to ROOT's Minuit 2 if available, otherwise standalone edition
    if(TARGET ROOT::Minuit2)
        target_link_libraries(${NAME} INTERFACE ROOT::Minuit2)
    else()
        target_link_libraries(${NAME} INTERFACE Minuit2::Minuit2)
    endif()
endfunction()

# Adding core libraries

add_library_link_def(WorkspaceCalculator)
add_library_link_def(ProfileLikelihoodTestStatEnhanced)
add_library_link_def(LimitCalculator)
add_library_link_def(LikelihoodCalculator)
add_library_link_def(PValueCalculator)
add_library_link_def(FitCalculator)


# Adding executables

add_executable(StandardHypoTestInv StandardHypoTestInv.cxx)
target_link_libraries(StandardHypoTestInv PUBLIC
    WorkspaceCalculator ProfileLikelihoodTestStatEnhanced LimitCalculator)

#add_executable(StandardFrequentistDiscovery StandardFrequentistDiscovery.cxx)
#target_link_libraries(StandardFrequentistDiscovery PUBLIC
#    WorkspaceCalculator ProfileLikelihoodTestStatEnhanced PValueCalculator)

#add_executable(LikelihoodScan LikelihoodScan.cxx)
#target_link_libraries(LikelihoodScan PUBLIC
#    WorkspaceCalculator ProfileLikelihoodTestStatEnhanced LikelihoodCalculator)

#add_executable(SimpleFit SimpleFit.cxx)
#target_link_libraries(SimpleFit PUBLIC
#    WorkspaceCalculator ProfileLikelihoodTestStatEnhanced FitCalculator)

