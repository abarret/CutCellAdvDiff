CCAD_ADD_EXAMPLE(
  TARGET_NAME
    "radial"
  OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/examples/radial"
  OUTPUT_NAME
    main2d
  EXAMPLE_GROUP
    examples
  SOURCES
    main.cpp QFcn.cpp RadialBoundaryCond.cpp SetLSValue.cpp
  LINK_TARGETS
    CCAD2d
  INPUT_FILES
    input2d
  )

