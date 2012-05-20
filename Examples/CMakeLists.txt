#
# Insert examples here.
#
option(LSTK_USE_VTK "Build visualization helper tools." OFF)


if (LSTK_USE_VTK)
  find_package(VTK REQUIRED)

  set(VERSION_MIN "5.9.20110419")
  if (${VTK_VERSION} VERSION_LESS ${VERSION_MIN})
    message(ERROR " LSTK_USE_VTK requires VTK version ${VERSION_MIN} or newer but the current version is ${VTK_VERSION}")
  endif()

  include(${VTK_USE_FILE})
 
  option(LSTK_USE_AIM "Allow AIM output of lesion segmentation results" OFF)

  if (LSTK_USE_AIM)
    find_package( Boost REQUIRED )
    include_directories(${Boost_INCLUDE_DIR})

    find_package( LibXml2 REQUIRED )
    if( LIBXML2_FOUND )
      include_directories(${LIBXML2_INCLUDE_DIR})
    endif()
  endif()

  # Send LSTK_USE_AIM as a preprocessor def
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/LesionSegmentationConfig.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/LesionSegmentationConfig.h)
  include_directories(${CMAKE_CURRENT_BINARY_DIR})

  if(LSTK_USE_AIM)
    set(ITK_VTK_LSTK_LIBRARIES
      ${ITKLesionSizingToolkit_LIBRARIES}
      ${VTK_LIBRARIES} 
      ${LIBXML2_LIBRARIES} )
  else()
    set(ITK_VTK_LSTK_LIBRARIES
      ${ITKLesionSizingToolkit_LIBRARIES}
      ${VTK_LIBRARIES} )
  endif()

  add_executable( LesionSegmentation
    LesionSegmentation.cxx
    itkLesionSegmentationCommandLineProgressReporter.cxx )

  target_link_libraries( LesionSegmentation ${ITK_VTK_LSTK_LIBRARIES} )

endif() # LSTK_USE_VTK
