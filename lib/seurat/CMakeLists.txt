#----------------------------------------------------------------------------
# USER! : SET THE NAME OF YOUR LIBRARY
# (Replace 'MyLib' by your own library name)

#############################
SET ( LIBRARY_NAME   seurat  )
#############################
set(CMAKE_BUILD_TYPE Debug)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# CREATES A USER OPTION IN CMAKE
OPTION ( BUILD_${LIBRARY_NAME}  "Build ${LIBRARY_NAME} library ?" ON)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
IF ( BUILD_${LIBRARY_NAME} )
#----------------------------------------------------------------------------

  #----------------------------------------------------------------------------
  # BUILD LIBRARY
  #----------------------------------------------------------------------------

  #----------------------------------------------------------------------------
  # LIBRARY HEADERS (TO BE INSTALLED)
  # EITHER LIST ALL .h, *.txx IN CURRENT DIR USING NEXT LINE:

  FILE(GLOB_RECURSE ${LIBRARY_NAME}_HEADERS "*.h" "*.txx")
  
  # OR MANUALLY LIST YOUR HEADERS WITH NEXT COMMAND
  #  SET ( ${LIBRARY_NAME}_HEADERS
  #
  #      )
  #----------------------------------------------------------------------------

  #----------------------------------------------------------------------------
  # LIBRARY SOURCES (TO BE COMPILED)
  # EITHER LIST ALL .cxx, *.cpp, *.cc IN CURRENT DIR USING NEXT LINE:

  FILE(GLOB_RECURSE ${LIBRARY_NAME}_SOURCES *.cxx *.cpp *.cc *.h *.hxx)

  # OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  #  SET ( ${LIBRARY_NAME}_SOURCES 
  #   
  #      )
  #----------------------------------------------------------------------------

  #----------------------------------------------------------------------------
  # LIBRARY DEPENDENCIES (LIBRARIES TO LINK WITH)
  #
  # USER! : Uncomment the Libraries you need
  #
  SET ( ${LIBRARY_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
  #    ${VTK_LIBRARIES}
  #    ${ITK_LIBRARIES}
  #    ${GDCM_LIBRARIES}
  #    ${BOOST_LIBRARIES}
    ${OPENGL_LIBRARIES}

  # If this library must link against other libraries 
  # USER! : Add here any extra Library you need

      )
  #----------------------------------------------------------------------------

  #----------------------------------------------------------------------------
  # MACRO WHICH DOES ALL THE JOB : BUILD AND INSTALL

  # USER! : The default is to create a Dynamic Library.
  # if you need to create a static library
  # comment out the following line :

  CREA_ADD_LIBRARY( ${LIBRARY_NAME} )

  # and uncomment the 2 lines hereafter:

  # ADD_LIBRARY(${LIBRARY_NAME} STATIC  ${${LIBRARY_NAME}_SOURCES})
  # TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${${LIBRARY_NAME}_LINK_LIBRARIES} )

  #
  #----------------------------------------------------------------------------

  #---------------------------------------------------------------------------
ENDIF ( BUILD_${LIBRARY_NAME} )
