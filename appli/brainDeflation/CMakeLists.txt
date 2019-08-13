IF(BUILD_TEST)

#########################
SET ( EXE_NAME   testStereographic  )
#########################

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      testStereographic.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering
      vtkHybrid
      vtkCommon      

      SRep
      SRepInterpolation
      SRepVisualization
      SRepCortexRepLib

      itkvnl
      itkvnl_algo


      register
      MinimalGeodesic
  #    ${ITK_LIBRARIES}
  #    ${GDCM_LIBRARIES}
  #    ${BOOST_LIBRARIES}

  # USER! : Add here those agmonst the various (?) PROJECT LIBRARIES
  # you need for the current executable
  # (If you created only one Library, don't forget it !...)

  )
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# USER! : UNCOMMENT NEXT LINE IF YOU WANT A CONSOLE ON WINDOWS
# NB : YOUR MAIN MUST BE ADAPTED ALSO
#      SEE THE MACRO CREA_WXMAIN_WITH_CONSOLE IN creaWx.h
#SET(${EXE_NAME}_CONSOLE TRUE)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
CREA_ADD_EXECUTABLE( ${EXE_NAME} )




#########################
SET ( EXE_NAME   testStereographicProjection1  )
#########################

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      testStereographicProjection1.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering
      vtkHybrid
      vtkCommon

      SRep
      SRepInterpolation
      SRepVisualization
      SRepCortexRepLib

      itkvnl
      itkvnl_algo

      register
      MinimalGeodesic
  #    ${ITK_LIBRARIES}
  #    ${GDCM_LIBRARIES}
  #    ${BOOST_LIBRARIES}

  # USER! : Add here those agmonst the various (?) PROJECT LIBRARIES
  # you need for the current executable
  # (If you created only one Library, don't forget it !...)

  )
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# USER! : UNCOMMENT NEXT LINE IF YOU WANT A CONSOLE ON WINDOWS
# NB : YOUR MAIN MUST BE ADAPTED ALSO
#      SEE THE MACRO CREA_WXMAIN_WITH_CONSOLE IN creaWx.h
#SET(${EXE_NAME}_CONSOLE TRUE)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
CREA_ADD_EXECUTABLE( ${EXE_NAME} )

#########################
SET ( EXE_NAME   testStereographicAnimation  )
#########################

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      testStereographicAnimation.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering
      vtkHybrid
      vtkCommon
      itkvnl

      SRep
      SRepInterpolation
      SRepVisualization
  #    ${ITK_LIBRARIES}
  #    ${GDCM_LIBRARIES}
  #    ${BOOST_LIBRARIES}

  # USER! : Add here those agmonst the various (?) PROJECT LIBRARIES
  # you need for the current executable
  # (If you created only one Library, don't forget it !...)

  )
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# USER! : UNCOMMENT NEXT LINE IF YOU WANT A CONSOLE ON WINDOWS
# NB : YOUR MAIN MUST BE ADAPTED ALSO
#      SEE THE MACRO CREA_WXMAIN_WITH_CONSOLE IN creaWx.h
#SET(${EXE_NAME}_CONSOLE TRUE)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
CREA_ADD_EXECUTABLE( ${EXE_NAME} )

#########################
SET ( EXE_NAME   testRegistration  )
#########################

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      testRegistration.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering
      vtkHybrid
      vtkCommon

      SRep
      SRepInterpolation
      SRepVisualization
      SRepCortexRepLib
      MinimalGeodesic

      itkvnl
      itkvnl_algo

      register
  #    ${ITK_LIBRARIES}
  #    ${GDCM_LIBRARIES}
  #    ${BOOST_LIBRARIES}

  # USER! : Add here those agmonst the various (?) PROJECT LIBRARIES
  # you need for the current executable
  # (If you created only one Library, don't forget it !...)

  )
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# USER! : UNCOMMENT NEXT LINE IF YOU WANT A CONSOLE ON WINDOWS
# NB : YOUR MAIN MUST BE ADAPTED ALSO
#      SEE THE MACRO CREA_WXMAIN_WITH_CONSOLE IN creaWx.h
#SET(${EXE_NAME}_CONSOLE TRUE)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
CREA_ADD_EXECUTABLE( ${EXE_NAME} )


#########################
SET ( EXE_NAME   testRegistrationToSrep  )
#########################

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# EXECUTABLE SOURCES (TO BE COMPILED)
# EITHER LIST ALL .cxx, *.cpp, *.cc IN CURRENT DIR USING NEXT LINE:

#FILE(GLOB ${EXE_NAME}_SOURCES *.cxx *.cpp *.cc)

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      testRegistrationToSrep.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering
      vtkHybrid
      vtkCommon

      SRep
      SRepInterpolation
      SRepVisualization
      SRepCortexRepLib

      itkvnl
      itkvnl_algo


      register
      MinimalGeodesic
  #    ${ITK_LIBRARIES}
  #    ${GDCM_LIBRARIES}
  #    ${BOOST_LIBRARIES}

  # USER! : Add here those agmonst the various (?) PROJECT LIBRARIES
  # you need for the current executable
  # (If you created only one Library, don't forget it !...)

  )
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# USER! : UNCOMMENT NEXT LINE IF YOU WANT A CONSOLE ON WINDOWS
# NB : YOUR MAIN MUST BE ADAPTED ALSO
#      SEE THE MACRO CREA_WXMAIN_WITH_CONSOLE IN creaWx.h
#SET(${EXE_NAME}_CONSOLE TRUE)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
CREA_ADD_EXECUTABLE( ${EXE_NAME} )


#########################
SET ( EXE_NAME   testRegistrationConstruct  )
#########################

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# EXECUTABLE SOURCES (TO BE COMPILED)
# EITHER LIST ALL .cxx, *.cpp, *.cc IN CURRENT DIR USING NEXT LINE:

#FILE(GLOB ${EXE_NAME}_SOURCES *.cxx *.cpp *.cc)

# OR MANUALLY LIST YOUR FILES WITH NEXT COMMAND (WITHOUT EXTENSION)
  SET ( ${EXE_NAME}_SOURCES
      testRegistrationConstruct.cxx
    )
#----------------------------------------------------------------------------

INCLUDE_DIRECTORIES (

# USER! : Add here the directories holding the extra .h files you need
# e.g.
 ../../lib/mySampleLib

)

#----------------------------------------------------------------------------
# DEPENDENCIES (LIBRARIES TO LINK WITH)
SET ( ${EXE_NAME}_LINK_LIBRARIES
  #    ${crea_LIBRARIES}
  #    ${WXWIDGETS_LIBRARIES}
  #    ${KWWidgets_LIBRARIES}
      vtkRendering
      vtkFiltering
      vtkHybrid
      vtkCommon

      SRep
      SRepInterpolation
      SRepVisualization
      SRepCortexRepLib

      itkvnl
      itkvnl_algo


      register
      MinimalGeodesic

  # USER! : Add here those agmonst the various (?) PROJECT LIBRARIES
  # you need for the current executable
  # (If you created only one Library, don't forget it !...)

  )
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# USER! : UNCOMMENT NEXT LINE IF YOU WANT A CONSOLE ON WINDOWS
# NB : YOUR MAIN MUST BE ADAPTED ALSO
#      SEE THE MACRO CREA_WXMAIN_WITH_CONSOLE IN creaWx.h
#SET(${EXE_NAME}_CONSOLE TRUE)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# CREATES AND INSTALLS THE EXE
CREA_ADD_EXECUTABLE( ${EXE_NAME} )


ENDIF(BUILD_TEST)
