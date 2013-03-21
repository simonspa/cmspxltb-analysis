# Install script for directory: /home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/eudaq/trunk/bin/libeudaq.so")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libEutelescope.so.0.9.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libEutelescope.so.0.9"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libEutelescope.so"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      FILE(RPATH_CHECK
           FILE "${file}"
           RPATH "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/lib:/home/spanns/ilcsoft/v01-16-02/Marlin/v01-04/lib:/home/spanns/ilcsoft/v01-16-02/lcio/v02-03-03/lib:/home/spanns/ilcsoft/v01-16-02/java/usr/lib/x86_64-linux-gnu:/home/spanns/ilcsoft/v01-16-02/gear/v01-02-02/lib:/home/spanns/ilcsoft/v01-16-02/CLHEP/2.1.1.0/lib:/home/spanns/ilcsoft/v01-16-02/ilcutil/v01-00/lib:/home/spanns/ilcsoft/v01-16-02/MarlinUtil/v01-05-03/lib:/home/spanns/ilcsoft/v01-16-02/CED/v01-07/lib:/home/spanns/ilcsoft/v01-16-02/RAIDA/v01-06-02/lib:/home/spanns/ilcsoft/v01-16-02/root/5.34.04/lib:/home/spanns/ilcsoft/v01-16-02/lccd/v01-02/lib")
    ENDIF()
  ENDFOREACH()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES
    "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/lib/libEutelescope.so.0.9.0"
    "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/lib/libEutelescope.so.0.9"
    "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/lib/libEutelescope.so"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libEutelescope.so.0.9.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libEutelescope.so.0.9"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libEutelescope.so"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      FILE(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/spanns/ilcsoft/v01-16-02/Marlin/v01-04/lib:/home/spanns/ilcsoft/v01-16-02/lcio/v02-03-03/lib:/home/spanns/ilcsoft/v01-16-02/java/usr/lib/x86_64-linux-gnu:/home/spanns/ilcsoft/v01-16-02/gear/v01-02-02/lib:/home/spanns/ilcsoft/v01-16-02/CLHEP/2.1.1.0/lib:/home/spanns/ilcsoft/v01-16-02/ilcutil/v01-00/lib:/home/spanns/ilcsoft/v01-16-02/MarlinUtil/v01-05-03/lib:/home/spanns/ilcsoft/v01-16-02/CED/v01-07/lib:/home/spanns/ilcsoft/v01-16-02/RAIDA/v01-06-02/lib:/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/eudaq/trunk/bin:/home/spanns/ilcsoft/v01-16-02/root/5.34.04/lib:/home/spanns/ilcsoft/v01-16-02/lccd/v01-02/lib:"
           NEW_RPATH "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/lib:/home/spanns/ilcsoft/v01-16-02/Marlin/v01-04/lib:/home/spanns/ilcsoft/v01-16-02/lcio/v02-03-03/lib:/home/spanns/ilcsoft/v01-16-02/java/usr/lib/x86_64-linux-gnu:/home/spanns/ilcsoft/v01-16-02/gear/v01-02-02/lib:/home/spanns/ilcsoft/v01-16-02/CLHEP/2.1.1.0/lib:/home/spanns/ilcsoft/v01-16-02/ilcutil/v01-00/lib:/home/spanns/ilcsoft/v01-16-02/MarlinUtil/v01-05-03/lib:/home/spanns/ilcsoft/v01-16-02/CED/v01-07/lib:/home/spanns/ilcsoft/v01-16-02/RAIDA/v01-06-02/lib:/home/spanns/ilcsoft/v01-16-02/root/5.34.04/lib:/home/spanns/ilcsoft/v01-16-02/lccd/v01-02/lib")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pede2lcio" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pede2lcio")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pede2lcio"
         RPATH "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/lib:/home/spanns/ilcsoft/v01-16-02/Marlin/v01-04/lib:/home/spanns/ilcsoft/v01-16-02/lcio/v02-03-03/lib:/home/spanns/ilcsoft/v01-16-02/java/usr/lib/x86_64-linux-gnu:/home/spanns/ilcsoft/v01-16-02/gear/v01-02-02/lib:/home/spanns/ilcsoft/v01-16-02/CLHEP/2.1.1.0/lib:/home/spanns/ilcsoft/v01-16-02/ilcutil/v01-00/lib:/home/spanns/ilcsoft/v01-16-02/MarlinUtil/v01-05-03/lib:/home/spanns/ilcsoft/v01-16-02/CED/v01-07/lib:/home/spanns/ilcsoft/v01-16-02/RAIDA/v01-06-02/lib:/home/spanns/ilcsoft/v01-16-02/root/5.34.04/lib:/home/spanns/ilcsoft/v01-16-02/lccd/v01-02/lib")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/bin/pede2lcio")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pede2lcio" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pede2lcio")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pede2lcio"
         OLD_RPATH "/home/spanns/ilcsoft/v01-16-02/Marlin/v01-04/lib:/home/spanns/ilcsoft/v01-16-02/lcio/v02-03-03/lib:/home/spanns/ilcsoft/v01-16-02/java/usr/lib/x86_64-linux-gnu:/home/spanns/ilcsoft/v01-16-02/gear/v01-02-02/lib:/home/spanns/ilcsoft/v01-16-02/CLHEP/2.1.1.0/lib:/home/spanns/ilcsoft/v01-16-02/ilcutil/v01-00/lib:/home/spanns/ilcsoft/v01-16-02/MarlinUtil/v01-05-03/lib:/home/spanns/ilcsoft/v01-16-02/CED/v01-07/lib:/home/spanns/ilcsoft/v01-16-02/RAIDA/v01-06-02/lib:/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/eudaq/trunk/bin:/home/spanns/ilcsoft/v01-16-02/root/5.34.04/lib:/home/spanns/ilcsoft/v01-16-02/lccd/v01-02/lib:/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/lib:"
         NEW_RPATH "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/lib:/home/spanns/ilcsoft/v01-16-02/Marlin/v01-04/lib:/home/spanns/ilcsoft/v01-16-02/lcio/v02-03-03/lib:/home/spanns/ilcsoft/v01-16-02/java/usr/lib/x86_64-linux-gnu:/home/spanns/ilcsoft/v01-16-02/gear/v01-02-02/lib:/home/spanns/ilcsoft/v01-16-02/CLHEP/2.1.1.0/lib:/home/spanns/ilcsoft/v01-16-02/ilcutil/v01-00/lib:/home/spanns/ilcsoft/v01-16-02/MarlinUtil/v01-05-03/lib:/home/spanns/ilcsoft/v01-16-02/CED/v01-07/lib:/home/spanns/ilcsoft/v01-16-02/RAIDA/v01-06-02/lib:/home/spanns/ilcsoft/v01-16-02/root/5.34.04/lib:/home/spanns/ilcsoft/v01-16-02/lccd/v01-02/lib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pede2lcio")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pedestalmerge" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pedestalmerge")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pedestalmerge"
         RPATH "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/lib:/home/spanns/ilcsoft/v01-16-02/Marlin/v01-04/lib:/home/spanns/ilcsoft/v01-16-02/lcio/v02-03-03/lib:/home/spanns/ilcsoft/v01-16-02/java/usr/lib/x86_64-linux-gnu:/home/spanns/ilcsoft/v01-16-02/gear/v01-02-02/lib:/home/spanns/ilcsoft/v01-16-02/CLHEP/2.1.1.0/lib:/home/spanns/ilcsoft/v01-16-02/ilcutil/v01-00/lib:/home/spanns/ilcsoft/v01-16-02/MarlinUtil/v01-05-03/lib:/home/spanns/ilcsoft/v01-16-02/CED/v01-07/lib:/home/spanns/ilcsoft/v01-16-02/RAIDA/v01-06-02/lib:/home/spanns/ilcsoft/v01-16-02/root/5.34.04/lib:/home/spanns/ilcsoft/v01-16-02/lccd/v01-02/lib")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/bin/pedestalmerge")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pedestalmerge" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pedestalmerge")
    FILE(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pedestalmerge"
         OLD_RPATH "/home/spanns/ilcsoft/v01-16-02/Marlin/v01-04/lib:/home/spanns/ilcsoft/v01-16-02/lcio/v02-03-03/lib:/home/spanns/ilcsoft/v01-16-02/java/usr/lib/x86_64-linux-gnu:/home/spanns/ilcsoft/v01-16-02/gear/v01-02-02/lib:/home/spanns/ilcsoft/v01-16-02/CLHEP/2.1.1.0/lib:/home/spanns/ilcsoft/v01-16-02/ilcutil/v01-00/lib:/home/spanns/ilcsoft/v01-16-02/MarlinUtil/v01-05-03/lib:/home/spanns/ilcsoft/v01-16-02/CED/v01-07/lib:/home/spanns/ilcsoft/v01-16-02/RAIDA/v01-06-02/lib:/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/eudaq/trunk/bin:/home/spanns/ilcsoft/v01-16-02/root/5.34.04/lib:/home/spanns/ilcsoft/v01-16-02/lccd/v01-02/lib:/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/lib:"
         NEW_RPATH "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/lib:/home/spanns/ilcsoft/v01-16-02/Marlin/v01-04/lib:/home/spanns/ilcsoft/v01-16-02/lcio/v02-03-03/lib:/home/spanns/ilcsoft/v01-16-02/java/usr/lib/x86_64-linux-gnu:/home/spanns/ilcsoft/v01-16-02/gear/v01-02-02/lib:/home/spanns/ilcsoft/v01-16-02/CLHEP/2.1.1.0/lib:/home/spanns/ilcsoft/v01-16-02/ilcutil/v01-00/lib:/home/spanns/ilcsoft/v01-16-02/MarlinUtil/v01-05-03/lib:/home/spanns/ilcsoft/v01-16-02/CED/v01-07/lib:/home/spanns/ilcsoft/v01-16-02/RAIDA/v01-06-02/lib:/home/spanns/ilcsoft/v01-16-02/root/5.34.04/lib:/home/spanns/ilcsoft/v01-16-02/lccd/v01-02/lib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/pedestalmerge")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/home/spanns/ilcsoft/v01-16-02/Eutelescope/HEAD/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
