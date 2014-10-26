
# DETECT_MSVC_VERSION
# -------------
#
# Detect current plaform and set following
# variables :
#   PLATFORM_COMPILER 
#   PLATFORM 

MACRO(DETECT_MSVC_VERSION)

  if(MSVC)
    if(MSVC_VERSION EQUAL 1400) # msvc8
      set(PLATFORM "msvc8")
      set(PLATFORM_COMPILER "vc80")
    endif()
    
    if(MSVC_VERSION EQUAL 1500) # msvc9
      set(PLATFORM "msvc9")
      set(PLATFORM_COMPILER "vc90")
    endif()
    
    if(MSVC_VERSION EQUAL 1600) # msvc10
      set(PLATFORM "msvc10")
      set(PLATFORM_COMPILER "vc100")
    endif()
	
	if(MSVC_VERSION EQUAL 1700) # msvc11
      set(PLATFORM "msvc11")
      set(PLATFORM_COMPILER "vc110")
    endif()
	
	if(MSVC_VERSION EQUAL 1800) # msvc12
      set(PLATFORM "msvc12")
      set(PLATFORM_COMPILER "vc120")
    endif()
  endif(MSVC)

ENDMACRO(DETECT_MSVC_VERSION)
