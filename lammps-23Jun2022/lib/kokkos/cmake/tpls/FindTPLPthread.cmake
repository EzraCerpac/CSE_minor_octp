# @HEADER
# ************************************************************************
#
#                        Kokkos v. 3.0
#       Copyright (2020) National Technology & Engineering
#               Solutions of Sandia, LLC (NTESS).
#
# Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Christian R. Trott (crtrott@sandia.gov)
#
# ************************************************************************
# @HEADER


SET(USE_THREADS FALSE)

IF(NOT TPL_Pthread_INCLUDE_DIRS AND NOT TPL_Pthread_LIBRARY_DIRS AND NOT TPL_Pthread_LIBRARIES)
  # Use CMake's Thread finder since it is a bit smarter in determining
  # whether pthreads is already built into the compiler and doesn't need
  # a library to link.
  FIND_PACKAGE(Threads)
  #If Threads found a copy of pthreads make sure it is one of the cases the tribits
  #tpl system cannot handle.
  IF(Threads_FOUND AND CMAKE_USE_PTHREADS_INIT)
    IF(CMAKE_THREAD_LIBS_INIT STREQUAL "" OR CMAKE_THREAD_LIBS_INIT STREQUAL "-pthread")
      SET(USE_THREADS TRUE)
    ENDIF()
  ENDIF()
ENDIF()

IF(USE_THREADS)
  SET(TPL_Pthread_INCLUDE_DIRS "")
  SET(TPL_Pthread_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
  SET(TPL_Pthread_LIBRARY_DIRS "")
ELSE()
  KOKKOS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( Pthread
    REQUIRED_HEADERS pthread.h
    REQUIRED_LIBS_NAMES pthread
      )
ENDIF()
