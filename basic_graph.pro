QMAKE_CXXFLAGS += -Werror

HEADERS = window.h \
	  hermite.h \
          spline.h

SOURCES = main.cpp \
	  window.cpp \
	  hermite.cpp \
	  spline.cpp

QT += widgets
