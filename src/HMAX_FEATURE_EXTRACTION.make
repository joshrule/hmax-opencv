C++FLAGS = -O3 -fPIC
C++ = /usr//bin/g++ $(C++FLAGS) -Wall $(INCLUDE)
CFLAGS = -O
CC = /usr//bin/gcc $(CFLAGS) -Wall $(INCLUDE)

HOME = /usr#USER PUT PATH TO DIRECTORY WHERE OPENCV INCLUDE AND LIBRARY ARE LOCATED
LD_LIBRARY_PATH = $LD_LIBRARY_PATH
export LD_LIBRARY_PATH
INCLUDE = -I$(HOME)/include/opencv -I$(HOME)/include -I$(HOME)/include/opencv2
CURRENT_DIR = $(PWD)/
LIBPATH = -L$(HOME)/lib
LIBS = -lm #-ljpeg -ltiff -lpng 

OPENCV = -D__STDC_CONSTANT_MACROS -lopencv_core -lopencv_imgproc -lopencv_highgui -lopencv_ml -lopencv_video -lopencv_features2d -lopencv_calib3d -lopencv_objdetect -lopencv_contrib -lopencv_legacy -lopencv_flann

GCC = gcc -Wall

EXECUTABLES  = HMAX_FE

all: $(EXECUTABLES)

.c.o:; $(CC) $(CFLAGS) $< -c
.cpp.o:; $(C++) $(C++FLAGS) -D__STDC_CONSTANT_MACROS $< -c

OBJECTS =basicFeatureExtraction.o cHMAX.o Filters.o utilities.o 

HMAX_FE: $(OBJECTS)
	$(C++) -o $@ $(C++FLAGS) $(OBJECTS) $(LIBPATH) $(OPENCV) $(LIBS)
	
clean:; -rm -rf *.o
very_clean_vc: clean
	-rm -f $(LIBRARY) $(EXECUTABLES)

#Dependencies
$(EXECUTABLES): $(LIBRARY)

#Remember to add library path to LD_LIBRARY_PATH if not done already in bashr.
