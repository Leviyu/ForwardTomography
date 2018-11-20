# This is my big Makefile which is used to compile all my c code


TOOL_DIR = ../01_cpp_lib
TOOL_CPP = $(wildcard ../01_cpp_lib/*.cpp)
TOOL_OBJ = $(TOOL_CPP:.cpp=.o)

TOMO_DIR = ../02_forward_tomo
TOMO_CPP = $(wildcard ../02_forward_tomo/*.cpp)
TOMO_OBJ = $(TOMO_CPP:.cpp=.o)

#ESF_DIR = ../../c_lib/c02_empirical_wavelet
#ESF_CPP = $(wildcard ../../c_lib/c02_empirical_wavelet/*.c)
#ESF_OBJ = $(ESF_CPP:.c=.o)

SACHOME = /usr/local/sac
SOURCE = ./main.cpp
#SOURCE = ./main_checker.cpp
#SOURCE = ./main_read_orignal_model.cpp
OBJ = $(SOURCE:.cpp=.o)

CC = c++ 
DEBUG = 
LIBS = libHongyucpp.a 

CXXFLAGS =  -c ${DEBUG} -I$(TOOL_DIR)  \
		 -I$(SACHOME)/include -I./  		\
		  -I$(TOMO_DIR)  -L$(TOOL_DIR)  \
		  -std=c++11
		  

LFLAGS =  ${DEBUG} -m64 -L$(TOOL_DIR) -L$(SACHOME)/lib -L./  -L. \
		 -I$(TOOL_DIR)  -I. -lsacio -lm  -L$(TOMO_DIR)  
		 #-L$(ESF_DIR) 


all: forward_tomography


forward_tomography: $(OBJ) libHongyucpp.a $(TOOL_DIR)/hongyulibcpp.h $(TOMO_DIR)/forward_tomography.h
	$(CC) -o $@ $^  $(LFLAGS)

$(OBJ): $(SOURCE) 
	$(CC) -c -o $@ $^ $(CXXFLAGS)

libHongyucpp.a: $(TOOL_OBJ)  $(TOMO_OBJ)
	ar cr $@ $^

$(TOOL_DIR)/%.o: $(TOOL_DIR)/%.cpp $(TOOL_DIR)/%.h
	$(CC) -c -o $@ $^ $(CXXFLAGS) 

$(TOMO_DIR)/%.o: $(TOMO_DIR)/%.cpp $(TOMO_DIR)/%.h
	$(CC) -c -o $@ $^ $(CXXFLAGS) 


clean :
	\rm *.o *~ forward_tomography *.a $(TOOL_OBJ) $(TOMO_OBJ) 
