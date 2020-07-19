python_include=/software/SOURCES.MODIFIED/dedalus/dedalus/include/python3.5m/
boost_include=/home/yesil/Boost3.5/include/
numpy_include=/software/SOURCES.MODIFIED/dedalus/dedalus/lib/python3.5/site-packages/numpy-1.11.0-py3.5-linux-x86_64.egg/numpy/core/include/
eigen_include=/usr/include/eigen3
nsolver_include=/home/yesil/Documents/channelflowNew/channelflow-master/
dedsi_include=/home/yesil/deDSI_SRC/src/

_IDIR = $(python_include) $(boost_include) $(numpy_include) $(eigen_include) $(nsolver_include) $(dedsi_include)
IDIR=$(foreach d, $(_IDIR), -I$d)

python_lib=/software/SOURCES.MODIFIED/dedalus/dedalus/lib/python3.5/config-3.5m/
boost_lib=/home/yesil/Boost3.5/lib/
dedalus_lib=/software/SOURCES.MODIFIED/dedalus/dedalus/lib/
nsolver_lib=/home/yesil/Documents/channelflowNew/build/nsolver

_LDIR = $(python_lib) $(boost_lib) $(nsolver_lib) $(dedalus_lib)
LDIR=$(foreach d, $(_LDIR), -L$d)

_LIBS = :libpython3.5m.a pthread dl util m :libboost_python3.so.1.64.0 :libboost_numpy3.so.1.64.0 nsolver
LIBS=$(foreach d, $(_LIBS), -l$d)

LFLAGS = -Wl,-rpath=${boost_lib} -Wl,-rpath=${nsolver_lib} -Xlinker -export-dynamic

CXX = mpicxx
CFLAGS = -std=c++11 -Wall

EXECS = ./execs/

all: $(EXECS)findeigenvalsx $(EXECS)findsolnx $(EXECS)continuationx


SRC_DIR=./src/
SRC_FILES = $(wildcard $(SRC_DIR)*.cpp)

OBJ_DIR = ./objs_src/
OBJ_FILES = $(patsubst $(SRC_DIR)%cpp,$(OBJ_DIR)%o,$(SRC_FILES))

HEADERS = $(wildcard $(dedsi_include)*.h)


PROG_SRC_DIR=./progs/
PROG_SRC_FILES = $(wildcard $(PROGS_SRC_DIR)*.cpp)

PROG_OBJ_DIR = ./objs_progs/
PROG_OBJ_FILES = $(patsubst $(PROG_SRC_DIR)%cpp,$(PROG_OBJ_DIR)%.o,$(PROG_SRC_FILES))

$(OBJ_DIR)deDSI.o: $(SRC_DIR)deDSI.cpp $(HEADERS) $(nsolver_lib)
	$(CXX) $(CFLAGS) $(IDIR) -c $(SRC_DIR)deDSI.cpp -o $(OBJ_DIR)deDSI.o

$(OBJ_DIR)error.o: $(SRC_DIR)error.cpp $(HEADERS) $(nsolver_lib)
	$(CXX) $(CFLAGS) $(IDIR) -c $(SRC_DIR)error.cpp -o $(OBJ_DIR)error.o

$(OBJ_DIR)newton.o: $(SRC_DIR)newton.cpp $(HEADERS) $(nsolver_lib)
	$(CXX) $(CFLAGS) $(IDIR) -c $(SRC_DIR)newton.cpp -o $(OBJ_DIR)newton.o

$(PROG_OBJ_DIR)findsoln.o: $(PROG_SRC_DIR)findsoln.cpp 
	$(CXX) $(CFLAGS) $(IDIR) -c $(PROG_SRC_DIR)findsoln.cpp -o $(PROG_OBJ_DIR)findsoln.o

$(PROG_OBJ_DIR)findeigenvals.o: $(PROG_SRC_DIR)findeigenvals.cpp
	$(CXX) $(CFLAGS) $(IDIR) -c $(PROG_SRC_DIR)findeigenvals.cpp -o $(PROG_OBJ_DIR)findeigenvals.o

$(PROG_OBJ_DIR)continuation.o: $(PROG_SRC_DIR)continuation.cpp
	$(CXX) $(CFLAGS) $(IDIR) -c $(PROG_SRC_DIR)continuation.cpp -o $(PROG_OBJ_DIR)continuation.o


$(EXECS)findeigenvalsx: $(PROG_OBJ_DIR)findeigenvals.o $(OBJ_FILES)
	$(CXX) $^ $(CFLAGS) $(LDIR) $(LFLAGS) $(LIBS) -o $@

$(EXECS)findsolnx: $(PROG_OBJ_DIR)findsoln.o $(OBJ_FILES)
	$(CXX) $^ $(CFLAGS) $(LDIR) $(LFLAGS) $(LIBS) -o $@

$(EXECS)continuationx: $(PROG_OBJ_DIR)continuation.o $(OBJ_FILES)
	$(CXX) $^ $(CFLAGS) $(LDIR) $(LFLAGS) $(LIBS) -o $@

