TOP= ./
BIN= $(TOP)bin/
SRC= $(TOP)src/
OBJDIR= $(TOP)obj/
OUT= $(TOP)output/

SYSLIB= -lm

BOOST_ROOT=/opt/homebrew/Cellar/boost/1.86.0
BOOST_INC= ${BOOST_ROOT}/include

EIGEN_ROOT= YOUR_EIGEN_PATH

VPATH= $(SRC)
#=============================================================================
CC=g++#icpc#
CFLAGS= -g -O2 -std=c++14 -I$(BOOST_INC) -I$(EIGEN_ROOT) -Wall -Wextra
#=============================================================================
OBJ= $(addprefix $(OBJDIR), \
	main.o \
	lane-emden-solver.o \
	conservative_tides.o \
	dissipative_tides_bulk.o \
	ode_steppers.o \
	alpha-beta-funcs.o \
	shear_funcs.o \
	Li2.o \
	Li3.o \
	local_expansions.o \
	sim_params.o \
	outputfiles.o \
	)

DEPS= $(subst .o,.cpp,$(SRCS))
#=============================================================================
all: default.run
test: $(TEST)
#=============================================================================
%.run: $(OBJ)
	$(CC) $^ -o $(BIN)$@ $(SYSLIB) $(CFLAGS)
#=============================================================================
$(OBJDIR)%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<
#=============================================================================
