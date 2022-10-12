#=======================================================================
#  MAKEFILE for modify_gravity
#=======================================================================
OBJDIR   = ./
EXEDIR   = ./
LDR         = gfortran
OPTS = -c 
# insert here if you need any flags
COMP_OPTIM = 
LDR_OPTIM =
DEBUG =
INCLIST =
LIBLIST = 


#	macro definition
.SUFFIXES: .o .f90

SRC = modify_gravity.f90 parameter.inc

SRC_GRAV = gravitytable_use_example.f90  gravitytable_subr.f90 

OBJ_GRAV = $(SRC_GRAV:.f90=.o) 

INC = gravitytable.inc

#--------------------  implicit compiling rules ---------------------------------
.f90.o: $(INC) $(SRC_GRAV)
	${LDR} ${DEBUG} ${OPTS} $(COMP_OPTIM) $(INCLIST) $*.f90



gtable: $(OBJ_GRAV)
	${LDR} ${LDR_OPTIM} ${DEBUG} -o gravy_table  $(INCLIST)  ${OBJ_GRAV} ${LIBLIST}  


clean:	
	rm *.o


#----------------------------- help ------------------------------------

help:
	@echo Type 'make' to generate executable
	@echo Type 'make clean' to remove object files, etc






