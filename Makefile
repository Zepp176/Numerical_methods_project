
INC_DIR := -I/home/victor/Documents/Etudes/Numerical_methods/project/lib_petsc/include
LIB_DIR := -L/home/victor/Documents/Etudes/Numerical_methods/project/lib_petsc/lib -Wl,-rpath=/home/victor/Documents/Etudes/Numerical_methods/project/lib_petsc/lib

LIB := -lpetsc

CXX_FLAGS := -O0 -Wall #-g

#Compilation
all :
	gcc -o project project.c poisson.c functions.c -lm -lmpi $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)

#Delete of executable files
clean :
	rm projet

#Delete of results
clean_txt :
	rm -vf results/*.txt results/P-* results/U-* results/V-* results/Reh-* results/Vtx-* results/Rehw-* results/Div-*

test :
	gcc -o utest unittest.c poisson.c functions.c -lm -lmpi $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)
