# Set and uncomment these variables
# GLPK_DIR=where glpk lives
# FTI2_DIR=where 4ti2 lives

export CFLAGS=-I $(GLPK_DIR)/include -I $(FTI2_DIR)/include 
export LDFLAGS=-L $(GLPK_DIR)/lib -L $(FTI2_DIR)/lib

all: module2 module3

module3: Py4ti2.cc setup.py
	python3 setup.py build_ext --inplace

module2: Py4ti2.cc setup.py
	python2 setup.py build_ext --inplace

install: install3

install3:
	python3 setup.py install --user

install2: module2
	python2 setup.py install --user

clean:
	rm Py4ti2*.so

