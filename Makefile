# Set and uncomment these variables here, or give them its value in the 
# make command line.
#
# GLPK_DIR=where glpk lives
# FTI2_DIR=where 4ti2 lives

.PHONY: all clean install module3 module2

ifndef GLPK_DIR
	$(info GLPK_DIR variable undefined)
endif
ifndef FTI2_DIR
	$(info FTI2_DIR variable undefined)
endif

export CFLAGS=-I $(GLPK_DIR)/include -I $(FTI2_DIR)/include 
export LDFLAGS=-L $(GLPK_DIR)/lib -L $(FTI2_DIR)/lib

ifneq ("$(MAKECMDGOALS)", "clean")
$(eval CONFIG_HEADER := $(FTI2_DIR)/include/4ti2/4ti2_config.h)
ifeq ($(shell grep -c 'define _4ti2_int32_t int32_t' $(CONFIG_HEADER)), 1) 
	export _4ti2_INT32=yes
endif
ifeq ($(shell grep -c 'define _4ti2_int64_t int64_t' $(CONFIG_HEADER)), 1)
	export _4ti2_INT64=yes
endif
ifeq ($(shell grep -c 'define _4ti2_HAVE_GMP' $(CONFIG_HEADER)), 1)
	export _4ti2_HAVE_GMP=yes
endif
endif

all: module3 module2

module3: Py4ti2.cc setup.py
	python3 setup.py build_ext --inplace

module2: Py4ti2.cc setup.py
	python2 setup.py build_ext --inplace

install: install3 install2

install3: module3
	python3 setup.py install --user --prefix=

install2: module2
	python2 setup.py install --user --prefix=

clean:
	rm -rf Py4ti2*.so build

