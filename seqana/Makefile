include ../Make.helper

TFMEXECUTABLES = tfm_index_construct.x tfm_index_invert.x 

all: $(TFMEXECUTABLES)

#TFM EXECUTABLES
tfm_index_construct.x: include/tfm_index.hpp lib/tfm_index_construct.cpp
	$(MY_CXX) -Wall -Wextra -g $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib lib/tfm_index_construct.cpp -o tfm_index_construct.x $(LIBS)

tfm_index_invert.x: include/tfm_index.hpp lib/tfm_index_invert.cpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib lib/tfm_index_invert.cpp -o tfm_index_invert.x $(LIBS)

#CLEANERS
clean:
	rm -f *.sdsl
	rm -f *.x
