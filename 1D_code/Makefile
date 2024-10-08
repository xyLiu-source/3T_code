#======================================================================
# Local Makefile 
#======================================================================

# directories
MAIN_DIR  = $(PWD)
BUILD_DIR = $(MAIN_DIR)/build
SRC_DIR   = $(MAIN_DIR)
OBJ_DIR   = $(BUILD_DIR)/obj
DEP_DIR   = $(BUILD_DIR)/dep
INC_DIR   = $(EIGEN_ROOT_DIR) $(MKLROOT)/include
LIB_DIR   = $(MKLROOT)/lib/intel64
BIN_DIR   = $(PWD)

# extension
SRC_EXT = .cpp
OBJ_EXT = .o
DEP_EXT = .d

# source, objective, dependency and executable files for unit tests
SRC_FILE  = $(wildcard $(SRC_DIR)/*$(SRC_EXT))
OBJ_FILE  = $(SRC_FILE:$(SRC_DIR)/%$(SRC_EXT)=$(OBJ_DIR)/%$(OBJ_EXT))
DEP_FILE  = $(SRC_FILE:$(SRC_DIR)/%$(SRC_EXT)=$(DEP_DIR)/%$(DEP_EXT))
EXEC_FILE = $(BIN_DIR)/test

# compile and link settings
CXX      := g++
CXXFLAGS  = -std=c++17 -Wall -O3 -g
CXXFLAGS += -MT $@ -MMD -MP -MF $(DEP_DIR)/$*$(DEP_EXT)
CXXFLAGS += $(INC_DIR:%=-I %)
CXXFLAGS += -fopenmp
CXXFLAGS += -DWITHOUT_NUMPY -DMKL_LP64 
LDFLAGS   = $(LIB_DIR:%=-L %)

# single thread
LDFLAGS  += -m64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

# multiple threads
# LDFLAGS  += -m64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

all: $(EXEC_FILE) | build

.PHONY: build
build:
	@mkdir -p $(BUILD_DIR) $(OBJ_DIR) $(DEP_DIR)

-include $(DEP_FILE)

.PHONY: clean
clean:
	@rm -f $(EXEC_FILE)
	@rm -rf $(BUILD_DIR)

.PHONY: info
info:
	@echo "Source directory for test:\n$(DIR)\n"
	@echo "Build directory for test:\n$(BUILD_DIR)\n"
	@echo "Object directory for test:\n$(OBJ_DIR)\n"
	@echo "Dependency directory for test:\n$(DEP_DIR)\n"
	@echo "Source files for test:\n $(SRC_FILE:%=%\n)"
	@echo "Object files for test:\n $(OBJ_FILE:%=%\n)"
	@echo "Dependency files for test:\n $(DEP_FILE:%=%\n)"
	@echo "Exectuable file for test:\n $(EXEC_FILE:%=%\n)"

# linked object files to get an executable file
$(EXEC_FILE): $(OBJ_FILE) 
	@echo "[Link all objective files to executable file $@]"
	@$(CXX) -o $@ $^ $(LDFLAGS)

# compile all source files in test directory to get object and dependency files
$(OBJ_DIR)/%$(OBJ_EXT): $(SRC_DIR)/%$(SRC_EXT) | build
	@echo "[Compile $<]"
	@$(CXX) $(CXXFLAGS) -c $< -o $@ 

