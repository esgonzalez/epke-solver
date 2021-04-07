exec       = epke-run
cc         = g++
opt        = -g
modules    = pugi parareal epke utility
source_dir = $(addprefix src/,$(modules))
build_dir  = $(addprefix build/src/,$(modules))
test_dir   = build/test/
cflags     = -std=c++17 $(opt) -I src/
main       = src/main.cpp

source     = $(foreach sdir,$(source_dir),$(filter-out $(main), $(wildcard $(sdir)/*.cpp)))
objects    = $(patsubst src/%.cpp,build/src/%.o,$(source))

vpath %.cpp $(source_dir)

define make-goal
$1/%.o: %.cpp
	$(cc) $(cflags) -c $$< -o $$@
endef

.PHONY : all

all : checkdirs $(objects) $(exec)

checkdirs : $(build_dir)

test : $(test_dir)
	@ mkdir -p $(test_dir)
	@ cd $(test_dir) && $(MAKE)

$(build_dir):
	@ mkdir -p $@

$(exec) : $(main)
	@ rm -f $(exec)
	@ $(cc) $(cflags) $(objects) $< -o $@

clean :
	@ rm -rf $(build_dir)
	@ rm -rf $(exec)*
	@ cd  $(test_dir) && $(MAKE) clean

$(foreach bdir,$(build_dir),$(eval $(call make-goal,$(bdir))))
