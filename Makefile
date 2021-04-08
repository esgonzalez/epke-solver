exec           = epke-run
test_exec      = epke-test
cc             = g++
opt            = -g
modules        = pugi parareal epke utility
source_dir     = $(addprefix src/,$(modules))
test_dir       = $(addprefix test/,$(modules))
build_src_dir  = $(addprefix build/src/,$(modules))
build_test_dir = $(addprefix build/test/,$(modules))
cflags         = -std=c++17 $(opt) -I src/
main           = src/main.cpp
test_main      = test/test.cpp

source         = $(foreach sdir,$(source_dir),$(filter-out $(main), $(wildcard $(sdir)/*.cpp)))
test_source    = $(foreach tdir,$(test_dir),$(filter-out $(test_main), $(wildcard $(tdir)/*.cpp)))
objects        = $(patsubst src/%.cpp,build/src/%.o,$(source))
test_objects   = $(patsubst test/%.cpp,build/test/%.o,$(test_source))

vpath %.cpp $(source_dir)
vpath %.cpp $(test_dir)

define make-goal
$1/%.o: %.cpp
	$(cc) $(cflags) -c $$< -o $$@
endef

.PHONY : all

all : checkdirs $(objects) $(exec)

checkdirs : $(build_src_dir) $(build_test_dir)

test : checkdirs $(test_objects) $(test_exec)

$(build_src_dir):
	@ mkdir -p $@

$(build_test_dir):
	@ mkdir -p $@

$(exec) : $(main)
	@ rm -f $(exec)
	@ $(cc) $(cflags) $(objects) $< -o $@

$(test_exec) : $(test_main)
	@ rm -f $(test_exec)
	@ $(cc) $(cflags) -I test/ $(test_objects) $< -o $@
	@ ./$(test_exec)

clean :
	@ rm -rf $(build_src_dir)
	@ rm -rf $(build_test_dir)
	@ rm -rf $(exec)*
	@ rm -rf $(test_exec)*

$(foreach bdir,$(build_src_dir),$(eval $(call make-goal,$(bdir))))
$(foreach bdir,$(build_test_dir),$(eval $(call make-goal,$(bdir))))
