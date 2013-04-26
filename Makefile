TOP:=.
.PHONY: polyselect sqrt tags etags
# This makefile is a placeholder. Please have a look to $(TOP)/scripts/call_cmake.sh,
# and (possibly) edit a file $(TOP)/local.sh to tweak your build preferences.
all polyselect sqrt: ; +@$(TOP)/scripts/call_cmake.sh $@
tags: ; grep -h '^[^#].*\.[ch]' files.dist files.nodist | xargs ctags
etags: ; grep -h '^[^#].*\.[ch]' files.dist files.nodist | xargs etags
show variables install tidy cmake dist: ; +@$(TOP)/scripts/call_cmake.sh $@
%: ; +@$(TOP)/scripts/call_cmake.sh $@
