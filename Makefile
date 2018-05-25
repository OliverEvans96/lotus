# Set "all" as default target
all:

# Actual make file is in include directory
%:
	cd include && $(MAKE) $@

.PHONY: test

test:
	cd include && $(MAKE) $@
