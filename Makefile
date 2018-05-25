# Set "all" as default target
all:

# Actual make file is in include directory
%:
	cd include && $(MAKE) $@

# Files in current directory should be treated as phony
# since there are no top-level targets.
phony=$(wildcard ./*)

.PHONY: $(phony)

$(phony):
	cd include && $(MAKE) $@
