# coding: utf-8
# Create api/*.rst files.
# Run from repo root dir.

import os
files = os.listdir('include')
names = [f[:-2] for f in files if f[-2:] == '.h' and f[0] != '.' and f != 'yaml.h']

template = """.. _`{name}.rst`:

{name}.h
{line}
.. doxygenfile:: {name}.h
"""

for name in names:
    # +2 for .h
    line = '=' * (len(name) + 2)
    txt = template.format(name=name, line=line)
    path = 'docs/source/api/{name}.rst'.format(name=name)
    with open(path, 'w') as fh:
       fh.write(txt)
