#!/usr/bin/env python

# Supply the mapping file (arg 1) and the suffix without the extension of the files you want to rename in bulk (arg 2)
# column 1 = TwinID, column 3 = BAM ID.

import os,sys

# A dict with keys being the old filenames and values being the new filenames
mapping = {}

# Read through the mapping file line-by-line and populate 'mapping'
with open(sys.argv[1]) as mapping_file:
    header=next(mapping_file)
    for line in (line.strip().split() for line in mapping_file):
        # Note: this fails if your filenames have whitespace
        if line[3] == 1:
        	pass
        else:
        	mapping[line[2]] = line[0]

#e.g. "_sorted"
suffix = sys.argv[2]

# List the files in the current directory
for filename in os.listdir('.'):
    root, extension = os.path.splitext(filename)
    if not root.endswith(suffix):
        # File doesn't end with this suffix; ignore it
        continue
    # Strip off the number of characters that make up suffix
    stripped_root = root[:-len(suffix)]
    if stripped_root in mapping:
        os.rename(filename, ''.join(mapping[stripped_root] + suffix + extension))
