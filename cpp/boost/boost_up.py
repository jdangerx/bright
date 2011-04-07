#! /usr/bin/env python

import os
import re


FILES = [
    "boost/math/special_functions/bessel.hpp",
    ]


include_pattern = """#include\s+[<'"](.+?)["'>]"""


def copy_path(path_list):
    """Reads in a path, replaces 'include <>' with 
    'include ""', writes the path to the same directory
    structure and and returns a set of the includes
    found in this file.
    """
    extra_paths = set()

    # Set up directories
    dirs = path_list[2:-1]
    if 0 < len(dirs):
        dirs_path = [os.path.join(*dirs[:i]) for i in range(1, len(dirs)+1)]
        [os.mkdir(p) for p in dirs_path if not os.path.exists(p)]

    # Make new file
    fnew = ""
    fpath = os.path.join(*path_list)
    with open(fpath, 'r') as f:
        for line in f:
            m = re.search(include_pattern, line)

            if (m is not None) and (m.group(1).split('/')[0] == 'boost'):
                extra_path = m.group(1)
                extra_paths.add(extra_path)

                line = line.replace('<', '"')
                line = line.replace('>', '"')

            fnew += line

    # Write the new file
    fnewpath = path_list[2:]
    fnewpath = os.path.join(*fnewpath)
    with open(fnewpath, 'w') as f:
        f.write(fnew)

    return extra_paths


def copy_boost(bdir):
    """Takes the root boost directory and copies all files 
    in the FILES list, and all dependent files, over to this 
    directory.
    """
    # Split 'boost' off the end of the path
    bdir = os.path.split(bdir)[0]

    queue = set(FILES)
    copied = set()

    # Keep going until all files have been copied over.
    while queue != copied:
        # Loop through everything in the queue
        for f in queue:
            # Skip previously done files
            if f in copied:
                continue

            print "Copying file {0}".format(f)

            # Get os version of the path
            fs = f.split('/')
            path_list = [bdir] + fs

            extra_paths = copy_path(path_list)
            queue = (queue | extra_paths)

            copied.add(f)


if __name__ == "__main__":
    path_to_boost = "/home/scopatz/Downloads/boost_1_46_1/boost"
    copy_boost(path_to_boost)
