#!/usr/bin/env python

import glob
import os

ROOT = 'http://nimble-dev.github.io/nimble-demos'
TOC_LINE = '## Table of contents'

lines = []
with open('README.md') as f:
    for line in f:
        line = line.strip()
        if line == TOC_LINE:
            break
        lines.append(line)

lines += ['', TOC_LINE]

for path in glob.glob('*/*.html'):
    dirname, basename = os.path.split(path)
    dirname = dirname.replace('_', ' ').capitalize()
    basename = basename[:-5].replace('_', ' ').capitalize()
    lines += ['', '[{}: {}]({}/{})'.format(dirname, basename, ROOT, path)]

with open('README.md', 'w') as f:
    f.write('\n'.join(lines))
