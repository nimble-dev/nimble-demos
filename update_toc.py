#!/usr/bin/env python

import glob
import os

lines = [
    '<h1>Nimble Demos</h1>',
    '',
    'Demos for the <a href="http://r-nimble.org">NIMBLE</a> project.',
    '',
    '<h2>Table of Contents</h2>'
    ,
    '',
]

for path in glob.glob('*/*.html'):
    dirname, basename = os.path.split(path)
    dirname = dirname.replace('_', ' ').capitalize()
    basename = basename[:-5].replace('_', ' ').capitalize()
    lines += ['', '<a href="{}">{}: {}</a>'.format(path, dirname, basename)]

with open('index.html', 'w') as f:
    f.write('\n'.join(lines))
