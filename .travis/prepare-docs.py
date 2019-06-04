#!/usr/bin/env python3

import os
import fileinput

SNIPPETS_ROOT = 'snippets'

lines = ['- [Home](/README.md#getting-started)\n']
paths = ['/']

for root, dirs, files in os.walk(SNIPPETS_ROOT, topdown=True):
  for file in files:
    if file.endswith('README.md'):
      level = root.count(os.sep) - SNIPPETS_ROOT.count(os.sep) - 1
      path = os.path.join(root, file)
      with open(path, 'rt') as md:
        title = md.readline()
        title = title.replace('#', '').strip()
      offset = '  ' * level
      line = f'{offset}- [{title}]({path})\n'
      lines.append(line)
      paths.append(path)

with open('_sidebar.md', 'wt') as f:
  f.writelines(lines)


for line in fileinput.FileInput('index.html', inplace=True):
    line = line.replace('paths: []', f'paths: {paths}')
    print(line, end='')
