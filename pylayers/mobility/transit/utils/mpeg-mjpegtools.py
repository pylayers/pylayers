"""

You need the netpbm-progs and mjpegtools package to create an MPEG.

TkWorld outputs PostScript files for each frame.  This script takes a
sequence of PostScript files (200 or so), converts them to PPM using
pstopnm.

Use 'pstopnm FILE.ps | pnmcrop -verbose' to get a rough idea of what
to pass for pnmcut parameters.

mkdir images-ps images-ppm
PYTHONPATH=. python sims/entry.py
python utils/mpeg-mjpegtools.py 0 310 132 1280 40 720 entry.mpeg
cat images-ppm/* \
  | ppmtoy4m | y4mscaler -O chromass=420mpeg2 \
  | mpeg2enc -M 4 -o entry.mpeg

"""

import os, sys, re

first_frame, last_frame, left, width, top, height, mpeg_name = sys.argv[1:]

os.system('rm images-ppm/*.ppm ' + mpeg_name)

for ii in range(int(first_frame), int(last_frame)):
    in_file_name = "images-ps/%05d.ps" % ii
    out_file_name = "images-ppm/%05d.ppm" % ii
    os.system('pstopnm -dpi=125 -landscape -stdout %s | pnmflip -lr | pnmcut -left %s -width %s -top %s -height %s >%s'
              % (in_file_name, left, width, top, height, out_file_name))
