"""

To create an MPEG, you need mpeg2encode from

  ftp://mm-ftp.cs.berkeley.edu/pub/multimedia/mpeg2/software/mpeg2vidcodec_v12.tar.gz

The PPMtoMPEG.par file included in this directory came from

  http://tns-www.lcs.mit.edu/manuals/mpeg2/PPMtoMPEG.par

The NetPBM tools, usually part of most Linux distributions, is used to
convert PostScript to PPM files.

TkWorld outputs PostScript files for each frame.  This script takes a
sequence of PostScript files (200 or so), converts them to PPM using
pstopnm, and runs mpeg2encode.

Use 'pstopnm FILE.ps | pnmcrop -verbose' to get a rought idea of what
to pass for pnmcut parameters.

"""

import os, sys, re

first_frame, last_frame, left, width, top, height, mpeg_name = sys.argv[1:]

os.system('rm images-ppm/*.ppm ' + mpeg_name)

par = open('utils/PPMtoMPEG.par.in', 'r').read()
par = re.sub('@WIDTH@', width, par)
par = re.sub('@HEIGHT', height, par)
open('utils/PPMtoMPEG.par', 'w').write(par)

jj = int(first_frame)
while jj < int(last_frame):
    for ii in range(0, 200):
        in_file_name = "images-ps/%05d.ps" % (jj + ii)
        out_file_name = "images-ppm/%05d.ppm" % ii
        os.system('pstopnm -stdout %s | pnmcut -left %s -width %s -top %s -height %s >%s'
                  % (in_file_name, left, width, top, height, out_file_name))
    os.system('mpeg2encode utils/PPMtoMPEG.par tmp.mpg >/dev/null')
    os.system('cat tmp.mpg >>%s; rm tmp.mpg' % mpeg_name)
    jj += 200
