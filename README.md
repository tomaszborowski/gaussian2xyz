# gaussian2xyz
Python3 script extracting geometry from Gaussian output (log) file

The script expects 2 or 3 arguments: \
    #1 log-file-name or, for irc_f option, the name of a file specifying filenames and directions of IRC calculations (optionally also the number of IRC points to be extracted)\
    #2 type of extraction - one from: scan, irc, irc_f, all, last, nr, lasts, nrs\
    #3 (for IRC) file name with SP/FREQ calculations for the TS from which IRC calculations started\
        or for NR 1-based number of the structure to be extracted

Type of extractions:\
scan - for scan calculations - only final geometries for all steps of the scan are extracted\
irc - for IRC calculations - only final geometries for all IRC steps are extracted\
irc_f - for IRC calculations when the IRC path is covered by several log files - only final geometries for all IRC steps are extracted\
all - for any calculations - all geometries are extracted\
last - for any calculations - the last geometry is extracted\
nr - for any calculations - the nr'th geometry is extracted\
lasts - for any calculations - the last geometry of the QM part of the system with Mulliken spin populations is extracted\
nrs - for any calculations - the nr'th geometry of the QM part of the system with Mulliken spin populations is extracted

Example usage:\
gaussian2xyz.py oh_h2o.scan.log scan > oh_h2o.scan.xyz\
gaussian2xyz.py oh_h2o.irc.log irc oh_h2o.ts_fq.log > test.irc.xyz\
gaussian2xyz.py oh_h2o.irc.info irc_f > test.irc_f.xyz


Example content of the input file (oh_h2o.irc.info) to be used with the irc_f option:\
oh_h2o.irc_1.log       reverse	N1\
oh_h2o.irc_1b.log      reverse  N1b\
oh_h2o.ts_fq.log       ts\
oh_h2o.irc_2.log       forward	N2\
oh_h2o.irc_2b.log      forward	N2b

The order of the files specified in the input file is important! 
For a given direction ("reverse" or "forward") they should be given in the sequence 
of calculations starting from the TS structure. 
N<x> is optional, if not given, the script assumes the 2nd log for a given direction 
resumes at the last point from the 1st log; 3rd from the last geom from the 2nd, etc. 

IF N<x> is present - the script reads first N1 points from the 1st log, assumes the 2nd log 
started from the last point read from the 1st log (the N1st one) and reads the first N1b points
from 2nd log, etc.


