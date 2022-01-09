# gaussian2xyz
Python3 script extracting geometry from Gaussian output file

The script expects 2 or 3 arguments: 
    #1 log-file-name or, for irc_f, name of a file specifying filenames and directions of irc calc.
    #2 type of extraction - one from amoung: scan, irc, irc_f, all, last, 
    #3 (for IRC) file name with SP/FREQ calculations for the TS from which IRC calculations started
    
Example usage:
gaussian2xyz.py oh_h2o.scan.log scan > oh_h2o.scan.xyz

gaussian2xyz.py oh_h2o.irc.log irc oh_h2o.ts_fq.log > test.irc.xyz

gaussian2xyz.py oh_h2o.irc.info irc_f > test.irc_f.xyz


Example content of the input file (oh_h2o.irc.info) to be used with the irc_f option:

oh_h2o.irc_1.log       reverse
oh_h2o.irc_1b.log      reverse
oh_h2o.ts_fq.log       ts
oh_h2o.irc_2.log       forward
oh_h2o.irc_2b.log      forward

The order of the files specified in the input file is important! 
For a given direction ("reverse" or "forward") they should be given in the sequance 
of calculations starting from the TS structure. 
The script assumes the 2nd log for a given direction resumes at the last point from the 1st log; 
3rd from the last geom from 2nd, etc. 


