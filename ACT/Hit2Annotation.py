
'''This file contains a parser that accepts a tab delimited set of regions as shown below.  It returns 
a region that begins in the center of the region and stretches for 1000bp.

accept:

chr7	244017	244020
chr7	244022	244023
chr7	244070	244071

convert to:

chr7    244018  +
'''

import StandardParser

class AnnotationParser(StandardParser.AnnotationParser):

  def parseLine(self, l):
    c, start, stop = l.rstrip().split('\t')[0:3]
    mid = (int(start) + int(stop))/2
    start = mid
    stop = mid+1000
    return (c, start, stop, '+')

