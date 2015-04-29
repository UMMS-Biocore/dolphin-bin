
'''This file contains a parser that accepts a tab delimited set of regions as shown below.  It returns 
a region that begins in the center of the region and stretches for 1000bp.

accept:

chr11	62364364	62366285	2075	306	6.8012742339671	1769	0
chr1	148123142	148127230	2913	319	9.15470364075437	2594	0
chr17	8016743	8017870	1829	113	16.2547048846013	1716	0
chr6	26264348	26267282	5688	334	17.061221202399	5354	0

convert to:

chr11   6236436474  +
'''

import StandardParser

class AnnotationParser(StandardParser.AnnotationParser):

  def parseLine(self, l):
    c, start, stop = l.rstrip().split('\t')[0:3]
    mid = (int(start) + int(stop))/2
    start = mid
    stop = mid+1000
    return (c, start, stop, '+')


