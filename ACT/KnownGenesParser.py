
'''
This file contains a parser that accepts the UCSC knowngenes annotation file, and returns a set of regions that 
begin at each TSS site.

It demonstates how to subclass the standard parser when the only difference is a new file format

accept:

uc001aaa.2	chr1	+	1115	4121	1115	1115	3	1115,2475,3083,	2090,2584,4121,		uc001aaa.2
uc009vip.1	chr1	+	1115	4272	1115	1115	2	1115,2475,	2090,4272,		uc009vip.1
uc001aab.2	chr1	-	4268	14764	4268	4268	10	4268,4832,5658,6469,6716,7095,7468,7777,8130,14600,	4692,4901,5810,6628,6918,7231,7605,7924,8242,14764,		uc001aab.2
uc009viq.1	chr1	-	4268	19221	4268	4268	7	4268,5658,6469,6720,7468,14600,19183,	4692,5810,6628,6918,7924,14754,19221,

convert to:

chr1	1115    2115   +	
'''

import StandardParser

class AnnotationParser(StandardParser.AnnotationParser):

  # only difference from standard parser is the order of the fields
  def parseLine(self, l):
    elems=l.rstrip().split('\t')
    chr, start, stop, strand = [elems[idx] for idx in [1,3,4,2]]
    return (chr, start, stop, strand)


