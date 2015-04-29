
import StandardParser
import os

'''
This parser accepts annotations of this form:

chr6 - 53471513 53517402
chr17 + 43544405 43553869
chr17 + 43544405 43553869
chr12 - 8985681 8990267
chr8 - 17201742 17315172
'''

class AnnotationParser(StandardParser.AnnotationParser):
  def parseLine(self, l):
    chr, strand, start, stop = l.rstrip().split()
    return (chr, start, stop, strand)

''' this parser accepts snp signals of this form:
chr1	45027
chr1	45162
chr1	48677
chr1	51325
chr1	52066

It assumes that the input is sorted by chr and position.  The dependency on that is a bit
subtle.
'''

class SignalParser(object):
  def __init__(self, flist):
    self.flist=flist

  def parse(self):
    signals = {}

    for f in self.flist:
      for l in open(f):
          c, pos = l.rstrip().split()
          signals.setdefault(c, {})[int(pos)]=1
          signals.setdefault(c, {})[int(pos)+1]=0
    return signals

