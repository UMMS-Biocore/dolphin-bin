
''' This class parses standard bed files

The file contains a tab-delimited set of regions of interest on the genome.  Each line describes
one annotation: chr<TAB>start<TAB>stop<TAB>strand.
chr22	20337315	20337337	+
chr22	20337643	20337667	+
chr22	29081918	29081941	+
chr22	29886045	29886078	-
chr22	30072214	30072236	+
chr22	34376215	34376238	+

'''

import sys

class AnnotationParser(object):
  def __init__(self, f):
    self.f=f

  def parseLine(self, l):
    chr, start, stop, strand = l.rstrip().split()
    return (chr, start, stop, strand)

  def parse(self):
    annotations = {}
    for l in open(self.f):
      chr, start, stop, strand = self.parseLine(l)
      chr = annotations.setdefault(chr, {})
      if strand=='+':
        chr[int(start)]=(start, stop, "+")
      else:
        chr[int(stop)]=(start, stop, "-")
    return annotations

'''
This class parses signal files containing a tab-delimited set of signals.  Each line describes the beginning of a particular
stair-stepped signal level:  chr<TAB>pos<TAB>level

chr22	20337221	7
chr22	20337261	8
chr22	20337285	7
chr22	20337288	6
chr22	20337311	5
chr22	20337325	6
'''

class SignalParser(object):
  def __init__(self, flist):
    self.flist=flist

  def parseLine(self, l):
    chr, pos, height = l.rstrip().split()
    return (chr, pos, height)

  def parseConservative(self):
    signals = {}
    for f in self.flist:
      firstRec=True
      for l in open(f):
        chr, pos, height = self.parseLine(l)
        if firstRec:
          if signals.has_key(chrm):
            raise Exception("Found %s in multiple signal files.  This is not allowed." % chrm)
          firstRec=False
          signals.setdefault(chr, {})[int(pos)]=int(round(float(height)))
    return signals

# We start out assuming that signals are sorted, and that we can remove redundant entries.  If we discover this isn't the case
# we'll give up and return everything.  This handles the most common case well.
# We'll also check that signals for a given chrom appear in only one file.  This is the simplest way to prevent overlap, and
# again handles the most common case efficiently.
  def parse(self):
    signals = {}
    lastchrm=None; lastheight=None; lastpos=None
    for f in self.flist:
      firstRec=True
      for l in open(f):
        chrm, pos, height = self.parseLine(l)
        if firstRec:
          if signals.has_key(chrm):
            raise Exception("Found %s in multiple signal files.  This is not allowed." % chrm)
          firstRec=False

        if lastchrm==chrm and int(pos)<=lastpos: # bail out
          print >>sys.stderr, "Found position out of sort order.  Falling back to conservative parser"
          return self.parseConservative()
        if lastchrm==chrm and lastheight==height: # no new info here
          continue
        lastchrm=chrm; lastheight=height; lastpos=int(pos)
        signals.setdefault(chrm, {})[int(pos)]=int(round(float(height)))
    return signals

