'''
Written by Robert Bjornson bjornson@yale.edu, based on previous versions written by Joel Rozowsky
and Justin Jee, also at Yale.

This code takes an annotation file and one or more signal files, and computes a histogram of the 
signal intensity around annotations.  

Inputs:

By default the code expects the annotations and signals to be in the following formats.  This is the form
parsed by StandardParser.py.  However, alternative parsers can be specified for other file types, as 
discussed below.

Annotations:  This file contains a tab-delimited set of regions of interest on the genome.  Each line describes
one annotation: chr<TAB>start<TAB>stop<TAB>strand.  Here are some example lines:
chr22	20337315	20337337	+
chr22	20337643	20337667	+
chr22	29081918	29081941	+
chr22	29886045	29886078	-
chr22	30072214	30072236	+
chr22	34376215	34376238	+

Signals:  This file contains a tab-delimited set of signals.  Each line describes the beginning of a particular
stair-stepped signal level:  chr<TAB>pos<TAB>level

chr22	20337221	7
chr22	20337261	8
chr22	20337285	7
chr22	20337288	6
chr22	20337311	5
chr22	20337325	6

The descriptors used for chromosomes can be arbitrary strings, e.g. 22, chr22, ch22, but they must be consistent 
within and between files.

Custom Parsers:

Custom parsers can be defined to allow for different specification of Annotations and Signals.  The --annotationparser
and --signalparser flags specify the module name.  The specified module must define a class called
AnnotationParser or SignalParser, respectively.  These classes must define a function parse().
See KnownGenesParser.py or SNP2Signal.py for examples.

AnnotationParser.parse() should accept a filename and return a dictionary of this form:
{'chr1':{<pos>: (<start>, <stop>, +/-), <pos>: (<start>, <stop>, +/-), ...},
{'chr2': ...
...
}

SignalParser.parse() accepts a list of filenames and returns a single dictionary of this form:
{'chr1':{<pos>: <signal>, <pos>: <signal>, ...},
{'chr2': ...
...
}

Please note that StandardParser.parse() can be run on multiple files.  However, the files must be
disjoint by chromosome.  For example, signals for chr1 can only appear in one file.  This feature was
added to allow splitting very large signal files into multiple smaller files by chromosome, to
avoid memory issues.  If you want to do something more clever, you can easily define your own signal parser.

In the sample parsers in the distribution, two techniques are used to create custom parsers.  
1) Subclass StandardParser and redefine parseLine().  This is a very easy way to create a custom parser if the only
change is that the input lines are in a different format (columns in a different order, for example).  See
KnowGenesParser.py for an example.
2) Create a standalone parser class.  This is useful when you're doing something more fundamentally different.
See SNP2Signal.py for an example.

Output:
The program reports the location and signal level for each of the bins.

Methods:

Basically, the program lays a set of bins around each annotation, calculates the mean or median signal level in 
each bin.  The final result is the accumulated data for all annotations.  

The bins can be laid out in two ways:

1) Fixed size bins.  2*nbins bins, each of length radius/nbins.  This is the default, but also corresponds to --point

2) A fixed number (mbins) of bins of variable size over the annotation proper, with nbins on each size of the annotation.
Use --region to choose this layout.

See act.gersteinlab.org

'''

# Change this if you want to use different parsers for the annoation and signal files
import StandardParser as Parser

import mystats
import sys, bisect, getopt, math

def median(l):
  sortL=sorted(l)
  elems=len(sortL)
  if elems%2==1:
    return sortL[elems/2]
  else:
    return (sortL[elems/2-1]+sortL[elems/2])/2.0
    
class Finished(Exception):
  pass

# This code comes from http://www.johndcook.com/standard_deviation.html
# the point of all this complexity is to allow incremental computation of mean and std in
# a numerically stable way.
class accumMean(object):
  def __init__(self):
    self.m_n=self.m_oldM=self.m_newM=0
  def fields(self):
    return ["mean", "stdev"]
  def add(self, x, binWidth):
    x=float(x)/binWidth
    self.m_n+=1
    if self.m_n==1:
      self.m_oldM = self.m_newM = x
      self.m_oldS = 0.0
    else:
      self.m_newM = self.m_oldM + (x - self.m_oldM)/self.m_n
      self.m_newS = self.m_oldS + (x - self.m_oldM)*(x-self.m_newM)
      # set up for next iteration
      self.m_oldM = self.m_newM
      self.m_oldS = self.m_newS
  def finishPosition(self):
    pass # nothing to do for mean
  def finalize(self, count):
    mean=float(self.m_newM)*self.m_n/count
    if self.m_n>1:
        stdev = math.sqrt(self.m_newS/(self.m_n - 1))
    else:
        stdev = 0.0
    self.val={"mean":mean, "stdev":stdev}

  def value(self):
    return self.val

class accumMedian(object):
  def __init__(self):
    self.valsThisPos=[]
    self.valsAllPos=[]
  def fields(self):
    return ["mean", "stdev", "min", "Q1", "median", "Q3", "max"]
  def add(self, signal, binWidth):
    self.valsThisPos.append(signal)
  def finishPosition(self):
    if self.valsThisPos:
      self.valsAllPos.append(median(self.valsThisPos))
      self.valsThisPos=[]
  def finalize(self, count):
    if not self.valsAllPos:
      vs=[0]
    else:
      vs=self.valsAllPos

    self.val={"mean":mystats.mean(vs),
              "median":mystats.percentile(vs, 50),
              "stdev":mystats.std(vs),
              "min":min(vs),
              "max":max(vs),
              "Q1":mystats.percentile(vs, 25),
              "Q3":mystats.percentile(vs, 75)}

  def value(self):
    return self.val

def accumFactory(method):
  if method=='mean':
    return accumMean()
  else:
    return accumMedian()

def makeRelPos(nbins, mbins, radius, intragenicp, genelength):
  if intragenicp:
    nw=int(float(radius)/nbins + 0.5)
    mw=int(float(genelength)/mbins + 0.5)
    t=mw*mbins
    return ([((i-nbins)*nw,(i-nbins+1)*nw-1)  for i in xrange(nbins)] +
    [(i*mw,(i+1)*mw-1)  for i in xrange(mbins)]+
    [(t+i*nw,t+(i+1)*nw-1)  for i in xrange(nbins)])
  else: 
    w=max(1,int(float(radius)/nbins + 0.5))
    return [((i-nbins)*w,(i-nbins+1)*w-1)  for i in xrange(nbins*2)]    

def contribution(signalStart, signalEnd, signalLevel, bucketStart, bucketEnd, method):
  if bucketEnd < signalStart: return None
  if signalEnd < bucketStart: raise Finished
  # if we got here, they overlap
  if method=='mean': 
    return (min(signalEnd, bucketEnd) - max(signalStart, bucketStart)+1)*signalLevel
  else:
    return signalLevel

def getFlankingSignals(signals, start, end):
  assert(start<=end)
  startidx = max(bisect.bisect_left(signals, start)-1, 0)
  # This actually gives me the first location not in range, but the array slice below
  endidx=bisect.bisect_right(signals, end)
  assert(startidx<=endidx)
  return signals[startidx:endidx]

def getFlankingPairs(l):
  numpairs=len(l)-1
  pairs = [(l[i],l[i+1]-1) for i in range(numpairs)]
  pairs.append((l[-1],sys.maxint))
  return pairs

USAGE="%s [--nbins=#] [--mbins=#] [--radius=#] [--median] [--mean] [--region] [--point] [--annotationparser=<name>] [--signalparser=<name>] [--mingenelen=#] [--output=<file>] annotationFile signalFile(s)"

if __name__ == '__main__':
  if len(sys.argv) < 3:
    print USAGE % sys.argv[0]
    sys.exit(1)

  opts, args = getopt.getopt(sys.argv[1:], "", ["nbins=", "mbins=", "radius=", "median", "mean", "region", "point", "signalparser=", "annotationparser=", "mingenelen=", "output="])

  # defaults
  nbins=34
  mbins=0
  regionp=False
  radius=nbins*90
  method="mean"
  mingenelen=0
  signalParserModule=annotationParserModule="StandardParser"
  outfp=sys.stdout

  for o,a in opts:
    if o == '--nbins':
      nbins = int(a)
    elif o == '--mbins':
      mbins = int(a)
    elif o == '--radius':
      radius = int(a)
    elif o == '--median':
      method='median'
    elif o == '--mean':
      method = 'mean'
    elif o == '--region':
      regionp = True
    elif o == '--point':
      regionp = False
    elif o == '--mingenelen':
      mingenelen = int(a)
    elif o == '--annotationparser':
      annotationParserModule = a
    elif o == '--signalparser':
      signalParserModule = a
    elif o == '--output':
      outfp = open(a, "w")
      
  error=False
  # now a little sanity checking
  if nbins < 10:
    print >>sys.stderr, "nbins must be >= 10"
    error=True
  if regionp and (mbins < 10):
    print >>sys.stderr, "mbins must be >= 10 for region runs"
    error=True
  if radius < 100:
    print >>sys.stderr, "radius must be >= 100"
    error=True
  if not regionp and mbins != 0:
    print >>sys.stderr, "mbins must be 0 for point runs"
    error=True

  if error:
    sys.exit(1)

  exec "import %s as annotationParserModule" % annotationParserModule
  exec "import %s as signalParserModule" % signalParserModule

  annotationFile=args[0]
  signalFiles=args[1:]

  annotations = annotationParserModule.AnnotationParser(annotationFile).parse()
  signals = signalParserModule.SignalParser(signalFiles).parse()

  print >>sys.stderr, "Done reading input files"
  print >>outfp, "# " + " ".join(sys.argv)
  #print >>outfp, "# annotations:%(annotationFile)s signals:%(signalFiles)s nbins:%(nbins)d mbins:%(mbins)d radius:%(radius)d region:%(regionp)s mingenelen:%(mingenelen)d method:%(method)s point:%(regionp)s signal Parser:%(signalParserModule)s annotation Parser:%(annotationParserModule)s" % locals()

  numbins = 2*nbins+mbins
  accums = [accumFactory(method) for i in xrange(numbins)]
  annotationCount = 0

  for chr, posTable in annotations.iteritems(): 
    annotationCount+=len(posTable)
    # process one chromosome's worth
    sortedpos = sorted(posTable.keys())
    try:
      signalTable = signals[chr]
      sortedsigs = sorted(signalTable.keys())
    except KeyError:
      sortedsigs=[]

    for pos in sortedpos:
      # process one annotation position on this chromosome
      genelength=int(posTable[pos][1])-int(posTable[pos][0])
      if genelength < mingenelen:
        continue
      bins=makeRelPos(nbins, mbins, radius, regionp, genelength)

      # convert relative bin annotations to absolute for this position
      absbins=[(s+pos, e+pos) for s,e in bins]
      
      # get flanking signal indices, i.e. all signals within my bins, plus one outside
      flankingSignals = getFlankingSignals(sortedsigs, absbins[0][0], absbins[-1][1]) # array of signal positions, sorted
      if not flankingSignals:
        continue
      flankingPairs=getFlankingPairs(flankingSignals)

      startBin=lastbidx=0
      for signalStart, signalEnd in flankingPairs: # signal end is one less than the next signal position
        for bidx in xrange(startBin, numbins):
          binStart, binEnd=absbins[bidx]
          try:
            # returns mean:overlap*signal, median:signal
            contrib = contribution(signalStart, signalEnd, signalTable[signalStart], binStart, binEnd, method)
            if contrib==None: continue
            binWidth=binEnd-binStart+1
            if binWidth > 0:
              if posTable[pos][2]=='+':
                accums[bidx].add(contrib,binWidth)
              else:
                accums[-bidx-1].add(contrib, binWidth)
          except Finished: # the current bin is past the current signal
            startBin = lastbidx
            break
          lastbidx=bidx

      for a in accums: 
        a.finishPosition()

  for a in accums: 
    a.finalize(annotationCount)

  print >>outfp, "# annotationCount: %d" % annotationCount

  # construct and print a header
  fields = accums[0].fields()

  if regionp:
    print >>outfp, "\t".join(["Bin"] + fields)
  else:
    print >>outfp, "\t".join(["Bin", "Center"] + fields)

  for i, a in enumerate(accums):
    val=a.value()
    if regionp:
      print >>outfp, "%d\t" % (i-nbins),
    else:
      center=(bins[i][0]+bins[i][1])/2.0
      print >>outfp, "%d\t%d\t" % (i-nbins, center),
    print >>outfp, '\t'.join([str(val[fn]) for fn in fields])

