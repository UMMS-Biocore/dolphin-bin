
import os
import get1000GSNPAnnotations

class SignalParser(object):
  def __init__(self, flist):
    self.flist=flist

  def parse(self):
    signals = {}

    for f in self.flist:
      h=get1000GSNPAnnotations.Handler(f)
      for c, pos, snp in h.getAllAnnotationsGenerator():
        c='chr'+c
        hetsnp = os.getenv("HETSNP")
        if not hetsnp or hetsnp=='0' or snp[7]: # hetSNP
            signals.setdefault(c, {})[int(pos)]=1
            signals.setdefault(c, {})[int(pos+1)]=0
    return signals

