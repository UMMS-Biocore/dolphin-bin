
''' This parser returns a set of regions 1000bp in length, randomly placed on the genome.  The number of regions is 
specified by an environment variable N
'''

import WeightedSample
import random, os

d={
        'chr1':[(1,247249690)],
        'chr2':[(1,242951120)],
        'chr3':[(1,199501798)],
        'chr4':[(1,191273034)],
        'chr5':[(1,18085783)],
        'chr6':[(1,170899963)],
        'chr7':[(1,158821395)],
        'chr8':[(1,146274797)],
        'chr9':[(1,140273223)],
        'chr10':[(1,135374708)],
        'chr11':[(1,134452355)],
        'chr12':[(1,132349505)],
        'chr13':[(1,114142951)],
        'chr14':[(1,106368556)],
        'chr15':[(1,100338886)],
        'chr16':[(1,88827225)],
        'chr17':[(1,78774713)],
        'chr18':[(1,76117124)],
        'chr19':[(1,63811622)],
        'chr20':[(1,62435935)],
        'chr21':[(1,46944294)],
        'chr22':[(1,49691403)],
        #'chrM':[(1,16542)],
        'chrX':[(1,154913725)],
        #'chrY':[(1,57772925)],
        }

class AnnotationParser(object):
  def parse(self):
    n=int(os.getenv('N'))
    annotations = {}

    weights = {}
    total=0
    for k,v in d.iteritems():
      total += v[0][1]
      for k,v in d.iteritems():
        weights[k] = float(v[0][1])/total

    sampler = WeightedSample.Sampler(weights)

    for i in xrange(n):
      c = sampler.sample()
      start = random.randint(1,d[c][0][1])
      stop = start+1000
      c = annotations.setdefault(c, {})
      c[int(start)]=(int(start), int(stop), '+')
    return annotations

