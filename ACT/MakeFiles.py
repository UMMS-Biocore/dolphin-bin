
import sys
import Hit2Annotation, SNP2Signal

annfp = open('hits.bed', 'w')
sigfp = open('sig.sgr', 'w')

a=Hit2Annotation.parseAnnotations(sys.argv[1])
for c, d in a.iteritems():
    ks = sorted(d.keys())
    for k in ks:
        e=d[k]
        print >> annfp, "%s\t%d\t%d\t%s" % (c, e[0], e[1], e[2])

s=SNP2Signal.parseSignals((sys.argv[2],))
for c, d in s.iteritems():
    ks = sorted(d.keys())
    for k in ks:
        e=d[k]
        print >> sigfp, "%s\t%d\t%d" % (c, k, e)
