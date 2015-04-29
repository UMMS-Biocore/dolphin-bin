
import bisect

class annotationHandler(object):
    def __init__(self, f=None, hasHeader=False):
        self.dPos={}
        self.dAnnotations={}
        if f:
            self.load(f, hasHeader)

    def getClasses(self):
        return self.dPos.keys()

    def getClassPositions(self, c):
        return self.dPos[c]

    def getAllAnnotations(self):
        return [(c, pos, self.dAnnotations[c][pos]) for c, poslist in self.dPos.iteritems() for pos in poslist]

    def getAllAnnotationsGenerator(self):
        for c, poslist in self.dPos.iteritems():
            for pos in poslist:
                yield (c, pos, self.dAnnotations[c][pos])

    def getAnnotations(self, c, frm=0, to=99999999999):
        frm=int(frm)
        to=int(to)
        assert(frm<=to)
        if not self.dPos.has_key(c):
            return ([],[])
        start = bisect.bisect_left(self.dPos[c], frm)
        # This actually gives me the first location not in range, but the array slice below
        # fixed that
        end=bisect.bisect_right(self.dPos[c], to)
        assert(start<=end)
        poss = self.dPos[c][start:end]
        annotations = [self.dAnnotations[c][pos] for pos in poss]
        return poss, annotations

    def getAnnotationsAsDict(self, c, frm, to):
        d={}
        frm=int(frm)
        to=int(to)
        assert(frm<=to)
        if not self.dPos.has_key(c):
            return ([],[])
        start = bisect.bisect_left(self.dPos[c], frm)
        # This actually gives me the first location not in range, but the array slice below
        # fixed that
        end=bisect.bisect_right(self.dPos[c], to)
        assert(start<=end)
        poss = self.dPos[c][start:end]
        for pos in poss:
            d[pos]=self.dAnnotations[c][pos]
        return d

    def getAnnotationsGenerator(self, c, frm=0, to=99999999999):
        frm=int(frm)
        to=int(to)
        assert(frm<=to)
        start = bisect.bisect_left(self.dPos[c], frm)
        # This actually gives me the first location not in range, but the array slice below
        # fixed that
        end=bisect.bisect_right(self.dPos[c], to)
        assert(start<=end)
        for i in xrange(start, end):
            pos=self.dPos[c][i]
            yield (c, pos, self.dAnnotations[c][pos])

    def load(self, f, hasHeader=False):
        self.f = f
        fp = open(f)
        if hasHeader:
            self.header = fp.readline()
        for l in fp:
            c, pos, a = self.process(l)
            if c:
                self.dPos.setdefault(c,[]).append(pos)
                self.dAnnotations.setdefault(c,{})[pos]=a
        for k,v in self.dPos.iteritems():
            v.sort()

    
class TestHandler(annotationHandler):
    def process(self, l):
        v=l.rstrip('\n').split('\t')
        c, type, s, e, genep = v[0:5]
        return c, int(s), v
        

        
