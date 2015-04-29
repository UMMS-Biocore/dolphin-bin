
import getAnnotations

iupac2code={
'A':1,
'C':2,
'G':4,
'T':8,
'R':1|4,
'Y':2|8,
'S':2|4,
'W':1|8,
'K':4|8,
'M':1|2
}

code2iupac={
1:'A',
2:'C',
4:'G',
8:'T',
1|4:'R',
2|8:'Y',
2|4:'S',
1|8:'W',
4|8:'K',
1|2:'M'
}

def homo(ltr):
    n=iupac2code[ltr]
    return n==1 or n==2 or n==4 or n==8

def hetero(ltr):
    return not homo(ltr)

def comp(a,b):
    return iupac2code[a]^iupac2code[b]

def mutant(f,m,c):
    '''check if child has any alleles not present in the parents, indicating
    a mutation or sequencing error'''
    return iupac2code[c]&~(iupac2code[f]|iupac2code[m])

class Handler(getAnnotations.annotationHandler):
    def process(self, l):

        e=l.rstrip().split()
        chrom=e[1]
        pos=int(e[2])
        m_genotype=e[6]
        f_genotype=e[11]
        c_genotype=e[16]
        hetSNP=c_genotype in 'RYSWKM'
        ref=e[3]
        try:
            if e[20]=='SNP-1':
                m_allele, f_allele, typ = self.phase(m_genotype, f_genotype, c_genotype)
                return chrom, pos, (m_genotype, f_genotype, c_genotype, m_allele, f_allele, typ, ref, hetSNP)
            else:
                return (None, None, None)
        except:
            return (None, None, None)
                
    def phase(self, m, f, c):
        if mutant(f,m,c):
        #print "child mutant, skipping"
            return None, None, "mutant"
        if homo(c):
        #print "child homo, skipping"
            return None, None, "homo"
        if hetero(m) and hetero(f):
        #print "all hetero, skipping"
            return None, None, "hetero"
        if homo(m):
        # since we're here, we know child and father are het, so 
            m_allele=m
            f_allele=code2iupac[comp(c, m_allele)]
        else:
            # since we're here, we know child and mother are het, so 
            f_allele=f
            m_allele=code2iupac[comp(c, f_allele)]

        return m_allele, f_allele, "phased"

    def toDict(self, v):
        fields=('m_genotype', 'f_genotype', 'c_genotype', 'm_allele', 'f_allele', 'typ', 'ref', 'hetSNP')
        assert len(v)==len(fields)
        d=dict(zip(fields, v))
        return d
