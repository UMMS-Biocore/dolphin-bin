
tests: testmean testmedian testregion

testmean:
	python ../ACT.py --mean --radius=1000 --nbins=40 --output=testmean.out ch22_genome.pos PolII.chr22.sgr

testmedian:
	python ../ACT.py --median --radius=1000 --nbins=40 --output=testmedian.out ch22_genome.pos PolII.chr22.sgr

testregion:
	python ../ACT.py --mean --region --mbins=40 --radius=1000 --nbins=40 --output=testregion.out ch22_genome.pos PolII.chr22.sgr

