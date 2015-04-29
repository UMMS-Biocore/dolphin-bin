
d=read.table("testextended_forward.out", header=T, skip=2, sep="\t")
dt=t(d[,1:4])
buckets=ncol(dt)
l=list()
nms=list()
for (i in seq(buckets)) {
  l[[i]]=dt[,i]
  nms[[i]]=d[i,1]
}
boxplot(l, names=nms, ylim=c(0,1000))


