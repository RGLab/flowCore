library(flowCore)

b08 = read.FCS(system.file("extdata","0877408774.B08",package="flowCore"),transformation="scale")
e07 = read.FCS(system.file("extdata","0877408774.E07",package="flowCore"),transformation="scale")
f06 = read.FCS(system.file("extdata","0877408774.F06",package="flowCore"),transformation="scale")


fs = new("filterSet")
fs[["filter1"]] = rectangleGate("FSC-H"=c(.2,.8),"SSC-H"=c(0,.8))
fs[["filter1"]]
fs["filter1"]
sum(b08 %in% fs["filter1"])
b08.result1 = filter(b08,fs["filter1"])
summary(b08.result1)

fs[[""]] = norm2Filter("FSC-H","SSC-H",scale.factor=2,filterId="Live Cells")

filter1 = rectangleGate("FSC-H"=c(.2,.8),"SSC-H"=c(0,.8))
filter2 = norm2Filter("FSC-H","SSC-H",scale.factor=2,filterId="Live Cells")

format(filter1)


as("filter1","filter")
as(as.name("filter1"),"filter")
as(~ filter1 %subset% filter2,"filter")

fs[["Combined"]] = ~ `Live Cells` %subset% filter1
as(fs,"list")
fs
sort(fs,dependencies=TRUE)
f = filter(b08,fs)
f
as.data.frame(f)
f@subSet[1:10,]
colSums(f@subSet)
summary(f[[1]])
summary(f[[2]])
summary(f[[3]])
rownames(f@dependency)[rowSums(f@dependency)==0]

split(b08,f,flowSet=TRUE)
split(b08,f,drop=TRUE,flowSet=TRUE)



