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

as("filter1","filter")
as(as.name("filter1"),"filter")
as(~ filter1 %subset% filter2,"filter")

fs[["Combined"]] = ~ filter1 %subset% `Live Cells`
as(fs,"list")
f = filter(b08,fs)
f
summary(f)



