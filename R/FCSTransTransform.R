######
#Transformation code from FCSTrans imported by Greg Finak
# finak@fhcrc.org

##FCSTrans Authors
## Authors: Yue Liu and Yu "Max" Qian
# Contact: yu.qian@utsouthwestern.edu or yliu0@yahoo.com
 


# set output to 0 when input is less than cutoff value
ipfloor <- function (x, cutoff = 0, target = 0) {
  y = x
  if (x <= cutoff) y = target
  y
}

# set output to 0 when input is less than cutoff value
ipceil <- function (x, cutoff = 0, target = 0) {
  y = x
  if (x >= cutoff) y = target
  y
}

# calculation core of iplogicle
iplogicore <- function (x, w, r, d, scale,rescale=TRUE) {
  tol = .Machine$double.eps^0.8
  maxit = as.integer(5000)
  d = d * log(10)
  scale = scale / d
  p = if (w == 0) 
    1
  else uniroot(function(p) -w + 2 * p * log(p)/(p + 1), c(.Machine$double.eps, 
                                                          2 * (w + d)))$root
  a = r * exp(-(d - w))
  b = 1
  c = r * exp(-(d - w)) * p^2
  d = 1/p
  f = a * (p^2 - 1)
  y = .Call("biexponential_transform", as.numeric(x), a, b, c, d, f, w, tol, maxit)
  if(!rescale){
	scale<-1
  }
  y = sapply(y * scale, ipfloor)
  y
}

# function for calculating w 
iplogiclew <- function (w, cutoff = -111, r = 262144, d = 4.5, scale = 1) {
  if (w > d) 
    w = d
  y = iplogicore(cutoff, w, r, d, scale) - .Machine$double.eps^0.6
  y
}

# immport logicle function - convert fluorescent marker values to channel output
# rescale parameter will transform the data to cover the range, otherwise it will be on the scale of the actual transformation.
iplogicle <- function (x, r=262144, d=4.5, range=4096, cutoff=-111, w=-1, rescale=TRUE) {
  if (w > d) 
    stop("Negative range decades must be smaller than total number of decades")
  if (w < 0)
    w = uniroot(iplogiclew, c(0, d), cutoff=cutoff)$root
  print(paste("params: r=", r, "d=", d, "range=", range, "cutoff=", cutoff, "w=", w))
  y = iplogicore(x, w, r, d, range,rescale=rescale)
  y
}

#FlowCore Transformation using FCSTrans
FCSTransTransform<-function(transformationId="defaultFCSTransTransform",
	channelrange=16777215,
	channeldecade=7.224719870049579,rescale=TRUE)
		{
			k<-new("transform",
			.Data=function(x){
			x<-iplogicle(x,channelrange,channeldecade,rescale)
			});
			k@transformationId<-transformationId;
			k
		}

