Code for Intrinsic Dimension Estimators
================
2025-05-06

``` r
library(webshot2)
library(rlist)
library(rgl)
options(rgl.useNULL = TRUE)
library(reticulate)
use_python("C:/PythonFiles/PythonInstallFolder")
```

The following code presents one implementation of the ideas from
Kleindessner and von Luxburg in their methods to estimate the intrinsic
dimension of a data set. Their algorithm is set up to work without the
need for evaluating distances between the data points, so that only
relative ordering is required (A \< B, B \< C, etc), as opposed to the
usual Euclidean distance often used in standard dimension estimators
(\|A-B\|=0.5, \|B-C\|=0.3, etc).

The authors explain how often this relative ordering is easier to obtain
in place of exact values. As an example, they present the statement
“Movie A is more similar to movie B than to movie C” as contrast to the
statement “The similarity between A and B is 0.9 and the similarity
between A and C is 0.5”, the latter being much less natural for a person
to make.

However, we generate data sets and impose the ordering on the points by
computing these distances. We stress that this does not fundamentally
change the algorithm, but rather gives us a way of imposing the relative
order on the points.

We make several functions to generate data on spheres of dimensions 2,
3, and 4:

``` r
#generate data on a 2-sphere
sphere2 = function(n){
  
  y = rnorm(3*n)
  cols = list()
  
  for(i in 1:n){
    point = c(y[i],y[2*i],y[3*i])/
      sqrt(y[i]**2+y[2*i]**2+y[3*i]**2)
    cols[[i]] = point
  }
  
  return(cols)
  
}

#generate data on a 3-sphere
sphere3 = function(n){
  
  y = rnorm(4*n)
  cols = list()
  
  for(i in 1:n){
    point = c(y[i],y[2*i],y[3*i],y[4*i])/
      sqrt(y[i]**2+y[2*i]**2+y[3*i]**2+y[4*i]**2)
    cols[[i]] = point
  }
  
  return(cols)
  
}

#generate data on a 4-sphere
sphere4 = function(n){
  
  y = rnorm(5*n)
  cols = list()
  
  for(i in 1:n){
    point = c(y[i],y[2*i],y[3*i],y[4*i],y[5*i])/
      sqrt(y[i]**2+y[2*i]**2+y[3*i]**2+y[4*i]**2+y[5*i]**2)
    cols[[i]] = point
  }
  
  return(cols)
}
```

Let’s generate two data sets, sphere2 and sphere4:

``` r
spheres2 = sphere2(100)

spheres4 = sphere4(100)
```

We can visualize the points of the 2-dimensional sphere:

``` r
x_vals <- sapply(spheres2, function(v) v[1])
y_vals <- sapply(spheres2, function(v) v[2])
z_vals <- sapply(spheres2, function(v) v[3])

plot3d(x_vals, y_vals, z_vals, type = "s", col = "blue", size = 1)
rglwidget()
```

<img src="../../blob/main/DimEst0.jpg" width="672" />
The following code sets up all the components of the intrinsic dimension
estimator $E_{DP}$. Although the method works on a data set with just
relative ordering given, the distances are computed to impose this
relative ordering on the data points so as to generate the k-NN graph
information.

``` r
#distance function, dimension is of the ambient space
dist = function(x,y,dimension){
  
  s = 0
  
  for(i in 1:dimension){
    s = s + (x[i]-y[i])**2
  }
  s = sqrt(s)
  
  return(s)
}

#check all distances from given element
distCheck = function(element, set, dimension){
  
  element = unlist(element)
  distMatrix = list()
  
  for(i in 1:length(set)){
    x = unlist(set[i])
    d = dist(x,element,dimension)
    distMatrix = list.append(distMatrix,x,d)
  }
  
  return(distMatrix)
}

#order all elements by distance away from fixed element
distOrder = function(element, set, dimension){
  
  distMatrix = distCheck(element, set, dimension)
  d2 = list()
  
  for(i in 1:length(set)){
    d2 = append(d2,distMatrix[2*i])
  }
  
  d2 = sort(unlist(d2))
  orderedDistMatrix = list()
  
  for(j in 1:length(d2)){
    for(i in 1:length(set)){
      if(unlist(distMatrix[2*i])==d2[j]){
        orderedDistMatrix = 
          list.append(orderedDistMatrix,
                      unlist(distMatrix[2*i-1]),unlist(distMatrix[2*i]))
        break
      }
    }
  }
  return(orderedDistMatrix)
}

#find k-nearest neighbours and distances for the given element
knn = function(element, set, dimension, kvalue){
  
  #add 1 to kvalue to account for first element being the given element
  kvalue2 = 2*(kvalue+1)
  kCutOff = distOrder(element, set, dimension)[1:kvalue2]
  
  return(kCutOff)
  
}
  
#find the volume of the ball B_{SP}(i,2); this finds the number of elements 
#which are k-nearest neighbours to a fixed elements, an then adds the elements
#which are k-nearest neighbours to these, without repetition
knn2 = function(element, set, dimension, kvalue){
  
  #the radius 1 ball B_{SP}(i,1) always has volume k+1, we only need to add the 
  #missing elements in B_{SP}(i,2)
  count = kvalue+1
  
  toCheck = knn(element, set, dimension, kvalue)
  l = kvalue+1
  toCheck2 = list()
  
  for(i in 1:l){
    toCheck2[[i]] = unlist(toCheck[2*i-1])
  }
  
  for(m in 1:l){
    
    neighbours = knn(unlist(toCheck2[m]),set,dimension,kvalue)
    neighbours2 = list()
    
    #creates a workable list of new elements to check
    for(i in 1:l){
      neighbours2[[i]] = unlist(neighbours[2*i-1])
    }
    
    #checks whether the new list shares elements with the original, and only
    #counts those which do not, adding them to the total to avoid double counting
    for(i in 1:l){
      if(!(is.element(neighbours2[i],toCheck2))){
        count = count+1
        toCheck2 = append(toCheck2,neighbours2[i])
      }
    }
  }
  return(count)
}

#construct the volume ratio function used for the intrinsic dimension estimate,
#averaging over the first 10% of the elements in the set (so we need a data size
#of at least 11 elements)
LDP = function(set,dimension,kvalue){
  
  i = 1
  total = 0
  
  while(10*i<length(set)){
    x = unlist(set[10*i-1])
    total = total + (kvalue+1)/knn2(x,set,dimension,kvalue)
    i = i+1
  }
  total = total/i
  
  return(total)
}

#construct the intrinsic dimension estimator E_{DP} for the doubling property
#the dimension is again the ambient space dimension
EDP = function(set,dimension,kvalue){
  
  intDim = -log2(LDP(set,dimension,kvalue))
  
  return(intDim)
}
```

Now let us check the performance: experimentally we expect it to be very
poor, underestimating the dimension with slow convergence to the true
value. We set $k$ to be equal to $15$, following the authors.

``` r
EDP(spheres2,3,15)
```

    ## [1] 1.481723

``` r
EDP(spheres4,5,15)
```

    ## [1] 2.08091

As we can see, the dimensions are quite awful. Moreover, the estimator
seems not well able to work with higher dimensions. (And my computer
cannot run this algorithm with much higher values). Let’s try with
$n=175$

``` r
spheres22 = sphere2(175)
spheres42 = sphere4(175)
```

which give the dimension estimates:

``` r
EDP(spheres22,3,15)
```

    ## [1] 1.654715

``` r
EDP(spheres42,5,15)
```

    ## [1] 2.137223

Now we move onto the second estimator $E_{CAP}$. Define here the
necessary functions:

``` r
#check elements in intersection B(i,1), B(j,1): the k-nearest neighbours in 
#common for both elements
intersect = function(element1,element2,set,dimension,kvalue){
  
  B1 = list()
  B2 = list()
  B3 = list()

  knn1 = knn(element1,set,dimension,kvalue)
  knn2 = knn(element2,set,dimension,kvalue)
  
  l = kvalue+1
  
  for(i in 1:l){
    B1[[i]] = unlist(knn1[2*i-1])
    B2[[i]] = unlist(knn2[2*i-1])
  }
  #new list for keeping track of intersection
  for(i in 1:l){
      if(is.element(B1[i],B2)){
        B3 = append(B3,B1[i])
    }
  }
  
  B3Vol = length(B3)
  return(B3Vol)
}

#minimize the intersect function volume over all k-nearest neighbours of a
#given element
minimize = function(element, set, dimension, kvalue){
  
  #set up the list of elements to check, and remove the distances from knn
  knnList = knn(element,set,dimension,kvalue)
  knnList2 = list()
  
  l = kvalue+1
  
  for(i in 1:l){
    knnList2[[i]] = unlist(knnList[2*i-1])
  }
  
  #initialize list of all intersect function values
  volumes = list()
  
  #compute the intersect volumes, store and sort
  for(i in 1:l){
    vol = intersect(element, knnList2[i],set,dimension,kvalue)
    volumes = append(volumes,vol)
  }
  volumes = sort(unlist(volumes))
  
  #take smallest volume
  
  minimum = volumes[1]
  
  return(minimum)
}

#construct the volume ratio of intersection to B_{SP}(i,1) used in the estimator
#E_{CAP}, and average over 10% of the data set (so again need data size
#n at least 11)

LCAP = function(set,dimension,kvalue){
  
  i = 1
  total = 0
  
  while(10*i<length(set)){
    x = unlist(set[10*i-1])
    total = total + minimize(x,set,dimension,kvalue)/(kvalue+1)
    i = i+1
  }
  total = total/i
  
  return(total)
}

#define the incomplete regularized beta function with the necessary 
#parameters fixed
IRBeta = function(z){
  pbeta(0.75,(z+1)/2,0.5)
}

#define the inverse of the IRBeta function over the interval [0,10], since we 
#will not need to check dimensions outside of this range
IRBetaInv = function(w){
  f = function(u){
    IRBeta(u)-w
  }
  uniroot(f,c(0,10))
}
```

Now we can take the above to create the estimator $E_{CAP}$:

``` r
#we define the E_{CAP} estimator using the inverse of the IRBeta function 
#alongside our L_{CAP} value
ECAP = function(set,dimension,kvalue){
  
  y = LCAP(set,dimension,kvalue)
  IntDim = IRBetaInv(y)[1]
  
  return(IntDim)
}
```

Finally, let us run the above tests with this estimator:

``` r
spheres2 = sphere2(100)
spheres4 = sphere4(100)
spheres22 = sphere2(175)
spheres42 = sphere4(175)
```

and the estimated dimensions are:

``` r
ECAP(spheres2,3,15)
```

    ## $root
    ## [1] 2.485068

``` r
ECAP(spheres4,5,15)
```

    ## $root
    ## [1] 3.600966

``` r
ECAP(spheres22,3,15)
```

    ## $root
    ## [1] 2.389047

``` r
ECAP(spheres42,5,15)
```

    ## $root
    ## [1] 3.79899

We see somewhat better results here, despite low values of $n$ and $k$.
We need to increase both $n$ and $k$ values to see improvement, but due
to technical limitations, this is as far as my estimates can currently
go. The occasional failure to get closer to the true result with higher
$n$ and $k$ values is likely due to the fact that they are still very
small. Indeed, the tables from Kleindessner and von Luxburg about the
implementation of the code show that much higher values are required for
useful results.

We can try the algorithm on other data sets, for example the S-curve
from Python’s sklearn library, to see how well we can approximate the
intrinsic dimension on a data set we did not generate.

``` python
import sklearn
```

We generate the S-curve and visualize it below:

``` python
from sklearn.datasets import make_s_curve
X, t = make_s_curve(noise=0.05, random_state=0)
```

``` r
X_array = as.array(py$X)
```

To make the data convenient to use with the functions above, we convert
it to the same format of a list:

``` r
X_list = list()
for(i in 1:100){
    point = c(X_array[i,1],X_array[i,2],X_array[i,3])
    X_list[[i]] = point
  }
```

``` r
x_vals <- sapply(X_list, function(v) v[1])
y_vals <- sapply(X_list, function(v) v[2])
z_vals <- sapply(X_list, function(v) v[3])

plot3d(x_vals, y_vals, z_vals, type = "s", col = "blue", size = 1)
rglwidget()
```

    ## file:///C:/Users/barts/AppData/Local/Temp/RtmpGG5swG/file381469a02759.html screenshot completed

<img src="../../Users/barts/AppData/Local/Temp/RtmpGG5swG/file38141b483d7.png" width="672" />
We will apply our EDP and ECAP dimension estimators, using a lower value
of k=10. Since our data set only contains 100 points, a lower k value
seems to produce more accurate results.

``` r
EDP(X_list,3,10)
```

    ## [1] 1.232202

``` r
ECAP(X_list,3,10)
```

    ## $root
    ## [1] 2.102429

Indeed, the better ECAP estimator accurately sees that the intrinsic
dimension of the data set is 2, despite the S shape curving in on itself
a little bit.
