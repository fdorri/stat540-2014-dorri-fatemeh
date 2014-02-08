Seminar3-Summary of R- Learning
========================================================

R data type:
----
**Vector**

```r
a <- c(1, 2, 5.3, 6, -2, 4)  # numeric vector
b <- c("one", "two", "three")  # character vector
c <- c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE)  #logical vector
length(c)  #length of a vector
```

```
## [1] 6
```


```r
a[c(2, 4)]  # 2nd and 4th elements of vector
```

```
## [1] 2 6
```

**Matrices**

All columns in a matrix must have the same mode(numeric, character, etc.) and the same length.

```r
# generates 5 x 4 numeric matrix
y <- matrix(1:20, nrow = 5, ncol = 4)
# another example
c <- c(1, 4, 15, 7, 81, 19, 23, 5, 36, 47)
cells <- c(1, 3, 6, 4)
rnames <- c("R1", "R2")
cnames <- c("C1", "C2")
mymatrix <- matrix(cells, nrow = 2, ncol = 2, byrow = TRUE, dimnames = list(rnames, 
    cnames))
mymatrix
```

```
##    C1 C2
## R1  1  3
## R2  6  4
```

```r
dimnames(mymatrix)
```

```
## [[1]]
## [1] "R1" "R2"
## 
## [[2]]
## [1] "C1" "C2"
```

```r
dim(mymatrix)  #size of a matrix
```

```
## [1] 2 2
```

```r
nrow(mymatrix)
```

```
## [1] 2
```

```r
ncol(mymatrix)
```

```
## [1] 2
```


```r
n <- 4
B <- 3
x <- matrix(rnorm(n * B), nrow = n)
rownames(x) <- sprintf("obs%02d", 1:n)
colnames(x) <- sprintf("samp%02d", 1:B)
x
```

```
##        samp01  samp02  samp03
## obs01 -0.2790  0.6884  0.1013
## obs02  0.2060 -0.3385 -0.3096
## obs03  1.6901  0.2715 -0.9239
## obs04  0.6962 -1.7896 -1.0511
```


**Arrays**

Arrays are similar to matrices but can have more than two dimensions.

**Data Frame**

```r
d <- c(1, 2, 3, 4)
e <- c("red", "white", "red", NA)
f <- c(TRUE, TRUE, TRUE, FALSE)
mydata <- data.frame(d, e, f)
names(mydata) <- c("ID", "Color", "Passed")  # variable names
```

**List**

An ordered collection of objects (components). A list allows you to gather a variety of (possibly unrelated) objects under one name.

```r
# example of a list with 4 components - a string, a numeric vector, a
# matrix, and a scaler
w <- list(name = "Fred", mynumbers = a, mymatrix = y, age = 5.3)

# example of a list containing two lists v <- c(list1,list2)
```

**Some useful functions!**
```
c(object,object,...)       # combine objects into a vector
cbind(object, object, ...) # combine objects as columns
rbind(object, object, ...) # combine objects as rows 

length(object) # number of elements or components
str(object)    # structure of an object 
class(object)  # class or type of an object
names(object)  # names

ls()       # list current objects
rm(object) # delete an object
```
**Numeric functions**
```

Function  Description
abs(x)  absolute value
sqrt(x)	square root
ceiling(x)	ceiling(3.475) is 4
floor(x)	floor(3.475) is 3
trunc(x)	trunc(5.99) is 5
round(x, digits=n)	round(3.475, digits=2) is 3.48
signif(x, digits=n)	signif(3.475, digits=2) is 3.5
cos(x), sin(x), tan(x)	also acos(x), cosh(x), acosh(x), etc.
log(x)	natural logarithm
log10(x)	common logarithm
exp(x)	e^x
```
Control structure
-----
**Logical Operators**

- and: &
- or: |
- not: !

**IF**

```r
# if(cond1=true) { cmd1 } else { cmd2 }
if (1 == 0) {
    print(1)
} else {
    print(2)
}
```

```
## [1] 2
```

**for**

```r
# for(variable in sequence) { statements }
x <- 1:10
z <- NULL
for (i in seq(along = x)) {
    if (x[i] < 5) {
        z <- c(z, x[i] - 1)
    } else {
        z <- c(z, x[i]/x[i])
    }
}
z
```

```
##  [1] 0 1 2 3 1 1 1 1 1 1
```

```r
# Example 2
n <- 3
B <- 10
x <- rnorm(n)
for (j in 1:(B - 1)) {
    x <- cbind(x, rnorm(n))
}
```

**while**

```r
# while(condition) statements
z <- 0
while (z < 5) {
    z <- z + 2
    print(z)
}
```

```
## [1] 2
## [1] 4
## [1] 6
```

**switch**
```
switch(expr, ...)
```
**ifelse**

```r
# ifelse(test, true_value, false_value)
x <- 1:10  # Creates sample data
ifelse(x < 5 | x > 8, x, 0)
```

```
##  [1]  1  2  3  4  0  0  0  0  9 10
```


```r
# Example transpose of a matrix a poor alternative to built-in t() function

mytrans <- function(x) {
    if (!is.matrix(x)) {
        warning("argument is not a matrix: returning NA")
        return(NA_real_)
    }
    y <- matrix(1, nrow = ncol(x), ncol = nrow(x))
    for (i in 1:nrow(x)) {
        for (j in 1:ncol(x)) {
            y[j, i] <- x[i, j]
        }
    }
    return(y)
}
z <- matrix(1:10, nrow = 5, ncol = 2)
tz <- mytrans(z)
```


Statistical functions
-------
```
Function  Description
mean(x, trim=0,
na.rm=FALSE)	mean of object x
# trimmed mean, removing any missing values and 
# 5 percent of highest and lowest scores 
mx <- mean(x,trim=.05,na.rm=TRUE)
sd(x)	standard deviation of object(x). also look at var(x) for variance and mad(x) for median absolute deviation.
median(x)	median
quantile(x, probs)	quantiles where x is the numeric vector whose quantiles are desired and probs is a numeric vector with probabilities in [0,1].
# 30th and 84th percentiles of x
y <- quantile(x, c(.3,.84))
range(x)	range
sum(x)	sum
diff(x, lag=1)	lagged differences, with lag indicating which lag to use
min(x)	minimum
max(x)	maximum
scale(x, center=TRUE, scale=TRUE)	column center or standardize a matrix.
```

List of distributions:
---------
**Uniform distribution**

```r
n <- 10
runif(n)
```

```
##  [1] 0.56303 0.24342 0.66645 0.02288 0.47225 0.01604 0.59917 0.47697
##  [9] 0.43123 0.26728
```

**Normal distribution**

```r
rnorm(n, sd = 1)
```

```
##  [1] -0.5330 -1.3024  0.2791 -1.5271 -0.5358 -0.6411  1.7031 -0.2225
##  [9] -0.6187  0.1321
```

**Poisson distribution**

```r
lambda <- 0.5
rpois(n, lambda)
```

```
##  [1] 1 0 0 1 1 1 3 1 2 0
```

**Binomial Distribution**

```r
size <- 3
p <- 0.5
rbinom(n, size, p)
```

```
##  [1] 1 2 3 0 3 1 2 1 2 2
```

```
name------description
dname-----density or probability function
pname-----cumulative density function
qname-----quantile function
rname-----random deviates
```

```r
x <- 6
dpois(x, lambda, log = FALSE)
```

```
## [1] 1.316e-05
```


