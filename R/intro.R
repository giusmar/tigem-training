# COMMENTS ----------------------------------------------------------------

# This is my first comment

hello

getwd() # get the working directory
setwd() # set the working directory

# MATHEMATICAL BASIC OPERATIONS -------------------------------------------

## addition
5 + 7 

## substraction
5 - 7 

## moltiplication
5 * 7

## division
5 / 7

## other operation
# ^ or ** exponentiation
# % percentage
# %% modulus
# %/% integer division

# CREATING OBJECTS --------------------------------------------------------

x <- 15
y <- 17
x + y

## number sequences
number_sequence <- 1:10

## character sequences
### WARNING! R is case sensitive!
char <- 'computer'
char_sequence <- c(char, 'monitor')

## mathematical operation
z = x + y

# BRACKETS ----------------------------------------------------------------

## In R we can use all types of brackets and each one works differently

## Round brackets
### They define a function
print()
view()
head()
getwd()
setwd()

## Square brackets
### They allows you to pick a specific portion of an object or dataset

number_sequence
number_sequence[1]
char_sequence[1]

View(iris)
head(iris)
iris[1,5]

## Curly brackets
### Used to write your own function (last section)

# DATA STRUCTURES ---------------------------------------------------------

## Types of variables
### Real number, Integer number, Character and Logical operators

## Vectors
number_vector <- c(1,2,4,6,10)
character_vector <- c('London', 'Paris', 'Roma', 'Madrid')
logical_vector <- c(TRUE, FALSE, TRUE, TRUE)
logical_vector <- c(T, F, T, T)

x <-seq(from = 0, to = 20, by = 2)
seq(1,20)

sum(x)
mean(x)
median(x+1)
x / 2

y <- c(1,2,3)
z <- c(4,5,6)

y + z
y * z
y / z

m <- c(1,2,3,4)
n <- c(5,6) #recycling
m + n

## error
m <- c(1,2,3,4)
n <- c(5,6,3)
m + n

sort(x)
character_vector[c(1,2)]

## Matrices
x <- 1:20
matrix(x)
matrix(x, nrow = 2)
matrix(x, nrow = 2, byrow = TRUE)
matrix(x, nrow = 20, ncol = 2) #recycling

mat_1 <- matrix(x, nrow = 4, byrow = TRUE)
mat_1
t(mat_1)

mat_2 <- matrix(x, 
               nrow = 4, 
               byrow = TRUE,
               dimnames = list(c("col_1", "col_2", "col_3", "col_4"),
                               c("row_1", "row_2", "row_3", "row_4", "row_5")))
mat_2

y <- 1:5
z <- 11:15
cbind(y,z)
cbind(y,z)[3,"z"]


## Array
array_1 <- array(1:24, c(2,3,4))
array_1

## Lists
list_1 <- list(
                11:15,
                c("T","T","F","T","T"),
                sum,
                matrix(1:9, nrow = 3))
list_1
list_1[[1]]
list_1[[1]][3]
list_1[[4]][2,2]

## Factors
ppl <- c("Carlo", "Marco", "Giulia", "Ettore", "Lucia")
factor(ppl)

val <- c("Positive", "Negative", "Positive", "Positive", "Negative")
factor(val)

factor(c("primary school","secondary school","high school","high school","secondary school","secondary school"),
       order = TRUE,
       levels = c("primary school","secondary school","high school"))

## Dataframes
iris
iris$Species
iris[,1:4]
iris

summary(iris)

library(tidyverse)
names(iris) <- tolower(names(iris))
virginica <- filter(iris, species == "virginica")
sepalLength6 <- filter(iris, species == "virginica", sepal.length > 6)
selected <- select(iris, sepal.length, sepal.width, petal.length)
selected2 <- select(iris, sepal.length:petal.length)
newCol <- mutate(iris, greater.half = sepal.width > 0.5 * sepal.length)
iris$greater.half <- iris$sepal.width > 0.5 * iris$sepal.length

# RELATIONAL OPERATORS ----------------------------------------------------

## They are helpful when you wanto to subset data with specific charateristic

## With number
1 < 2
1 <= 2
1 > 2
1 >= 2

## With character
"Hello" == "hello"
"Hello" != "hello"

## More complex situation
1 > 2 & 2 > 3 # true x true
1 < 2 & 2 < 3 # false x false
1 > 2 & 2 < 3 # true x false
1 > 2 | 2 < 3 # at least one true

ToF <- c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE)

sum(ToF)

# CONTROL STRUCTURES ------------------------------------------------------

## if and else: to test conditions
x <- 5

if (x < 0) {
  print("The number is smaller than 0")
}

if (x < 0) {
  print("The number is smaller than 0")
} else {
  print("The number is bigger than 0 ")
}

## for: repeat the instructions for a certain number
vect <- c(13,16,45,28,5,7,8,39)

for (i in vect) {
  print(i)
}

for (i in 1:length(vect)) {
  print(i)
}

for (i in 1:length(vect)) {
  print(vect[i])
}

## while: repeat the instructions until a condition 
while (1+1 == 2) {
   print("ok")
 }

z <- 1

while (z < 5) {
  print("add 1")
  z <- z +1
  if (z == 5){
    print("stop")
  }
}

## break: break the execution of instructions 
x <- 1

while (x <= 10) {
  print(x*2)
  if ((x*2) > 20) {
    break
  }
  x <- x + 1 
}

## return: exit from a function and give a specific output

# FUNCTIONS ---------------------------------------------------------------

print()
print
sum
args(print)
help(print)
?print


# my_fun <- function(arg_1,arg_2) {
#   instruction
# }

fun_sum <- function(x,y) {
  x + y
}

fun_sum(1,5)

fun_1 <- function(x,y) {
  z <- x + y
  if (z < 20) {
    print(paste(z, "< 20"))
  } else {
    print(paste(z, "> 20"))
  }
}

fun_1(13,15)

## function on dataset

head(iris)
iris[1:10,1:4]

fun_dat <- function(dataset,col_seq) {
  for ( i in col_seq) {
    printsum <- sum(dataset[,i])
    print(paste("Sum of", names(iris)[i], "column is", sum(iris[,i])))
    }
}

fun_dat(iris,1:4)
