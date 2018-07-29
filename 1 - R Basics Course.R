19.71 * .15
install.packages("dslabs")
library(dslabs)
library(tidyverse)
library(dplyr)
a = 1
b = 1
c = -1
#ls() shows all objects in workspace
sol_1 = (-b +sqrt(b^2-4*a*c))/(2*a)
#R does not impute "4c" as multiplication, need 4*c
sol_2 = (-b - sqrt(b^2-4*a*c))/(2*a)
#you ask for the arguments for a function with args()
args(log)
#list all datasets provided in R
data()
# print value with print()
print(sol_1)
# exp() and log() default to natural
exp(1)
#class() gives you object type
data("murders")
class(murders)
#str() gives details of obj
str(murders)
#names() gives names of columns
names(murders)
#length() gives # of entries in vector
length(murders$population)
#factors store categorical data (repeated characters)
#levels() gives factors in column
#using [[""]] rather than $ returns a data.frame rather than a vector
# determine if two things are identical with identical()
#table() returns the frequency of unique elements in a vector

#example for using vectors
codes <- c(380,124, 818)
country <- c("italy", "canada", "egypt")
codes <- c(italy = 380, canada = 124, egypt = 818)
codes
class(codes)
#names() can be used to assign names to entries in a vector
names(codes) = country
codes

#Subsetting
codes[2]

#Vector Coercion
#creating a vector with a string in it will lead R to coerce everything else to a string
#can force coercion with as.numeric(), as.character(), and as.factor
#failure to coerce => NA
#NA is a special data type for "Not Available"

#Integer Class
#a distinct class from numeric is the integer, can be created by adding "L" after a #
#takes up less space than a number

#Sorting
sort(murders$total)
#sort() orders the column
#order() produces an index which, when applied to the vector, would sort it
#e.g. index = order(x)
#x[index] returns x in ascending order

#one can order one column by the order of another by using the index of the first
index <-order(murders$total)
murders$abb[index]
#want to find which state has the most murders?
max(murders$total)
i_max = which.max(murders$total)
print(i_max)
murders$state[i_max]
#which.min accomplishes the same for the minimum
#rank() gives the rank of each number in vector (where the smallest element has rank 1)
# sort() == x[order(x)]

#Creating a data.frame
data.frame()
#data.frame(NameOfColumn1 = VectorObject1, NameofColumn2 = VectorObject2)
murders$state[which.max(murders$population)]

#arithmetic operations occur element-wise

#remember: indexing refers to ordering based on a certain criterion

# a useful way to create a new column is to create a var from manipulation of columns
#e.g. col1/col2 = col3; df$col3 = col3

#one can create an index by defining a var e.g. "index" as a logical

#you can then subset based on the logical with data[index]

#sum coerces TRUE to numeric i.e. 1, so sum() of logical gives total TRUE

#want to satisfy multiple conditions? create a variable for each, and then use
#the and operator "&" to create an index that contains each

#Indexing Functions
which() #gives entries for which something is true
match() #give it a vector of names and a column, and it returns the row # where they appear
# %in% checks to see if the elements of a first vector are in a second vector

#here's an example of using which() and %in% in concert to identify which components of one
#vector are present in another:
abbs <- c("MA", "ME", "MI", "MO", "MU")
ind <- which(!abbs%in%murders$abb)
abbs[ind]


#Dplyr Package
install.packages("dplyr")
library(dplyr)

#Important functions: mutate(), filter(), select()

mutate()

murders <-mutate(murders, rate= total/population*100000)
#note that total and population are not defined; mutate knows to look to column names

filter()

filter(murders, rate < 0.71)
#arg: data, conditional

select()

new_table <- select(murders,state, region,rate)
#args: data, var1, var2, var3

##The Pipe

#A simpler way to accomplish the following:
new_table <- select(murders,state, region,rate)
filter(new_table, rate < 0.71)

murders %>%  select(state, region, rate) %>% filter(rate > 0.71)

#The Data Frame Function

date.frame()
grades <- data.frame(names = c("John", "Juan", "Jean", "Yao"),
                     exam_1 = c(95, 80, 85, 92),
                     exam_2 = c(93, 86, 95, 93), stringsAsFactors = FALSE)
#stringsAsFactors used to prevent data.frame from turning characters into factors

#rank() gives ranks of x from lowest to highest, -x gives highest to lowest

genome %>% mutate(frequency = mutant/total) %>% select() %>% filter() 

library(dplyr)
require(dplyr)

#Basic Plots

population_in_millions <- murders$population/10^6
total_gun_murders <- murders$total
#basic plot
plot(population_in_millions, total_gun_murders)
#basic histogram
hist(murders$rate)
#which "state" has the highest murder rate?
murders$state[which.max(murders$rate)]
#how do the distributions of murder rates compare across regions?
boxplot(rate~region, data = murders)

#Super Basic Programming in R
a <- 2
if(a!=0){
  print(1/a)
}else{
  print("No reciprocal for 0.")
}

#General Form
#if(boolean_condition{
#  expressions
#}else{
# expressions
#})

#alternative function:
ifelse()
a <- 0
ifelse(a>0, 1/a, NA)
data(na_example)
require(dslabs)

any() #takes a logical and returns true if any are true

all() #takes a logical and returns true if all are true

#Basic Functions

avg <- function(x){
  s<- sum(x)
  n <- length(x)
  s/n
}

avg(murders$rate)

#variables defined inside a function are not defined outside, so
# I can define s within avg and then define s <- 3 outside avg without conflict

#General Form
my_function <- function(x){
  operations_on_x
  value_of_final_line_is_returned
}

#General Form with Argument
my_function <- function(x, arithmetic = TRUE){
  n <-length(x)
  ifelse(arithmetic = TRUE, sum(x)/n, prod(x)^(1/n))
}
#returns arithmetic mean if TRUE and geometric if false

#Using For Loops

compute_s_n <- function(n){
  x <- 1:n
  sum(x)
}
compute_s_n(4)


#General Form
for (i in range_of_values) {
  operations_that_use_i_which_is_changing_across_values
}
m = 25
#create an empty vector
s_n <- vector(length = m)
#create for loop to find product for 1:25
for(n in 1:m){
  s_n[n] <- compute_s_n(n)
}
n <- 1:25
plot(n, s_n)
lines(n, n*(n+1)/2)
args(lines)



#Conditionals
#any() and all() return logicals as you would predict

#More Functions
sum_n <- function(n){
  x <- 1:n
  y <- sum(x)
  y
}
#the above creates a function called sum(n) that creates a series from 1 to the input,
#and returns a sum


#If you want to fill a vector from a for-loop, you need to first create an empty vector
# with vector()

results <- vector("numeric", 10)
n <- 10
for(i in 1:n){
  x <- 1:i
  results[i] <- sum(x)
}

