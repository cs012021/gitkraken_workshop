x <- 3
y <- 5
y
x+y
number <- x+y
x<-5
number
y<-10
number
number <- x+y
number

# Create a numeric vector and store the vector as a variable called 'glengths'
glengths <- c(4.6, 3000, 50000)
glengths

# Create a character vector and store the vector as a variable called 'species'
species <- c("ecoli", "human", "corn")
species

# Forget to put quotes around corn
species <- c("ecoli", "human", corn)

# Create a character vector and store the vector as a variable called 'species'
species <- c("ecoli", "human", "corn")

# Create a character vector and store the vector as a variable called 'expression'
expression <- c("low", "high", "medium", "high", "low", "medium", "high")

# Create a data frame and store it as a variable called 'df'
df <- data.frame(species, glengths)
# Turn 'expression' vector into a factor
expression <- factor(expression)

samplegroup <- c("CTL1", "CTL2", "CTL3", "KO1", "KO2", "KO3", "OE1", "OE2", "OE3")

samplegroup <- factor (samplegroup)

df1 <- data.frame(titles, pages)

list1 <- list(species, df, number)

list2 <- list(species, glengths, number)

glengths <- c(glengths, 90) # adding at the end	
glengths <- c(30, glengths) # adding at the beginning

sqrt(81)

sqrt(glengths)

round(3.14159)

?round

args(round)

example("round")

round(3.14159, digits=2)

round(3.14159, 2)

average(glengths)

mean(glengths)
?mean

test <- c(1, NA, 2, 3, NA, 4)
mean(test)
mean(test, na.rm=TRUE)

multiply_it <- function(x,y) {
  multiply <- x * y
  return(multiply)
}
multiply_it(4,6)

metadata <- read.csv(file="data/mouse_exp_design.txt")

# First column as row names
proj_summary <- read.table(file.choose(), header=T, row.names=1)

head(metadata)
str(metadata)
summary(metadata)
class(glengths)
class(metadata)
summary(proj_summary)
length(samplegroup)
dim(proj_summary)
rownames(metadata)
colnames(proj_summary)
length(colnames(proj_summary))

f_to_c <- function(temp_f) {
  temp_c  <- (temp_f - 32) * 5 / 9
  return(temp_c)
}
 
f_to_c(32)

age <- c(15, 22, 45, 52, 73, 81)  

age > 50 | age < 18

age

age[age > 50 | age < 18]  ## nested

# OR

## create a vector first then select
idx <- age > 50 | age < 18
age[idx]

which(age > 50 | age < 18)

age[which(age > 50 | age < 18)]  ## nested

# OR

## create a vector first then select
idx_num <- which(age > 50 | age < 18)
age[idx_num]

expression[expression == "high"]

sessionInfo()

install.packages("ggplot2")

library(ggplot2)

# Extract value 'Wt'
metadata[1, 1]

# Extract value '1'
metadata[1, 3] 

# Extract third row
metadata[3, ] 

# Extract third column
metadata[ , 3]   

# Extract third column as a data frame
metadata[ , 3, drop = FALSE] 

# Dataframe containing first two columns
metadata[ , 1:2] 

# Data frame containing first, third and sixth rows
metadata[c(1,3,6), ] 

# Extract the celltype column for the first three samples
metadata[c("sample1", "sample2", "sample3") , "celltype"] 

# Check column names of metadata data frame
colnames(metadata)

# Check row names of metadata data frame
rownames(metadata)

# Extract the genotype column
metadata$genotype 

# Extract the first five values/elements of the genotype column
metadata$genotype[1:5]

metadata[c(4,9), "replicate"]

metadata[ , "replicate", drop = FALSE] 

metadata$celltype == "typeA"

logical_idx <- metadata$celltype == "typeA"

metadata[logical_idx, ]

which(metadata$celltype == "typeA")

idx <- which(metadata$celltype == "typeA")
metadata[idx, ]

idx <- which(metadata$replicate > 1)

metadata[idx, ]

metadata[which(metadata$replicate > 1), ]

sub_meta <- metadata[which(metadata$replicate > 1), ]

metadata$genotype=="KO"
ko_meta <- metadata[which(metadata$genotype=="KO"),]
ko_meta

list1[[2]]

comp2 <- list1[[2]]
class(comp2)

list1[[1]][1]

random <- list(metadata, age, list1, samplegroup, number)
random[[4]]

names(list1) <- c("species", "df", "number")

names(list1)

names(random) <- c("metadata", "age", "list1", "samplegroup", "number")
random$age

rpkm_data <- read.csv("data/counts.rpkm.csv")
head(rpkm_data)

ncol(rpkm_data)
nrow(metadata)

A <- c(1,3,5,7,9,11)   # odd numbers
B <- c(2,4,6,8,10,12)  # even numbers

# test to see if each of the elements of A is in B	
A %in% B

A <- c(1,3,5,7,9,11)   # odd numbers
B <- c(2,4,6,8,1,5)  # add some odd numbers in 

A %in% B

intersection <- A %in% B
intersection

A[intersection]

any(A %in% B)
all(A %in% B)

B %in% A

intersection2 <- B %in% A

intersection2

B[intersection2]


A <- c(10,20,30,40,50)
B <- c(50,40,30,20,10)  # same numbers but backwards 

# test to see if each element of A is in B
A %in% B

# test to see if each element of A is in the same position in B
A == B

# use all() to check if they are a perfect match
all(A == B)

x <- rownames(metadata)
y <- colnames(rpkm_data)

all(x %in% y)

all(rownames(metadata) %in% colnames(rpkm_data))

x == y
all(x == y)

important_genes <- c("ENSMUSG00000083700", "ENSMUSG00000080990", "ENSMUSG00000065619", "ENSMUSG00000047945", "ENSMUSG00000081010", "ENSMUSG00000030970")

all(important_genes %in% rpkm_data)

rpkm_data[important_genes,]

teaching_team <- c("Jihe", "Mary", "Meeta", "Radhika")

# Extracting values from a vector
teaching_team[c(2, 4)] 

teaching_team

# Extracting values and reordering them
teaching_team[c(4, 2)] 

# Extracting all values and reordering them
teaching_team[c(4, 2, 1, 3)]

# Saving the results to a variable
reorder_teach <- teaching_team[c(4, 2, 1, 3)] 




reorder_second [c(4,1,5,2,3)]

match(first,second)

# Saving indices for how to reorder `second` to match `first`
reorder_idx <- match(first,second) 

# Reordering the second vector to match the order of the first vector
second[reorder_idx]

# Reordering and saving the output to a variable
second_reordered <- second[reorder_idx] 

first <- c("A","B","C","D","E")
second <- c("D","B","A")

match(first,second)
second[match(first, second)]

# Check row names of the metadata
rownames(metadata)

# Check the column names of the counts data
colnames(rpkm_data)

genomic_idx <- match(rownames(metadata), colnames(rpkm_data))
genomic_idx

# Reorder the counts data frame to have the sample names in the same order as the metadata data frame
rpkm_ordered  <- rpkm_data[ , genomic_idx]

mean(rpkm_ordered$sample1)

library(purrr)  # Load the purrr

samplemeans <- map_dbl(rpkm_ordered, mean) 

# Named vectors have a name assigned to each element instead of just referring to them as indices ([1], [2] and so on)
samplemeans

# Check length of the vector before adding it to the data frame
length(samplemeans)

# Create a numeric vector with ages. Note that there are 12 elements here
age_in_days <- c(40, 32, 38, 35, 41, 32, 34, 26, 28, 28, 30, 32) 

# Add the new vector as the last column to the new_metadata dataframe
new_metadata <- data.frame(metadata, samplemeans, age_in_days) 

# Take a look at the new_metadata object
View(new_metadata)

library(ggplot2)

ggplot(new_metadata) +
  geom_point (aes(x=age_in_days,y=samplemeans, color=genotype, shape=celltype))

# Make the data points bigger

ggplot(new_metadata) +
  geom_point (aes(x=age_in_days,y=samplemeans, color=genotype, shape=celltype), size=2.25) +
  theme_bw() +
  theme(axis.title=element_text(size=rel(1.5))) +
  xlab("Age(days)") +
  ylab("Mean expression") +
  ggtitle("New plot") +
  theme (plot.title=element_text(hjust=0.5))+
  scale_color_manual(values=c("purple","orange"))

theme_bw() +
  theme(axis.title=element_text(size=rel(1.5))) +
  theme(plot.title=element_text(size=rel(1.5), hjust=0.5))

personal_theme <- function(){
  theme_bw() +
    theme(axis.title=element_text(size=rel(1.5))) +
    theme(plot.title=element_text(size=rel(1.5), hjust=0.5))
}

ggplot(new_metadata) +
  geom_point(aes(x=age_in_days, y=samplemeans, color=genotype, fill=genotype, shape=celltype), size=rel(3.0)) +
  xlab("Age (days)") +
  ylab("Mean expression") +
  ggtitle("Expression with Age") +
  personal_theme()

ggplot(new_metadata) +
  geom_boxplot(aes(x=genotype, y=samplemeans, fill=celltype)) +
  xlab("Genotype") +
  ylab("Mean expression") +
  ggtitle("Genotype differences in average gene expression") +
  theme_bw() +
  theme(plot.title=element_text(size=rel(1.5), hjust=0.5)) +
  scale_fill_manual(values=c("red","magenta"))

write.csv(sub_meta, file="data/subset_meta.csv")



# Save a vector to file as a single column
write(glengths, file="data/genome_lengths.txt", ncolumns = 1)

## Open device for writing
pdf("figures/scatterplot.pdf")

## Make a plot which will be written to the open device, in this case the temp file created by pdf()/png()
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,shape=celltype), size=rel(3.0)) 

## Closing the device is essential to save the temporary file created by pdf()/png()
dev.off()

random_numbers <- c(81, 90, 65, 43, 71, 29)

mean(random_numbers) %>% round (digits=3)

library(tidyverse)

# Read in the functional analysis results
functional_GO_results <- read_delim(file = "data/gprofiler_results_Mov10oe.csv", delim = "\t" )


animals_tb <- animals %>%
  rownames_to_column(var = "animal_names") %>%
  as_tibble()
                                    
ggplot(animals_tb) +
  geom_boxplot(aes(x=genotype, y=samplemeans, fill=celltype)) +
  xlab("Genotype") +
  ylab("Mean expression") +
  ggtitle("Genotype differences in average gene expression") +
  theme_bw() +
  theme(plot.title=element_text(size=rel(1.5), hjust=0.5)) +
  scale_fill_manual(values=c("red","magenta"))