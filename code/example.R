## PALMER LAB AIL BREEDER SELECTION CODE
## written by Dr. Andrew Skol (askol at uchicago.edu) 

#######################################################################

## This script selects breeders from an advanced intercross line 
## of mice using kinship data based on a pedigree. This is the script
## you will use to obtain a list of breeder pairs for each generation. 

## It calls on two other scripts: kinship.R, which calculates pairwise
## kinship coefficients for members of the current generation, and 
## find_mates.R, which pairs breeders with the smallest average inbreeding
## coefficients. find_mates.R also tries to pair animals that will produce 
## the minimum average inbreeding coefficient in the offspring. 

#########################################################################

## This is an example used to select breeders for generation 10 of 
## the LGxSM AIL. 

## First set the working directory (see function "setwd") to the
## location where the code and scripts are.
source("kinship.R")
source("find_mates.R")
data.dir = "../data/"
ped = c()

## Construct the matrix ped by looking through all files in this
## directory. All files have the same name format as described below,
## the loop is reading through files for generations 2-10. The file
## contains information on potential breeders.
for (i in 2:10){
  file.name = paste(data.dir,"ped",i,".csv",sep="")
   
  ## Here we assume that missing values are written as question marks.
  ped.tmp = read.table(file=file.name,sep=",", skip = 1,
                       as.is=TRUE, na.strings="?")
  
  ## the next step assumes that col 1 in your csv file is the
  ## mouse id. col 8 = dam id, coln 7 = sire id, and col 3 = sex.
  ped.tmp = ped.tmp[,c(1,8,7,3)]
  dim(ped.tmp)
  if (i == 2){
    ped.tmp[,c(2,3)] = 0}
 
  ## This step adds a decimal onto all IDs in generation 34 and beyond, which
  ## was necessary because after inheriting the LGxSM AIL in generation 34, 
  ## we realized that some of our ear tag IDs overplapped with Jim Chevreud's. 
  ## Since this example is for F2-4, it is commented out, but may be useful for
  ## similar situations.
  # if (i == 34){
  #   idx = (ped.tmp[,1] - ped.tmp[,1]%/%1) == 0
  #   ped.tmp[idx,1] = ped.tmp[idx,1]+.1 }
  
  ## Check that parents are in previous generation. In this example it
  ## will return five warning messages from generation 3.
  if (i != 2){
    tmp = check.ped(ped.tmp , ped.prev)           
    if (length(tmp) > 0){
       print(paste("Generation:",i,"Parents not present in previous generation:",tmp))
       }
   }
    ped.prev = ped.tmp	
  ped = rbind(ped,ped.tmp)
}

## Change sex coding from alphabetical to numeric.
ped[,4] = as.character(ped[,4])
ped[ped[,4]=="F",4] = 0
ped[ped[,4]=="M",4] = 1
ped[,4] = as.numeric(ped[,4])
ped = data.matrix(ped)

## Call kinship function & build the covariance matrix. This can take several
## minutes if the pedigree is large. 
k = kinship(ped)

N = dim(k)[1]
no = dim(ped.tmp)[1]
idx = (N-no+1):N
k.prev = k[idx,idx]
ped.prev = ped[idx,]

## If you need it again and don't want to run it twice:
save.image("F10.RDataStart")

## Call the mate pairing script. Alternatively, you may use the
## function find.mates.res, which uses only one individual per sibship. 
out = find.mates(ped.prev,k.prev)
colnames(out) = c("dam", "sire", "kinship")
write.table(file="F10pairs.txt", out)

## This usually outputs the 30 best pairs or so. Normally
## you will want more than this. Run it 4 more times or
## until it stops giving you new pairs. The resulting file will contain 
## three unlabeled columns - Dam, Sire, Kinship.

##  Select a second set of breeders.
breed.rslts = out
first.os = max(ped.prev[,1])+1000
ped.next= cbind(first.os:(first.os+(dim(breed.rslts)[1]-1)),
  breed.rslts[,1:2], rep(1,dim(breed.rslts)[1]))
ped.prev = as.data.frame(ped.prev, )
ped.next = as.data.frame(ped.next, )
names(ped.next)[1:4] = c("ID","Mother","Father","Sex")
names(ped.prev) = c("ID","Mother","Father","Sex")
ped.comb = rbind(ped.prev, ped.next)


## Update kinship matrix 
k.comb = kinship.update(ped.prev , k.prev, ped.next)
pos.prev = 1:dim(ped.prev)[1]
pos.next = dim(ped.prev)[1]+(1:dim(ped.next)[1])
not.avail.idx = ped.prev[,1] %in% c(ped.next[,2],ped.next[,3])
not.avail.pos = c(1:dim(ped.prev)[1])[not.avail.idx]
new.mates = find.mates.given.pop(ped.comb, k.comb, 
                                 pos.prev, pos.next, not.avail.pos)

write.table(rbind(out, new.mates[,1:3]), file="F10pairs.txt", 
            quote=FALSE, row.names=FALSE)

## If you get the following error - Error: (list) object cannot be coerced
## to type 'double' - it means one or more ids are duplicated

## Note that every time after the second set is selected, you must
## change a few things in the code to make it work. I've marked these lines.

## Select a third set of breeders. 
breed.rslts = rbind(out , new.mates[,1:3])    ## here
ped.next= cbind(first.os:(first.os+(dim(breed.rslts)[1]-1)),
  breed.rslts[,1:2], rep(1,dim(breed.rslts)[1]))
ped.prev = as.data.frame(ped.prev, )
ped.next = as.data.frame(ped.next, )
names(ped.next)[1:4] = c("ID","Mother","Father","Sex")
names(ped.prev) = c("ID","Mother","Father","Sex")
ped.comb = rbind(ped.prev, ped.next)

## update kinship matrix 
k.comb = kinship.update(ped.prev , k.prev, ped.next)
pos.prev = 1:dim(ped.prev)[1]
pos.next = dim(ped.prev)[1]+(1:dim(ped.next)[1])
not.avail.idx = ped.prev[,1] %in% c(ped.next[,2],ped.next[,3])
not.avail.pos = c(1:dim(ped.prev)[1])[not.avail.idx]

## here the df has a new name. you could also overwrite the old one.
newer.mates = find.mates.given.pop(ped.comb, k.comb, 
                                   pos.prev, pos.next, not.avail.pos)

all.breeders=rbind(out, new.mates[,1:3], newer.mates[,1:3])   ## here

write.table(all.breeders, file="F10pairs.txt", 
            quote=FALSE, row.names=FALSE)   ## here

## .... and so on. repeat until you're satisfied.

