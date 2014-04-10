kinship <- function(ped){

  ## ped should be a matrix with 4 rows in the follwing order:
  ## ID, Father, Mother, Sex and should be ordered from from founding
  ## generation to current.  Sex may be coded as "M"/"F" or as 1/0.
  ## 

  m.idx = ped[,4] == "M"
  f.idx = ped[,4] == "F"
  N = dim(ped)[1]
  
  if (sum(c(m.idx,f.idx)) == N){
    ped[m.idx,4] = 1
    ped[f.idx,4] = 0
  }else if (sum(c(m.idx,f.idx)) != 0) {
    sex.warning()
    return()
  } else if (sum(ped[,4] %in% c(0,1))<N){
    sex.warning()
    return()
  }

  ## Recode ids from 1:N ##
  ped.tmp = matrix(0,N,4)
  ped.tmp[,1] = 1:N
  for (i in 1:N){
    id.tmp = ped[i,1]
    idx = ped[,2]==id.tmp
    ped.tmp[idx,2] = i
    idx = ped[,3] == id.tmp
    ped.tmp[idx,3] = i
  }
  ped = ped.tmp
  
  ## Begin calculation ##
  kinship = diag(rep(1,N))
  for (i in 1:N){
    
    i.M = ped[i,2]
    i.F = ped[i,3]
    
    for (j in i:1){

      if (j==i){
        ## if j's parents are in the pedigree ##
        if (i.M!=0 & i.F!=0){          
          kinship[i,j] = 1 + .5*kinship[i.M,i.F]
        }
      } else {
        if (i.M!=0 & i.F!=0){
          kinship[i,j] = .5*kinship[i.M,j] + .5*kinship[i.F,j]
          kinship[j,i] = kinship[i,j]
        }
      }

    }
  }
  return(kinship)
}


kinship.update<-function(ped.old, k.ped.old, ped.new){

  #m.idx = ped.new[,4] == "M"
  #f.idx = ped.new[,4] == "F"
  N.old = dim(ped.old)[1]
  N.new = dim(ped.new)[1]
  N = N.old + N.new
  #if (sum(c(m.idx,f.idx)) == N.new){
  #  ped.new[m.idx,4] = 1
  #  ped.new[f.idx,4] = 0
  #}else if (sum(c(m.idx,f.idx)) != 0) {
  #  sex.warning()
  #  return()
  #} else if (sum(ped.new[,4] %in% c(0,1))<N.new){
  #  sex.warning()
  #  return()
  #}

  ped = rbind(ped.old,ped.new)
  
  ## Recode ids from 1:N ##
  ped.tmp = matrix(0,N,4)
  ped.tmp[,1] = 1:N
  for (i in 1:N){
    id.tmp = ped[i,1]
    idx = ped[,2]==id.tmp
    ped.tmp[idx,2] = i
    idx = ped[,3] == id.tmp
    ped.tmp[idx,3] = i
  }
  ped = ped.tmp
  
  ## Begin calculation ##
  kinship = matrix(0,N,N)
  kinship[1:N.old,1:N.old] = k.ped.old
  
  for (i in (N.old+1):N){
    
    i.M = ped[i,2]
    i.F = ped[i,3]
    
    for (j in i:1){

      if (j==i){
        ## if j's parents are in the pedigree ##
        if (i.M!=0 & i.F!=0){          
          kinship[i,j] = 1 + .5*kinship[i.M,i.F]
        }
      } else {
        if (i.M!=0 & i.F!=0){
          kinship[i,j] = .5*kinship[i.M,j] + .5*kinship[i.F,j]
          kinship[j,i] = kinship[i,j]
        }
      }

    }
  }
  return(kinship)
}




## Calculate empirical kinship ##
kinship.emp <- function(geno.v1 , geno.v2, map, method = "simple"){

  similar = 2 - abs(geno.v1 - geno.v2)
  no.comps = sum(is.na(similar)==FALSE)
  sum.sim = sum(similar , na.rm=TRUE)
  kin.emp = sum.sim / (2 * no.comps)
  return(kin.emp)
}







sex.warning <- function(){
  cat("There were problems found in the sex coding. Coding should be M/F or 1/0.\n")
}
