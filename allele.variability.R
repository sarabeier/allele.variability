####################################################
#
#TEMPORAL ALLELE VARIBILITY
#Last Version September 2019
#Sara Beier, Leibniz Institute for Baltic Sea Research
#
####################################################

rm(list=ls())
library(vegan)
library(car)
library(foreach)
library(doMC)

n.min <- 500 #minimum number of counts for each gene ortholog
n.boot <- 100 #number of permutations (bootstrapping)
n.thread <- 1 #number of threads
registerDoMC(n.thread) # define the default number of cpus used by R

print(paste("Minimum number of counts for each gene ortholog: n.min =",n.min))
print(paste("Number of permutations for bootstrapping: n.boot =",n.boot))
print(paste("Number of threads: n.thread =", n.thread))
cat("For demnostration the number of permutations was set to n.boot=100, while for publication n.boot=1000 was used.
Warning: high values for n.boot may require a long computing time!
We recommend to to adjust the number of threads (n.thread) for high values for n.boot" ,sep="\n")


#upload and format input data
cat("upload and format input data
countdata-file: tab-separated file with first column containing ortholog ids, second column containing allele.ids and all following columns containg the number of reads per allele.id in each metagenome:
example: counts.tab
info-file: containing information about the average per cell copy number, subcellular location and occurrence of orthologs in prokaryote genomes:)
example: TableS5.tab as published in Beier et al." ,sep="\n")
#counts <- read.table("/Users/sara/Documents/Manuscripts/LTG/countdata/counts.men.tab", header=T, sep ='\t') #upload count data as provided in supplmentary data x/y/z
counts <- read.table("counts.tab", header=T)
names(counts)[1:2] <- c('ortholog.id','allele.id')
rich.allele <- as.data.frame(table (counts$ortholog.id)) #total allel richness

info <- read.table("TableS5.tab", header=T, sep='\t')[,c(1,9,10,11)]#upload information about subcellular location, average gene copy nr and the frequency (occurence) of gene orthologs in prokaryotes as provided in Table S3
colnames (info) <- c("ortholog.id", "loc","copy", "occurence")
info <- info[info$occurence<0.67,] #keep auxiliary genes
info <- merge(info,rich.allele,by.x='ortholog.id', by.y='Var1')
info$meta <- info$Freq/info$copy #estimate metacommunity size for genes

counts.K <- aggregate(. ~ ortholog.id, data=counts[,c(1,3:dim(counts)[2])], FUN=sum) #KEGG counts
counts.K$min <- apply(counts.K[,2:dim(counts.K)[2]],1,min)
counts.Kmin <- counts.K[counts.K$min>n.min,] #selects ortholog.ids with n>n.min reads in each sample
counts.min <- counts[counts$ortholog.id %in% counts.Kmin$ortholog.id,] #select alleles that are annotated to ortholog.id>n.min counts
counts.Kmin$ortholog.id <- factor(counts.Kmin$ortholog.id)
minK <- min(counts.Kmin$min) #minimal number of reads found for the ortholog.id with lowest read number

############################################
#loop for bootstrapping
print("Loop for bootstrapping")
mylist <- foreach(j=1:n.boot) %dopar% {
  
  #loop to create betadispersion values for each individual gene ortholog
  datalist = list()
  for (i in 1:length(levels(counts.Kmin$ortholog.id))){
    #for (i in 1:5){
    ortholog.id.counts <-counts.min[counts.min$ortholog.id==levels(counts.Kmin$ortholog.id)[i],]
    row.names(ortholog.id.counts) <-ortholog.id.counts$allele.id
    sub <- t(rrarefy(t(ortholog.id.counts[,3:dim(ortholog.id.counts)[2]]),minK)) #subsampling
    dis <- vegdist(t(sub),method="bray")
    #dis <- vegdist(t(sub),method="jaccard", binary=T) #alternatve to estimate distance using the Jaccard index
    #dis <- vegdist(t(sub),method="horn") #alternatve to estimate distance using the Morisita-Horn index
    group <- rep(levels(counts.Kmin$ortholog.id)[i],dim(sub )[2]) #define groups to estimate betadispersion
    mod1 <- betadisper(dis, group) #estimate betadispersion
    d <- data.frame(mod1$distances, mod1$group)
    datalist[[i]] <- d
  }
  mod.ds <- do.call(rbind, datalist)
  colnames(mod.ds) <- c("allele.variability","ortholog.id")
  beta.j <- aggregate(. ~ ortholog.id, data=mod.ds, FUN=mean) #aggregate sample distance to centroid for each gene ortholog by mean
  beta.j <- merge (beta.j, info, by.x="ortholog.id",by.y="ortholog.id") #add information about subcellular location, average gene copy number and occurence to dataframe
  beta.sub <- beta.j[beta.j$loc %in% c('cytoplasmatic', 'non-cytoplasmatic'),] #only include cases with defined subcellular location for statistics
  
  #type III ANCOVA
  #randomize data for NULL model
  beta.rand <- beta.sub
  beta.rand$loc <- sample(beta.rand$loc)
  beta.rand$copy <- sample(beta.rand$copy)
  beta.rand$occurence <- sample(beta.rand$occurence)
  beta.rand$meta <- sample(beta.rand$meta)
  
  #NULL model (input: beta.rand; after shuffling the order of values in beta.sub)
  fit.null.iii <- lm(allele.variability ~ scale(log(copy),scale=F) + scale(log(meta),scale=F) + loc, data=beta.rand )
  NULLM <- Anova(fit.null.iii, type="III")
  
  #FULL model (input: beta.sub)
  fit.full.iii <- lm(allele.variability ~ scale(log(copy),scale=F) + scale(log(meta),scale=F) + loc, data=beta.sub )
  ANCOVA <- Anova(fit.full.iii, type="III")
  
  cop.j <- cbind.data.frame(ANCOVA[3:4])[2,]
  cop.j$test <- cbind.data.frame(ANCOVA[3])[2,]-cbind.data.frame(NULLM[3])[2,] #substract F-value from NULL model from FULL model F-value
  met.j <- cbind.data.frame(ANCOVA[3:4])[3,]
  met.j$test <- cbind.data.frame(ANCOVA[3])[3,]-cbind.data.frame(NULLM[3])[3,] #substract F-value from NULL model from FULL model F-value
  loc.j <- cbind.data.frame(ANCOVA[3:4])[4,]
  loc.j$test <- cbind.data.frame(ANCOVA[3])[4,]-cbind.data.frame(NULLM[3])[4,] #substract F-value from NULL model from FULL model F-value
  
  return(list(met.j,cop.j,loc.j,beta.j))
}

met <- do.call(rbind,sapply(mylist,`[`,1))
cop <- do.call(rbind,sapply(mylist,`[`,2))
loc <- do.call(rbind,sapply(mylist,`[`,3))


#bootstrapped statistics 
#p-value: fraction of F(full-model)-F(null-model) >0 after n=n.boot permutations
#F-value: mean F-value after 1000 permutations
print('bootstrapped statistics')

Fm <- mean(met[,1]) #F-value metacommunity size
pm <- length(met[,3][met[,3] < 0])/length(met[,3]) #bootstrapped p-value metacommunity size

Fc <- mean(cop[,1])#F-value copy number
pc <- length(cop[,3][cop[,3] < 0])/length(cop[,3]) #bootstrapped p-value copy number

Fl <- mean(loc[,1])#F-value subcellular location
pl <- length(loc[,3][loc[,3] < 0])/length(loc[,3]) #bootstrapped p-value subcellular location


print(paste("bootstrapped p-value for metacommunity size with", n.boot,"permutations:",ifelse(pm==0, paste("<",1/n.boot),pm)))
print(paste("average F-value for metacommunity size with", n.boot,"permutations:",Fm))
print(paste("bootstrapped p-value for gene copy number with", n.boot,"permutations:",ifelse(pc==0, paste("<",1/n.boot),pc)))
print(paste("average F-value for gene copy number with", n.boot,"permutations:",Fc))
print(paste("bootstrapped p-value for subcellular location with", n.boot,"permutations:",ifelse(pl==0, paste("<",1/n.boot),pl)))
print(paste("average F-value for subcellular location with", n.boot,"permutations:",Fl))

#create overview file with average allele variability for each gene ortholog
print(paste('create overview file with average allele variability for each gene ortholog obtained after',n.boot,'permutations'))
beta.all = do.call(rbind,sapply(mylist,`[`,4))
beta.all$loc <- as.character(beta.all$loc) #as.numeric
beta.all[is.na(beta.all)] <- "undefined" #replace NA by 'undefined'
beta <- aggregate(. ~ortholog.id+loc,beta.all,FUN=mean) #aggregate by KEGG-ID (and subcellular location): input file for Fig. 1

l1<-length(beta[,2][beta[,2] == 'cytoplasmatic'])
print(paste("number of genes identified as cytoplasmatic:",l1))
l2<-length(beta[,2][beta[,2] == 'non-cytoplasmatic'])
print(paste("number of genes identified as non-cytoplasmatic:",l2))

