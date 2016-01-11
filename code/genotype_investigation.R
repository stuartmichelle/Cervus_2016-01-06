

###############################################
## compare genotypes at pairs of individuals
###############################################
setwd('/Users/michelle/Documents/Rutgers/Philippines/Genetics/parentage/Cervus_2016-01-06/')
source('code/readGenepop_space.R')
library(RCurl)

# read in list of troublesome samples to analyze
comparisons <- ____

# genfile = '../lax-rxstacks_2015-06-17/batch_1.r75m5.genepop' # original
# genfile = 'laxnor-oldcat/r75m5.genepop' # small set of 23, no -r in process_radtags, original catalog
# genfile = 'laxnor/batch_1.genepop' # small  set of 23, no -r in process_radtags, de novo catalog
# genfile = 'default/batch_1.genepop' # small  set of 23, no -r in process_radtags, de novo catalog
genfile = 'DP20_edited_genepop/DP20g95_edited.genepop' #all of seq03-09 after dDocent filtering, no rad_haplotyper, and minus 3 non-APCL samples - "pop" must be lower case in order for readGenepop to work

dat = readGenepop(genfile)

ncol(dat)-2 # number of loci

# calculate number of pairwise comparisons
numpairwisecomparisons <- 

	# set up individuals to compare
totest = data.frame(ind1 = character(numpairwisecomparisons), ind2 = character(numpairwisecomparisons))
row <- 1
for(i in 1:nrow(comparisons)){
	if(!is.na(comparisons$Ligation_1) & !is.na(comparisons$Ligation_2)){
		totest$ind1[row] <- comparison$Ligation_1[i]
		totest$ind2[row] = comparisons$Ligation_2[i]
		row = row+1
	}
	
}


totest[[1]] = data.frame(ind1 = c('APCL_14126L1321','APCL_14126L1321'),  ind2 = c('APCL_14126L1701','APCL_14126L1701'))
totest[[2]] = data.frame(ind1 = c('APCL_13333L0357','APCL_13333L0357'),  ind2 = c('APCL_13333L1721','APCL_13333L1721'))
totest[[3]] = data.frame(ind1 = c('APCL_13565L0670','APCL_13565L0670'),  ind2 = c('APCL_13565L1735','APCL_13565L1735'))
totest[[4]] = data.frame(ind1 = c('APCL_14049L0726','APCL_14049L0726'),  ind2 = c('APCL_14049L1740','APCL_14049L1740'))
totest[[5]] = data.frame(ind1 = c('APCL_13584L0677','APCL_13584L0677'),  ind2 = c('APCL_13584L1736','APCL_13584L1736'))
totest[[6]] = data.frame(ind1 = c('APCL_14255L0737','APCL_14255L0737'),  ind2 = c('APCL_14255L1742','APCL_14255L1742'))
totest[[7]] = data.frame(ind1 = c('APCL_14452L0744','APCL_14452L0744'),  ind2 = c('APCL_14452L1754','APCL_14452L1754'))
totest[[8]] = data.frame(ind1 = c('APCL_14425L1343','APCL_14425L1343'),  ind2 = c('APCL_14425L1729','APCL_14425L1729'))
totest[[9]] = data.frame(ind1 = c('APCL_13651L0980','APCL_13651L0980'),  ind2 = c('APCL_13651L1725','APCL_13651L1725'))
totest[[10]] = data.frame(ind1 = c('APCL_13207L0444','APCL_13207L0444'),  ind2 = c('APCL_13207L1024','APCL_13207L1024'))
totest[[11]] = data.frame(ind1 = c('APCL_13087L0423','APCL_13087L0423'),  ind2 = c('APCL_13087L1719','APCL_13087L1719'))
totest[[12]] = data.frame(ind1 = c('APCL_13130L0432','APCL_13130L0432'),  ind2 = c('APCL_13130L1720','APCL_13130L1720'))
totest[[13]] = data.frame(ind1 = c('APCL_13587L0678','APCL_13587L0678'),  ind2 = c('APCL_13587L1737','APCL_13587L1737'))
totest[[14]] = data.frame(ind1 = c('APCL_14310L0597','APCL_14310L0597'),  ind2 = c('APCL_14310L1748','APCL_14310L1748'))
totest[[15]] = data.frame(ind1 = c('APCL_14424L1342','APCL_14424L1342'),  ind2 = c('APCL_14424L1747','APCL_14424L1747'))
totest[[16]] = data.frame(ind1 = c('APCL_14006L0722','APCL_14006L0722'),  ind2 = c('APCL_14006L1738','APCL_14006L1738'))
totest[[17]] = data.frame(ind1 = c('APCL_14032L0947','APCL_14032L0947'),  ind2 = c('APCL_14032L1739','APCL_14032L1739'))
totest[[18]] = data.frame(ind1 = c('APCL_14423L1341','APCL_14423L1341'),  ind2 = c('APCL_14423L1746','APCL_14423L1746'))
totest[[19]] = data.frame(ind1 = c('APCL_13082L1060','APCL_13082L1060'),  ind2 = c('APCL_13082L1688','APCL_13082L1688'))
totest[[20]] = data.frame(ind1 = c('APCL_13057L1040','APCL_13057L1040'),  ind2 = c('APCL_13057L1718','APCL_13057L1718'))
totest[[21]] = data.frame(ind1 = c('APCL_14557L1350','APCL_14557L1350'),  ind2 = c('APCL_14557L1717','APCL_14557L1717'))

	# to hold the results
a = rep(NA, length(totest))
out = data.frame(indivs = a, matches=a, mismatches=a, perc=a, hetmatch=a, hetmism=a, perchet=a)

for(i in 1:length(totest)){
	datrow = which(as.character(dat$names) %in% c(as.character(totest[[i]]$ind1), as.character(totest[[i]]$ind2)))
	out$indivs[i] = paste(dat$names[datrow], collapse = ', ')

	genosone = dat[datrow[1], 3:ncol(dat)]
	genostwo = dat[datrow[2], 3:ncol(dat)]
	matches = genosone == genostwo # where the two genotypes match or not
	matches[genosone == '0000' | genostwo == '0000'] = NA # remove missing data from calculations

	out$matches[i] = sum(matches, na.rm=TRUE) # number of matching loci
	out$mismatches[i] = sum(!matches, na.rm=TRUE) # number of mismatching loci
	out$perc[i] = 100*signif(sum(!matches, na.rm=TRUE)/(sum(matches, na.rm=TRUE) + sum(!matches, na.rm=TRUE)),2) # proportion mismatching


	alone1 = substr(genosone, 1,2) # first allele in individual one
	alone2 = substr(genosone, 3,4) # second allele in individual one
	altwo1 = substr(genostwo, 1,2) # first allele in individual one
	altwo2 = substr(genostwo, 3,4) # second allele in individual one

	hets = (alone1 != alone2) | (altwo1 != altwo2)
	hets[alone1 == '00' | alone2 == '00' | altwo1 == '00' | altwo2 == '00'] = NA

	out$hetmatch[i] = sum(hets & matches, na.rm=TRUE) # number of matching heterozygote loci
	out$hetmism[i] = sum(hets & !matches, na.rm=TRUE) # number of mismatching loci where at least one indiv is het
	out$perchet[i] = 100*signif(sum(hets & !matches, na.rm=TRUE)/(sum(hets & !matches, na.rm=TRUE)+sum(hets & matches, na.rm=TRUE)),2)

}

out

write.csv(out, file="regenotyped_mismatch.csv")

onematch = (alone1 == altwo1 & alone2 != altwo2) | (alone1 == altwo2 & alone2 != altwo1) | (alone2 == altwo1 & alone1 != altwo2) | (alone2 == altwo2 & alone1 != altwo1) # does one allele match but not the other?
homvhet = ((alone1 == altwo1 & alone2 != altwo2) | (alone1 == altwo2 & alone2 != altwo1) | (alone2 == altwo1 & alone1 != altwo2) | (alone2 == altwo2 & alone1 != altwo1)) & (alone1 == alone2 | altwo1 == altwo2) # a onematch where one of the genotypes is a homozygote (hom vs. het mismatch)
sum(onematch)
sum(homvhet) # the same, if all one allele matches are hom vs het mismatches
sum(!onematch)

rbind(genosone[which(!matches)], genostwo[which(!matches)]) # visually inspect
rbind(genosone[which(!matches)][onematch], genostwo[which(!matches)][onematch]) # visually inspect cases where they match on one allele
rbind(genosone[which(!matches)][!onematch], genostwo[which(!matches)][!onematch]) # visually inspect where no alleles match




