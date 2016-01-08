

###############################################
## compare genotypes at pairs of individuals
###############################################
setwd('/Users/mpinsky/Documents/Rutgers/Philippines/Genetics/genotyping/stacks_sensitivity_2015-07-17')
source('../../functions/readGenepop.R')
library(RCurl)

genfile = '../lax-rxstacks_2015-06-17/batch_1.r75m5.genepop' # original
genfile = 'laxnor-oldcat/r75m5.genepop' # small set of 23, no -r in process_radtags, original catalog
genfile = 'laxnor/batch_1.genepop' # small  set of 23, no -r in process_radtags, de novo catalog
genfile = 'default/batch_1.genepop' # small  set of 23, no -r in process_radtags, de novo catalog

dat = readGenepop(genfile)

ncol(dat)-2 # number of loci

	# set up individuals to compare
totest = vector('list', 5)
totest[[1]] = data.frame(ind1 = 'APCL13_048L1032', ind2 = c('APCL13_048L0881','APCL13_048L881'))
totest[[2]] = data.frame(ind1 = c('APCL13_060L0420', 'APCL13_060L420'), ind2 = c('APCL13_060L0421', 'APCL13_060L421'))
totest[[3]] = data.frame(ind1 = c('APCL13_119L0432', 'APCL13_119L432'), ind2 = c('APCL13_119L0433', 'APCL13_119L433'))
totest[[4]] = data.frame(ind1 = c('APCL13_226L0261', 'APCL13_226L261'), ind2 = 'APCL13_640L1231')
totest[[5]] = data.frame(ind1 = 'APCL13_120L1084', ind2 = 'APCL13_131L1092')

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


onematch = (alone1 == altwo1 & alone2 != altwo2) | (alone1 == altwo2 & alone2 != altwo1) | (alone2 == altwo1 & alone1 != altwo2) | (alone2 == altwo2 & alone1 != altwo1) # does one allele match but not the other?
homvhet = ((alone1 == altwo1 & alone2 != altwo2) | (alone1 == altwo2 & alone2 != altwo1) | (alone2 == altwo1 & alone1 != altwo2) | (alone2 == altwo2 & alone1 != altwo1)) & (alone1 == alone2 | altwo1 == altwo2) # a onematch where one of the genotypes is a homozygote (hom vs. het mismatch)
sum(onematch)
sum(homvhet) # the same, if all one allele matches are hom vs het mismatches
sum(!onematch)

rbind(genosone[which(!matches)], genostwo[which(!matches)]) # visually inspect
rbind(genosone[which(!matches)][onematch], genostwo[which(!matches)][onematch]) # visually inspect cases where they match on one allele
rbind(genosone[which(!matches)][!onematch], genostwo[which(!matches)][!onematch]) # visually inspect where no alleles match




