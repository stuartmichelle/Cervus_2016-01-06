##############
## prep

setwd('/Users/macair/Documents/Philippines/Genetics/parentage/Cervus_2016-01-06/')

dat = read.csv('DP20g95/DP20g95_ID.csv', stringsAsFactors=FALSE)
nrow(dat) # 65

	# add year of the sample - skipped this because dDocent IDs are different and I'm not sure what this code is doing
dat$First.year = as.numeric(paste('20', gsub('APCL_|.*', '', dat$First.ID), sep=''))
dat$Second.year = as.numeric(paste('20', gsub('APCL|_.*', '', dat$Second.ID), sep=''))

	# fix ID to always have 4-digit ligation IDs (needed for matching against Sample_Data google sheet) - skipped - ID's already have 4-digit ligation IDs
ind = grep('L[[:digit:]]{3}$', dat$First.ID) # rows with 3-digit ligation IDs
dat$First.ID[ind] = gsub('L([[:digit:]]{3})$', 'L0\\1', dat$First.ID[ind])
ind = grep('L[[:digit:]]{3}$', dat$Second.ID)
dat$Second.ID[ind] = gsub('L([[:digit:]]{3})$', 'L0\\1', dat$Second.ID[ind])

	# add sampleid
dat$First.SampleID = gsub('L[[:digit:]]{1,}$', '', dat$First.ID)
dat$Second.SampleID = gsub('L[[:digit:]]{1,}$', '', dat$Second.ID)

	# add lat/lon from our Google Sheet
require(googlesheets)
# gs_auth(new_user = TRUE) # run this if having authorization problems
mykey = '1Rf_dFJ5WK-vTTsIT_kHHOcFrKzQtMFtKiuXiFw1lh9Y' # for Sample_Data sheet
gssampdat <- gs_key(mykey)
sampledata <- gs_read(gssampdat, ws='Samples')

m1 = sampledata
names(m1) = paste('First.', names(m1), sep='')
dat = merge(dat, m1, by.x='First.SampleID', by.y = 'First.Sample_ID', all.x=TRUE)
m2 = sampledata
names(m2) = paste('Second.', names(m2), sep='')
dat = merge(dat, m2, by.x='Second.SampleID', by.y = 'Second.Sample_ID', all.x=TRUE)
	
	# distance between samples
require(fields)
# source('greatcircle_funcs.R') # alternative, probably faster
alldists = rdist.earth(as.matrix(dat[,c('First.Lon', 'First.Lat')]), as.matrix(dat[,c('Second.Lon', 'Second.Lat')]), miles=FALSE, R=6371) # see http://www.r-bloggers.com/great-circle-distance-calculations-in-r/ # slow because it does ALL pairwise distances, instead of just in order
dat$distkm = diag(alldists)

	# add ligation ID
dat$First.LigID = as.numeric(gsub('APCL[[:digit:]]{2}_[[:digit:]]{1,}L', '', dat$First.ID))
dat$Second.LigID = as.numeric(gsub('APCL[[:digit:]]{2}_[[:digit:]]{1,}L', '', dat$Second.ID))
	
	# add ligation pool
ligdata <- gs_read(gssampdat, ws='Ligations')

l1 = ligdata
names(l1) = paste('First.', names(l1), sep='')
dat = merge(dat, l1[,c('First.Ligation_ID', 'First.Pool')], by.x='First.ID', by.y = 'First.Ligation_ID', all.x=TRUE)

l2 = ligdata
names(l2) = paste('Second.', names(l2), sep='')
dat = merge(dat, l2[,c('Second.Ligation_ID', 'Second.Pool')], by.x='Second.ID', by.y = 'Second.Ligation_ID', all.x=TRUE)

	# add mismatch rate
	dat$Mismatch.prop = dat$Mismatching.loci/dat$Matching.loci
	
	# order by number of matching loci
dat = dat[order(dat$Matching.loci),]

#####################
## investigation
hist(dat$Matching.loci, col='grey', breaks=50) # most match at <1000 SNPs

hist(dat$Mismatch.prop, col='grey', breaks=50) # interesting second hump < 0.02

plot(dat$Matching.loci, dat$Mismatching.loci)
abline(0, 0.005) # 5% error rate

plot(dat$Matching.loci, dat$Mismatch.prop)
abline(h=0.015)

colnms = c('First.ID', 'Second.ID', 'Loci.typed', 'Loci.typed.1', 'Matching.loci', 'Mismatching.loci', 'First.Size', 'Second.Size', 'First.Tail_color', 'Second.Tail_color', 'distkm', 'First.Notes', 'Second.Notes')

	# exact matches
sum(dat$Status == 'Exact match') # 0 exact matches

	# high quality fuzzy matches
ind = dat$Matching.loci>=2000 # seem like good matches. out of ~11k loci
sum(ind) # only 35
dat[ind,colnms]
table(dat$First.year[ind], dat$Second.year[ind]) # figure out the years of the matches: most are 2013->2014

plot(dat$First.Lat[ind], dat$Second.Lat[ind]) # most are close
plot(dat$First.Lon[ind], dat$Second.Lon[ind])

hist(dat$distkm[ind], breaks=40, col='grey') # 3 outliers > 10km
plot(dat$Matching.loci[ind], dat$distkm[ind], log='y'); abline(h=5) # 3 outliers

dat[ind & (dat$distkm > 1 & !is.na(dat$distkm)),] # fairly distant matches

dat[ind & (dat$distkm > 5 & !is.na(dat$distkm)),colnms] # really distant matches. seem like false positives.

suspect = c('APCL13_584', 'APCL14_155', 'APCL13_130', 'APCL14_032') # see notes. unlikely based on distance and sizes

ind = dat$Matching.loci>=2000 & !(dat$First.SampleID %in% suspect) & !(dat$Second.SampleID %in% suspect) # seem like good matches, excluding the impossible/suspect ones based on distance
	sum(ind) # 33

dat[ind, colnms]
table(dat$First.year[ind], dat$Second.year[ind]) # figure out the years of the matches

plot(dat$First.Size[ind], dat$Second.Size[ind]); abline(0,1) # most of the time, second size is bigger, but not always

dat[ind & (dat$First.Size > dat$Second.Size) & dat$First.year == 2013 & dat$Second.year==2014, colnms] # shrunk = unlikely, except the 9.4 to 9.3 one. works because samples ordered by name/year (2013 before 2014)

dat[ind & (abs(dat$First.Size - dat$Second.Size)>1 & dat$First.year == dat$Second.year), colnms] # different sizes in same year = unlikely


suspect = c('APCL13_584', 'APCL14_155', 'APCL13_130', 'APCL14_032', 'APCL13_402', 'APCL14_452', 'APCL13_587', 'APCL14_255', 'APCL13_565', 'APCL14_049', 'APCL13_087', 'APCL14_006', 'APCL13_615', 'APCL14_392', 'APCL13_333', 'APCL13_334', 'APCL14_310', 'APCL14_318', 'APCL14_309L1018 APCL14_317', 'APCL13_639', 'APCL13_651') # see notes. unlikely bases on sizes, plus ones excluded before
	indbad = dat$Matching.loci>=2000 & (dat$First.SampleID %in% suspect | dat$Second.SampleID %in% suspect)
		sum(indbad) # 10

ind = dat$Matching.loci>=2000 & !(dat$First.SampleID %in% suspect) & !(dat$Second.SampleID %in% suspect) # seem like good matches, excluding the impossible/suspect ones based on distance
	sum(ind) # 25, out of 35 (71%) that looked OK just based on # matches
	
dat[ind, colnms]

boxplot(dat$Matching.loci[indbad], dat$Matching.loci[ind], names=c('suspect?', 'good?')) # "good" matches do not appear to have more matching loci

boxplot(abs(dat$First.LigID[indbad]-dat$Second.LigID[indbad]), abs(dat$First.LigID[ind]-dat$Second.LigID[ind]), names=c('suspect', 'good?'), ylab='Difference in ligation ID')

	t.test(abs(dat$First.LigID[indbad]-dat$Second.LigID[indbad]), abs(dat$First.LigID[ind]-dat$Second.LigID[ind]))
	wilcox.test(abs(dat$First.LigID[indbad]-dat$Second.LigID[indbad]), abs(dat$First.LigID[ind]-dat$Second.LigID[ind]))

	par(mfrow=c(2,1)) # histogram of ligation ID differences
	hist(abs(dat$First.LigID[indbad]-dat$Second.LigID[indbad]), col='grey', xlim=c(0,1000), breaks=seq(0,1000,length.out=40), main='Suspect', xlab='Difference in ligation ID', ylim=c(0,3))
	hist(abs(dat$First.LigID[ind]-dat$Second.LigID[ind]), col='grey', xlim=c(0,1000), breaks=seq(0,1000,length.out=40), main='Good?', xlab='Difference in ligation ID', ylim=c(0,3))
	
dat[indbad & abs(dat$First.LigID-dat$Second.LigID) < 100, colnms]

dat[indbad, c('First.ID', 'Second.ID', 'First.Pool', 'Second.Pool', 'distkm')]


