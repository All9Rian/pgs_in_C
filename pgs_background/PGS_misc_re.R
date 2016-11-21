loadData<-function(inputFilePath){
# @inputFilePath: The folder where the new batch is
# @return libData, A list which length equals to the number of new seq libdata folder
	'"
	:type inputFilePath: string
	:rtype: list(data.frame)
	"'
	# the structure will be kept, with lib name as folder
    # in which there is a sample_100k.stat.gz file providing a 6 columns window based summary
    libList <- dir(inputFilePath)
    libList = libList[ !libList %in% c('rmdup_job_ids', 'QC_CNV.report')]
    libStat <- lapply(libList, function(lid) {
        ssf<-dir(file.path( inputFilePath, lid), pattern = 'sample_100k.stat')
        if(length(ssf)!=1) {
            warning('Sample_100k.stat.(gz) file doesnot exist or not unique ', lid)
            return(NULL)
        }
		#print( file.path( inputFilePath, lid, ssf[1] ))
        dat <- read.table( file.path( inputFilePath, lid, ssf[1] ))[,1:6] # all chr included
        # the dat has a fixed structure currently, so shall be named accordinglly
        if(ncol(dat) == 6) {
            colnames(dat) <- c('chr', 'win_id', 'reads', 'u_reads', 'map_rate', 'gc')
        } else { 
            warning('The input data table is not supported!', lid) # due to column no mismatch
            return(NULL)
        }
        chrs <- paste0('chr', c(1:22, 'X', 'Y'))
        dat[, 'chr'] <- paste0('chr', gsub('CHR', '', toupper(dat[, 'chr'])))
        dat[, 'chr'] <- factor(dat[, 'chr'], chrs, chrs, ordered=TRUE)
        return(dat)
    })
    names(libStat)<-libList
    
    #load QC stat
    # qc stat in a df
    # the qc stuff sits in the parent folder of the new batch
    # the only overhead might be that the qc contains all the info for the whole batch, some extra readin
    qcLabel<-c("Sample_ID", "Reads_Mb", "GC_Content", "Q30_Ratio", "Align_Ratio", "UR_Ratio", "Duplication_Ratio", "Coverage", "UR_Upper150_Ratio")
    qcRepFile<-paste0(gsub('(.+)analysis', '\\1', inputFilePath), 'QC_CNV.report')
    
    qcTab<-data.frame(do.call(rbind, lapply(libList, function(x) rep(NA, length(qcLabel)))))
	#print(qcTab)
    colnames(qcTab)<-qcLabel
    rownames(qcTab)<-qcTab$Sample_ID<-libList
    
    if(file.exists(qcRepFile)) {
        qcTmp<-read.table(qcRepFile, stringsAsFactors = F) # for historical reasons, the first column can be longer, but the first item before '_' is assumed to be unique
        if(ncol(qcTmp)==9) {
            colnames(qcTmp)<-qcLabel
            qcTab[qcTmp$Sample_ID, ]<-qcTmp
        } else if (ncol(qcTmp)==8) {
            qcTab[qcTmp$Sample_ID, qcLabel[-6]]<-qcTmp
        } else {
            warning('QC_CNV.report formate mismatch, skipped & NA set!', inputFilePath)
        }
    } else {
        warning('NO QC_CNV.report found in the parent folder, NA set!', inputFilePath)
    }
    libData<-list(winStat=libStat, qcStat=qcTab)
    return(libData)
}

updateDB <- function(dbPath, gather = 40, useVal = "residD", depth = 11438000, family = "poisson") {
# update background data base in the same process method as every new incoming samples
# @dbPath: The path of background data base
# @gather: Parameter in function: feedBack
# @useVal: Name of used column in sample_100K
# @depth: uniform depth
# @return: NULL
	'"
	:type dbPath: string
	:type gather: int
	:type useVal: string
	:type depth: int
	:type family: string
	:rtype: void
	"'
    sampleList = loadData( dbPath )$winStat
	db = lapply(sampleList, function(dat) {
		dat = tagDat(dat, depth)
		dat = eliminateDiff( dat, family="poisson")
		dat = feedBack(dat, gather = 40, useVal="residD")
	})
	save( db, file = file.path(dirname( dbPath ), "db.Rda") )
}

tagDat <- function(dat, depth = 10^7, mpfilt.p = c(rep(0.1, 22), 0.2, 0.3), mpfilt.cut = 0.65) {
# @dat: sample100k data
# @depth: average depth of unique_reads
# @mpfilt.p: quantile in every chromosome
# @mpfilt.cut: mapbility filter cut at
# @return dat: a data.frame contain a new column 'tag' stand for the windows filtered or not
	'"
	:type dat: data.frame
	:type depth: int
	:type mpfilt.p: float
	:type mpfilt.cut: float
	:rtype: data.frame
	"'
    index = cumsum( table( dat$chr))
	dat$u_reads = round( dat$u_reads / sum(dat$u_reads) * depth)
	dat = do.call( rbind, mapply( function( dati, mpfilt) {
        mpbility = dati[,"u_reads"] / (dati[,"map_rate"] + 1)
        dati$tag = mpbility >= min(quantile(mpbility, mpfilt, na.rm = T), mpfilt.cut)
        if( unique(dati$chr) == "chrY") {
			set.seed(16)
			section = c(145:171, 211:249, 254:294, 350:372, 379:400)
			Val = dati[ ,"u_reads"][ section]
			if( median(Val) >= 10) {
				thre = quantile(Val, c(0, 0.8))
				Val = Val[Val > thre[1] & Val < thre[2]]
			}	
			ind = round( runif( length(dati[ ,"u_reads"][ -section]), min = 1, max = length( Val)))
	        dati[ ,"u_reads"][ -section] = Val[ ind]
            dati$tag = T
        }
        list(dati)
    }, split(dat, dat$chr), mpfilt.p) )
	return(dat)
}

eliminateDiff <- function(dat, family="poisson") {
# Eliminate the difference in every data
# @dat: batch sample100k data
# @family: glm method parameter
# @return Dat: a data.frame contains some new columns 'residS', 'residT', 'residT'
	'"
	:type dat: data.frame
	:type family: string
	:rtype: data.frame
	"'
	Dat = dat
	dat = Dat[Dat$tag, ]
	dat$residS = do.call(c, lapply(split(dat, dat$chr), function(dati) {
        if( family == "poisson") {
            modS = glm(formula = (dati[,'u_reads'] - min(dati[,'u_reads'])) ~ dati$gc , family=poisson(link = "log"))
            return( modS$residuals / sd(modS$residuals) * 0.1)
        }else {
            modS = glm(formula = (dati[,'u_reads'] - min(dati[,'u_reads'])) ~ dati$gc , family=gaussian(link = "identity"))
            return( modS$residuals / sd(modS$residuals) * 50)
        }
    }) )
    mod = glm( formula = (dat[,'u_reads'] - min(dat[,'u_reads'])) ~ dat$gc, family=poisson(link = "log"))
    dat$residT = mod$residuals
    dat = do.call(rbind, lapply(split(dat, dat$chr), function(dati) {
        if( family == "poisson") {
            dati$residD = ((dati$residS + median(dati$residT)) + 1) * 2
        }else dati$residD = (dati$residS + median(dati$residT)) / 100 + 2
        if( unique(dati$chr) == "chrY" & median(dati$u_reads) < 10) {
            dati$residD = dati$residD - median(dati$residD)
        }
        return(dati)
    }) )
    Dat$residD = NA
    Dat$residD[ Dat$tag] = dat$residD
    return(Dat)
}


feedBack <- function(dat, gather = 30, useVal = "residD") {
# Correct low-quality points to an acceptable simulated value, set.seed = 16
# @dat: a dataframe, load from sample100k
# @gather: the length of gather bin
# @useVal: colname of used value in "dat"
# @return dataframe, sample100k data which the useVal column is corrected
	'"
	:type dat: data.frame
	:type gather: int
	:type useVal: string
	"'
    set.seed( 16)
    dat = do.call( rbind, lapply( split( dat, dat$chr), function( dati) {#message(unique(dati$chr))
        if( sum( dati$tag) < 1) return( dati)
        if( unique( dati$chr) == "chrY") return(dati)
        highInd = which( dati$tag)
        mat = matrix( rep( which( !dati$tag), sum( dati$tag)), ncol = sum(dati$tag))
        tagMat = matrix( rep( highInd, sum( !dati$tag)), ncol = sum(dati$tag), byrow = T)
        distMat = mat - tagMat
        dati[,useVal][ !dati$tag] = apply( distMat, 1, function(x) {
            median( dati[ ,useVal][ highInd[ order( order( abs(x))) <= gather]])# + rnorm(1, 0, highQsd)
        })
        return(dati)
    }) )
	return(dat)
}

dbCorrention <- function(dat, db, useVal="residD") {
# Using background data base to correcte system error
# @dat: sample_100K data 
# @dbMat: background data base
# @useVal: colname of used value in "dat"
# @return dataframe, batch data after correction
	'"
	:type dat: data.frame
	:type dbMat: matrix
	:type useVal: string
	:rtype: data.frame
	"'
	dbMat = as.data.frame(do.call(cbind, lapply(db, function(dat) dat[,"residD"] )))
	index = cumsum( table( dat$chr ))
	x = paste0(do.call(c, lapply( 1:ncol(dbMat), function(i) paste0("dbVal[,", i, "]"))), collapse = "+")
	dat$result = do.call(c, lapply(list(1:index[22], (index[22] +1):index[23], (index[23] + 1):index[24]), function(ind) {
           	dbVal <- as.data.frame(dbMat[ ind,])
           	lm(eval(parse(text = paste0( "dat[,'residD'][ind] ~ ", x))))$residuals + median(dat[,"residD"][ind])
    }) )
	return(dat)
}

saveData <- function(dat, libn, corrDB, outputFilePath, useVal = "result") {
# Saving data including: "Bxplot.pdf"
#                        "Hist_of_Residuals.pdf"
#                        "result.csv"
#						 "Ctrl_info.csv"
#                        "Median.pdf"
# @dat: sample_100K data after data base correction
# @db: background data base
# @pathout: saving path
# @useVal: colname of used value in "dat"
# @return void
	'"
	
	"'
	dir.create( file.path( outputFilePath, libn))
	# plot 'boxplot' in "Boxplot.pdf"
	pdf( file.path( outputFilePath, libn, "Boxplot.pdf"), width = 10, height = 4)
	boxplot( split( dat[,useVal], dat$chr), main=libn, las=2)
	lapply( 1:24, function(i) {
           x = median( split( dat[,useVal], dat$chr)[[i]])
           text( i, x, round(x, digits=2), col=2)
    })
    dev.off()
	# plot 'hist' in "Hist_of_Residuals.pdf"
	pdf( file.path( outputFilePath, libn, "Hist_of_Residuals.pdf"), width = 10, height = 4)
	index = cumsum( table( dat$chr))
	layout( matrix(1:3, 1))
    showHist( dat[,useVal][ 1:index[22] ], main="1:22")
    showHist( dat[,useVal][ (index[22] + 1):index[23] ], main="23")
    showHist( dat[,useVal][ (index[23] + 1):index[24] ], main="24")
    dev.off()
	# write 'dat' into "result.csv" document.
	write.csv(dat, file=file.path(outputFilePath, libn, "result.csv"))
    # write 'Ctrl_info' into "Ctrl_info.csv") document.
	n = nrow(dat)
    x = dat[, useVal]
    Ctrl_info = data.frame("ID" = libn,
                           "p value to Norm" = chisq.test(x - min(x))$p.value,
                           "sd" = sd(x),
                           "Skewness" = n / ((n - 1) * (n - 2) * sd(x)^3) * sum( (x - mean(x))^3 ))
	write.csv(Ctrl_info, file=file.path(outputFilePath, libn, "Ctrl_info.csv"))
    # DB corresbond
	med = sapply( split( dat[ ,useVal], dat$chr), median)
	# plot 'median' in "Median.pdf"
	pdf( file.path( outputFilePath, libn, "Median.pdf"), width = 10, height = 4)
	plot( med, type="l", col = 3, main = libn)
	lapply(corrDB, function(dat) points(sapply( split(dat[ ,useVal], dat$chr), median)))
	dev.off()
}

# Plot histograme of modeling results(usually residuals of statistic model)
# @param dat, a dataframe, load from sample100k
# @param main, the main title of histograme
# @return null
showHist <- function(dat, main) {
    # hist(dat, breaks=1000, main=libn)
    n = length(dat)
    skew = n / ( (n - 1) * (n - 2) * sd(dat)^3) * sum( (dat - mean(dat))^3 )
    hist( dat, breaks = 1000, main=main )
    legend("top", legend = c('Skewness=', round(skew, digits=4)), bty='n')
}

aneuNormRatio2Karyo <- function(dat, pids=2:3, cps=0:4, dist='norm', useVar='u_reads', epsilon=0.001) {
# using a normalized coverage to represent the copy state of the whole chr, so the resulting ratio is not subject to common factors like chr length or gc content,
# @x: a list of sample_100k dataframe to represent all samples
# @pids: possible ploidy number to be evaluated
# @cps: possible copy number states to be considered
# @dist: the distribution used to evaluate the matching of copy states
# @return cns: a dataframe containing chr-wise ploidy assignement and copy states, the matching probablity and the original ratio
    "
    :type x: list(data.frame)
    :rtype: list(data.frame)
    "
	v <- tapply(dat[ ,useVar], dat[ ,'chr'], function(sdt) { median(sdt) })
    #v[is.na(v)]<-0 # force 0 when sometime the whole y is filtered out
    perVarRatio<-data.frame( names( v ), v / sum( v ), stringsAsFactors=F) #FIXME make another kind of visualization
    colnames(perVarRatio) <- c('chr', 'ratio')
    rownames(perVarRatio) <- names( v )
    # v2 more general
    cnp<-do.call(cbind, lapply(pids, function(pd) {
        if(dist == 't') {
            sapply(cps, function(cn) dt(perVarRatio$ratio/(1/23/pd), ncp=cn, df=22))
        } else{
            sapply(cps, function(cn) dnorm(perVarRatio$ratio/(1/23/pd), mean=cn, sd=1))
        }
    }))
    colnames(cnp)<-paste0(rep(paste0('ploid', pids), each=length(cps)), rep(paste0('cn', cps), length(pids)))
    bestIdx<-apply(cnp, 1, which.max)
    cnv<- bestIdx %% length(cps)
    cnv[cnv==0]<-length(cps)
    cnv<-cps[cnv]
    ploidy<-pids[ceiling(bestIdx / length(cps))]
    ratio = data.frame(ploidy=ploidy, CN=cnv, DP=apply(cnp, 1, max), Ratio=perVarRatio$ratio, stringsAsFactors = F)
	#print(ratio)
	Karyotype = chrCN2Karyo(ratio, epsilon=0.001)
	return(Karyotype)
}

chrCN2Karyo<-function(cndp, epsilon=0.001){
# this function take the output from aneuNormRatio and do the adjustment
# a, take the minimum of ploidy and cn if it is not diploid
# b, adjust the autosome by the the ploidy of sex chr if they are equal and more than 2
#@param cndp, a 24 rows x 4 columns matrix / df, to represent the predicted Karotype status
#@return allAneu, a string to represent the karyotype, e.g. 46,XX
	"	
	"
    chrs<-c(1:22, 'X', 'Y')
    if(cndp[24, 'CN']==0 & cndp[23, 'CN'] > 0 ) cndp[24, 'ploidy'] <- cndp[23, 'ploidy']
    sexPloidy<-(cndp[23:24, 'ploidy'])
    #a background DB check to adjust funny auto chrs, e.g. 17 19 22
    chrRatioAll<-NULL
    if(is.null(chrRatioAll)) data('cnvRatioAll', package='AnnoroadPD', envir = environment())
    inRange<-sapply(1:22, function(x) cndp[x, 'Ratio'] - range(chrRatioAll[[x]]))
	#print(cndp)
	#print(inRange)
    inRange[1,]<-inRange[1,]+epsilon
    inRange[2,]<-inRange[2,]-epsilon
    inRangeIdx<-which(sapply(1:22, function(i) {
        (inRange[1, i]+epsilon) >= 0 & (inRange[2, i]-epsilon) <= 0 & (!(cndp[i, 'ploidy'] %in% unique(sexPloidy)) |(i==19))
    }))
    cndp[inRangeIdx,'ploidy']<-cndp[inRangeIdx,'CN']<-2
    copyRate<-cndp[, 'CN'] / cndp[, 'ploidy'] # not yet used, intended for moasic
    # an artifical pushing to the normal state
    if(length(unique(sexPloidy)) == 1){
        copyCall<-apply(cndp[, c('ploidy', 'CN')], 1, function(x) ifelse(x[1]==unique(sexPloidy) | x[1]==2, x[2],  min(x)))
    }else{
        copyCall<-apply(cndp[, c('ploidy', 'CN')], 1, function(x) ifelse(x[1]==2, x[2], min(x)))
    }
    # consider tri or tera when both sex chr are in ploidy 3 while there is y present,
    if( (length(unique(sexPloidy)) == 1 & all(sexPloidy > 2) ) ){ #& (cndp[24, 'CN']>0) #y present is checked already
        adjIdx<-cndp[, 'ploidy'] != unique(sexPloidy)
        copyCall[adjIdx] <- copyCall[adjIdx] + unique(sexPloidy) - 2
        allAneu<-paste(ifelse(copyCall>unique(sexPloidy), paste0('+', chrs), ifelse(copyCall<unique(sexPloidy), paste0('-', chrs), ''))[1:22])
    } else {
        allAneu<-paste(ifelse(copyCall>2, paste0('+', chrs), ifelse(copyCall<2, paste0('-', chrs), ''))[1:22])
    }
    allAneu<-paste(allAneu[allAneu!=''], collapse='')
    totalChr<-sum(copyCall)
    # set NULL is nothing found
    if(allAneu=='') {
        allAneu<-NULL
    }
    allAneu<-paste(c(totalChr, paste(rep(chrs[23:24], copyCall[23:24]), collapse=''), allAneu), collapse=',')
    return(allAneu)
}

chrPValueEstimate <- function( dat, libn, corrDB, outputFilePath, fractile=0.01, useVal = "result") {
    resMat = do.call(rbind, lapply( corrDB, function(dat) dat$result ) )
    range = as.data.frame( t( apply(resMat, 2, function(v) c(quantile(v, 1-fractile), quantile(v, fractile)))))
    # copy number probbility estimate
	range$chr = dat$chr
    range$result = dat$result
    cn.p = sapply(split(range, range$chr), function(rangei) {
        sum( rangei$result > rangei$`1%` & rangei$result < rangei$`99%`) / dim(rangei)[1]
    })
    write.csv( cn.p, file = file.path( outputFilePath, libn, "copyNumber_p.csv"))
	# FIXME give the result whether more or less
	return(cn.p)
}

cnvSearch <- function( dat, libn, minWin=20, minRun=10, maxGap=3, useVal="result", outputFilePath) {
	libPathout = file.path( outputFilePath, libn)
	dir.create( file.path( libPathout, "Graphs"))
	
	M = median(dat[, useVal])
	cnv = lapply( split( dat, dat[,'chr']), function( dati) {
        # calculate median
        m = median(dati[, useVal])
        if (unique(dati$chr) != "chrY" & unique(dati$chr) != "chrX") {
            if ( m > (M-0.3) & m < (M+0.3) | m > (M/2-0.3) & m < (M/2+0.3)) { med = m } else { med = M }
        } else { med = m }
        v = dati[ ,useVal] - med
        SegTable = HaarSeg::haarSeg( v, breaksFdrQ = 1e-30, haarStartLevel = 15*sd(v), haarEndLevel = qnorm(0.95, sd = sd(v)) * 30)$SegmentsTable
        haarStartLevel = quantile(apply(SegTable, 1, function(segRow) sd( dati[,useVal][segRow[1]:(segRow[1] + segRow[2] - 1)])), 0.95, na.rm = T)
        segRes = HaarSeg::haarSeg( v, breaksFdrQ = 1e-30, haarStartLevel = 15*haarStartLevel, haarEndLevel = qnorm(0.95, sd = sd(v)) * 30)
        segSD = sd(segRes$SegmentsTable[,3])
        dati$seg = segRes$Segmented
        #temp = apply(segRes$SegmentsTable, 1, function(segRow) rep( segSD + sd(dati[,useVal][segRow[1]:(segRow[1]+segRow[2]-1)]), segRow[2]))
		temp = apply(segRes$SegmentsTable, 1, function(segRow) rep( 0.5*segSD + 1.5*sd(dati[,useVal][segRow[1]:(segRow[1]+segRow[2]-1)]), segRow[2]))
        if (!is.list(temp)) dati$segCut = temp else dati$segCut = do.call(c, temp )
        # save data: Results/CNV_report.csv
        posi = cnvRevise( dati, med, minWin, minRun, maxGap)
        # sample CNV plot: Graphs/chr( i )plot.pdf
        chrPlot( dati, posi, med, haarStartLevel, useVal="result", pathout=file.path( libPathout, "Graphs") )
        return( posi)
    })

	cnv = cnv[ lapply(cnv, length) > 0]
	cnvDat = do.call( rbind, lapply( cnv, function( chri_cnv) do.call( rbind, chri_cnv)))
	write.csv( cnvDat, file.path( libPathout, 'CNV_report.csv' ), quote = FALSE, row.names = FALSE)
	data('cnvDB', package='AnnoroadPD', envir = environment())
	CNV = cnvstat(cnv, libn)
    return(CNV)
}

cnvRevise <- function( dati, med, minWin=8, minRun=3, maxGap=3, p.val = 0.6) {
    # calculate SegCut
    dati[, 'segCut'][dati[, 'segCut'] < 0.1] = 0.1
    win = which( abs( dati[, 'seg']) > dati[, 'segCut'])
    seg = dati[ ,'seg'][ win]
    if( length( win) < minRun ) return( list())
    cnvList = list(); ind = vector(); j = 1
    for( i in 1: ( length(win) - 1)) {
        if( win[i + 1] - win[i] <= maxGap & seg[i + 1] * seg[i] > 0 ) {
            ind = unique( c( ind, win[i], win[i + 1]))
            if( length(ind) >= minRun & ind[length(ind)] - ind[1] >= minWin) {
				cnvList[j] = list(dati[ind, ]) }
        } else { j = j + 1; ind = vector() }
    }
    cnvList = cnvList[which(lapply(cnvList, length) > 0)]
    #CNV.p model
    p = sapply(cnvList, function(cnv) {
        temp = cnv$result - med
        v = (dati$result - med)[dati$win_id[-cnv$win_id]]
        sum(sapply(temp, function(i) sum(abs(v) < abs(i)) / length(v))) / length(temp)
    })
    cnvList = cnvList[p > p.val]
    # cnvDat = do.call(rbind, cnvList)
    # write.csv( cnvDat, file.path( pathout, 'CNV_report.csv' ), quote = FALSE, row.names = FALSE)
    return(cnvList)
}

chrPlot <- function( dati, posi, med, haarStartLevel, useVal="result", pathout){
    pdf( file.path( pathout, paste0( unique( dati$chr), '.pdf')), width = 14)
    if( is.null( posi)) {
        col = 1
    }else {
        col = rep( 1, dim(dati)[1])
        posi = do.call(rbind, posi)
        col[posi$win_id][posi$seg >= posi$segCut] = "orange"
        col[posi$win_id][posi$seg <= posi$segCut] = "yellow"
    }
    plot(dati[, useVal], pch=1, xlab='Window', ylab=useVal,
         main=paste('Chromosome',gsub("chr", "", unique(dati$chr))), col=col)

    col2 = rep( "red", dim(dati)[1])
    col2[posi$win_id][posi$seg >= posi$segCut] = "orange"
    col2[posi$win_id][posi$seg <= posi$segCut] = "yellow"
    points(dati$seg + med, type="l", col = col2)
    points(med + dati$segCut, type="l", col = 4)
    points(med - dati$segCut, type="l", col = 4)
    abline(h = med + haarStartLevel, col  = 3)
    abline(h = med - haarStartLevel, col  = 3)
    abline(h = med, lty = 2)
    # add quality tag info
    par(new=T)
    plot(x=seq_along(dati$tag), y=dati$tag, axes=F, xaxs='r',ylim=c(0, 1), xlab="",ylab="", type='p', col='lightgreen', cex.lab=0.2)
    axis(4, at = seq(from=0, to=1, 10), line=5)
    mtext("confidence",side=4)
    dev.off()
}

cnvstat <- function( cnvList, libn){
    # i = sapply(cnvList, length) > 0
    cnvRes = lapply( cnvList, function( chri_cnv) {
        if( all( sapply( chri_cnv, length) == 0) ) return("")
        cnvres = sapply( chri_cnv, function(cnv) {
			chr = as.numeric(gsub("chr", "", unique(cnv[, 'chr'])))
            if(unique(cnv[, 'chr']) == "chrX") chr = 23 
		    if(unique(cnv[, 'chr']) == "chrY") chr = 24
  	        area = cnv[, 'win_id'][c(1, length(cnv[, 'win_id']))]
			trans = tran( chr, area, round(mean(cnv[, 'seg']), digits=2) )
      	    return( trans[length(trans)]) # FIXME need standard value
        })
        return(paste(cnvres, collapse = ", "))
    })
    return( paste0("seq ", paste(cnvRes[cnvRes != ""], collapse = ", ") ) )
}

tran <- function(chr, area, segVal){
    winlap <- NULL
    if(is.null(winlap)) data('cnvDB', package='AnnoroadPD', envir = environment())

    c2 <- winlap[winlap[, 'chr'] == chr, ] # FIXME, win_xx must be present in the current or parent.envir frame
    c1 <- win_k[win_k[, 'chr'] == chr, ]
    mxc1 <- max( c1[, c('start1','end1', 'start2', 'end2')]) #fixed 6 columns
    cytoband <- cytoband[cytoband[, 1] == paste0("chr", cnvDB$chrLabel[chr]), ] #FIXME, label also need to be present

    site.s <- c2[area[1], 'start1']
    site.e <- max(c2[area[2], c('end1', 'end2')]) #FIXED 6col
    if (site.e > mxc1) site.e <- mxc1
    #find win_k ends based on win100overlap
    vs <- which((c1[, 'start1'] <= site.s) & (apply(c1[, c('end1','end2')],1,max) >= site.s)) #FIXED 6col
    ve <- which((c1[, 'start1'] <= site.e) & (apply(c1[, c('end1','end2')],1,max) >= site.e)) #FIXED 6col
    #rescaled cnv size
    nl <- round((site.e - site.s) / 1000)
    #find cytoband position ends based on win100overlap
    val.cytoband <- c( which(cytoband[, 'start'] <= site.s & cytoband[, 'end'] >= site.s)[1], which(cytoband[, 'start'] <= site.e & cytoband[, 'end'] >= site.e)[1])
    #find cytoband band ends based on val.cytoband
    val.position <- paste0(cnvDB$chrLabel[chr], cytoband[val.cytoband[1], 'cytoBand'], ifelse(val.cytoband[1] == val.cytoband[2], '', paste0("-", cnvDB$chrLabel[chr], cytoband[val.cytoband[2], 'cytoBand'])))
    result <- paste0(val.position, "(", site.s, "-", site.e, ")", " @", segVal)
    ash <- c(vs, ve, site.s, site.e, nl, result)
	return(ash)
}

cnv_surf <- function(sampleInfo, threshold = 0.2){
	data('cnvDB', package='AnnoroadPD', envir = environment())
	CNV = gsub( "seq ", "", sampleInfo[['CNV']])
    if( CNV == "") return( data.frame("Sample" = as.character(sampleInfo[['Sample']]), 
									  "Karyotype" = as.character(sampleInfo[['Karyotype']]), 
									  #"Chr.p" = as.character(sampleInfo[['Chr.p']]), 
									  "Chr.Site" = "", "Area" = "", "Size" = "", "Gain.Loss" = "", "Copy.Number" = "", 
									  "Syndrome" = "-", "Gene" = "", "DGV" = "-", "OMIM" = "-", "Combined" = "-", stringsAsFactors = F))
    CNV = do.call( rbind, lapply( strsplit( CNV, ", ")[[1]], function(cnv) {
        seq = seqex( cnv )
        len = length(seq); start = as.numeric(seq[length(seq) - 2]); end = as.numeric(seq[length(seq) - 1])
        indel = ifelse(as.numeric(seq[len]) > 0, 'duplication', 'deletion')
        cnvSurf = data.frame("Sample" = as.character(sampleInfo[['Sample']]),
                             "Karyotype" = as.character(sampleInfo[['Karyotype']]),
                             #"Chr.p" = as.character(sampleInfo[['Chr.p']]),
                             "Chr.Site" = paste0(seq[1], ':', start, '-', end),
                             "Area" = paste0(seq[1:3], collapse = ""),
                             "Size" = paste0( round((end - start) / 1000), 'kb'),
                             "Gain.Loss" = indel,
							 "Copy.Number" = seq[length(seq)], 
                             "Syndrome" = '-',
                             "Gene" = getGene(c(seq[c(1, len)], start, end), genelist=genelist, noOverlap=F), stringsAsFactors = F)
		pp = cnvref[(cnvref[, 'Chromosome'] == seq[1]) & grepl(cnvSurf[, "Gain.Loss"], cnvref[, "Type"]), ]
		pp2 = pp[(start < pp[,'Start']) & (end < pp[,'End']) & (end > pp[,'Start']) &
                         ((end - pp[,'Start']) / (pp[,'End'] - pp[,'Start']) >= threshold) |
                         (start >= pp[,'Start']) & (end <= pp[,'End']) &
                         ((end - start) / (pp[,'End'] - pp[,'Start']) >= threshold) |
                         (start >= pp[,'Start']) & (start < pp[,'End']) & (end > pp[,'End']) &
                         ((pp[,'End'] - start) / (pp[,'End'] - pp[,'Start']) >= threshold) |
                         (start < pp[,'End']) & (end > pp[,'End']), ]	
        #add symptoms to cnvSurf
        cnvSurf[, "Syndrome"] = paste0(unique(as.character(pp2[, 1])), collapse = ",")
        return( cnvSurf)
    }) )

	syndrome = CNV	
	data('dgvOmimDB', package='AnnoroadPD', envir = environment())
   	syndrome$DGV = apply(syndrome[, c("Chr.Site", "Gain.Loss")], 1, function(x) {
       	matchDatabase(seg = x[1], cnv.type = x[2], chr = dgvDB$Tab$chr, db.type = dgvDB$Tab$variantsubtype, st = dgvDB$Tab$start, end = dgvDB$Tab$end, p = .8)})
    syndrome$OMIM = apply(syndrome[, c("Chr.Site", "Gain.Loss")], 1, function(x) {
   	    matchDatabase(seg = x[1], cnv.type = x[2], chr = omimDB$Tab$Chromosome, db.type = omimDB$Tab$type, st = omimDB$Tab$Start, end = omimDB$Tab$End, Target = omimDB$Tab$Syndrome, p = .5)})
    syndrome$Combined = apply(syndrome[, c("Syndrome", "OMIM")], 1, function(x) gsub("^,|,$", "", paste(x, collapse = ",")))
   	filterOutIdx = as.numeric(gsub('(\\d+)kb', '\\1', syndrome[, 'Size'])) < 1000 & syndrome[,'Syndrome'] == '' & syndrome[,'Combined'] == '-'  & syndrome[,'OMIM'] == '-'
    syndrome <- syndrome[!filterOutIdx, ]
	return(syndrome)
}

seqex = function( ktype ){
    arm = gsub('.+([p|q]{1 }).+', '\\1', ktype)
    cn = gsub('.+@|.+ @', '\\1', ktype)
    #cn = gsub('.+\u00D7|.+ \u00D7', '\\1', ktype)
    #rest = strsplit(gsub('\u00D7.+| \u00D7.+', '', gsub('seq', '', ktype)), 'q|p|\\(|\\-|\\)')[[1]]
    rest = strsplit(gsub('@.+| @.+', '', gsub('seq', '', ktype)), 'q|p|\\(|\\-|\\)')[[1]] #DaN to a X = \u00D7
    return( c(rest[1], arm, rest[-1], cn) )
}

getGene <- function(ash, genelist=NULL, noOverlap=T){
    gene <- genelist[genelist[, 'Chromosome'] == paste0("chr", ash[1]), ]
    ash <- as.numeric( ash[ 3:4 ])
    # make it more resourceful, gene within cnv + cnv within gene or partilly overlapping
    valIdx<-(gene[, 'Start'] <= ash[1] & gene[, 'End'] >= ash[2]) | (gene[, 'Start'] >= ash[1] & gene[, 'End'] <= ash[2])
    if(!noOverlap){
        valIdx <- valIdx | (gene[, 'Start'] <= ash[1] & gene[, 'End'] >= ash[1]) | (gene[, 'Start'] <= ash[2] & gene[, 'End'] >= ash[2])
    }
    res<-ifelse(sum(valIdx) == 0, '', paste0(unique(as.character((pp<-gene[valIdx, ])[sort.list(pp[, 2]), 4])), collapse = ","))
    return(res)
}

graghToReport <- function( syndrome, pathout, libn){
    if( !file.exists( file.path(pathout, 'toReport', libn) ) ) dir.create( file.path(pathout, 'toReport', libn) )
    if( libn %in% syndrome[, "Sample"]) {
        chrs <- gsub('^(.+):\\d+-\\d+$', '\\1', syndrome[, "Chr.Site"][syndrome[, "Sample"] == libn])
        file.copy( file.path(pathout, libn, "Graphs", paste0('chr', chrs, '.pdf')), file.path(pathout, "toReport", libn, paste0(chrs, "号染色体Z值图.pdf")))
    }
    file.copy(file.path(pathout, libn, 'Boxplot.pdf'), file.path(pathout, "toReport", libn, '染色体倍数Z值箱线图.pdf'))
    file.copy(file.path(pathout, libn, 'copyNumber_p.csv'), file.path(pathout, "toReport", libn, '拷贝数P值表.csv'))
}

makeSHEET2 <- function(syndrome) {
	data('dgvOmimDB', package='AnnoroadPD', envir = environment())
	
	if(syndrome[, "Chr.Site"] == "") {
        syndrome$DGV = ""
        syndrome$OMIM = ""
        syndrome$Combined = ""
    } else {
		syndrome$DGV <- apply(syndrome[, c("Chr.Site", "Gain.Loss")], 1, function(x) {
			matchDatabase(seg = x[1], cnv.type = x[2], chr = dgvDB$Tab$chr, db.type = dgvDB$Tab$variantsubtype, st = dgvDB$Tab$start, end = dgvDB$Tab$end, p = .8)})
        syndrome$OMIM <- apply(syndrome[, c("Chr.Site", "Gain.Loss")], 1, function(x) {
            matchDatabase(seg = x[1], cnv.type = x[2], chr = omimDB$Tab$Chromosome, db.type = omimDB$Tab$type, st = omimDB$Tab$Start, end = omimDB$Tab$End, Target = omimDB$Tab$Syndrome, p = .5)})
		syndrome$Combined <- apply(syndrome[, c("Syndrome", "OMIM")], 1, function(x) gsub("^,|,$", "", paste(x, collapse = ",")))
        filterOutIdx<-as.numeric(gsub('(\\d+)kb', '\\1', syndrome[, 'Size'])) < 1000 &
            syndrome[,'Syndrome'] == '' & syndrome[,'Combined'] == '-'  & syndrome[,'OMIM'] == '-'
		syndrome <- syndrome[!filterOutIdx, ]
	}
	return(syndrome)
}

matchDatabase <- function(seg, cnv.type, chr, db.type, st, end, Target = NULL, p = .8) {
    iChr <- gsub(":+.+", "", as.character(seg))
    type <- c("loss|deletion", "duplication|gain|insertion")
    type <- type[grep(tolower(as.character(unlist(cnv.type))), type)]
    pos <- as.numeric(strsplit(gsub(".+:+", "", seg), "-")[[1]])
    indi <- which(chr == iChr & grepl(type, db.type))
    out <- c()
    for (iI in indi) {
        s1 <- max(pos[1], st[iI])
        s2 <- min(pos[2], end[iI])
        this.p <- ifelse(is.null(Target), (s2 - s1) / max(abs(pos[2] - pos[1]), 1), (s2 - s1) / max(abs(end[iI] - st[iI]), 1))
        if (this.p >= p)
            out <- c(out, iI)
    }
    if (length(out) == 0) {
		return("-")
    } else if (is.null(Target)) {
    	return("+")
    } else {
        return(paste(Target[out], collapse = ";"))
    }
}

makeAllInOneExcel <- function(SHEET1=NULL, SHEET2=NULL, SHEET3=NULL, reportName, outpath) {
	excelPathOut = file.path(outpath, "/toReport/", paste0(reportName, ".xlsx"))

    try( require(xlsx))
	# creating work book and saveing into a '.xlsx' file
	wb <- xlsx::createWorkbook()
	xlsx::saveWorkbook(wb, excelPathOut)

	if(!is.null(SHEET1)) {
		write.xlsx(SHEET1, excelPathOut, sheetName="Analysis Result", append=TRUE, row.names=FALSE)
		
	} else warning('SHEET1 not present!')
	#print("Analysis Result is OK")

	if(!is.null(SHEET2)) {
		write.xlsx(SHEET2, excelPathOut, sheetName="CNV in Detail", append=TRUE, row.names=FALSE)
    } else {
		#SHEET2 <- data.frame( Sample=character(0), Karyotype=character(0), Chr.Site=character(0),  Area=character(0), Size=character(0), Gain.Loss = character(0), Syndrome=character(0),  Gene=character(0), DGV=character(0),  OMIM=character(0), Combined=character(0))
		write.xlsx(SHEET2, excelPathOut, sheetName="CNV in Detail", append=TRUE, row.names=FALSE)
	}
	#print("CNV in Detail is OK")

    if(!is.null(SHEET3)) {
		write.xlsx(SHEET3, excelPathOut, sheetName="Qc Info", append=TRUE)
    } else warning('SHEET3 not present!')
	#print("Qc Info is OK")
}

