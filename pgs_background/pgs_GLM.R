source("/annoroad/data1/bioinfo/PMO/shenqingrui/PGS/PGS_RD/bin/PGS_misc_re.R")
#library( "AnnoroadPD")
library( "HaarSeg")

# PGS-GLM model:MAIN FUNCTION
PGS <- function(libDat, db, outputFilePath, reportName, useVal="u_reads") {
    "
    :type inputFilePath: list(\'winStat\' = list(data.frame), )#FIXME
    :type dbPath: string
    :type outputFilePath: string or NULL
    :rtype: void
    "
	sampleList = libDat$winStat
	#print(str(sampleList))
	#dbMat = do.call(cbind, lapply(db, function(dat) dat[,"residD"] ))
	# DB corresbond
	corrDB = list()
    for( i in 1:length( db)) corrDB[[i]] = dbCorrention(dat=db[[i]], db=db[-i], useVal = "residD")

	## MODEL 1:
	print("Modeling... ... ...")
	corrResult = lapply(names(sampleList), function(libn) {
		message(c("Processing... ", libn))
		dat = sampleList[[libn]]
		
		#tagDat: filtering low quality windows
		dat = tagDat(dat, depth = 11438000)
		#eliminateDiff: using GLM method return residuals
		dat = eliminateDiff( dat, family="poisson")
		#feedBack: correct low quality windows
		dat = feedBack(dat, gather = 40, useVal="residD")
		#dbCorrention: using data base to correcte residuals in MLM method
		corrDat = dbCorrention(dat, db, useVal = "residD")
		#saveData: including saveing:
		#
		#
		#
		saveData( corrDat, libn, corrDB, outputFilePath, useVal = "result")
		return(corrDat)
	})
	names(corrResult) = names(sampleList)

    ## MODEL 2:
	print("Analysising... ... ...")
	dir.create(file.path(outputFilePath, 'toReport'))
	SampleInfo = list()
	Syndrome = list()
	for( libn in names(corrResult)) {
		dat = corrResult[[libn]]
		#aneuNormRatio2Karyo: calculate CN and ploidy and determin karyotype
		Karyotype = aneuNormRatio2Karyo(dat, pids=2:3, cps=0:4, dist='t', useVar="result", epsilon=0.001)
		#chrPValueEstimate: estimate p value by chromosome
		#Chr.p = chrPValueEstimate(dat, libn, corrDB, outputFilePath, fractile=0.01, useVal = "result")
		#cnvSearch: search CNV in chromosome
		CNV = cnvSearch(dat, libn, minWin=25, minRun=10, maxGap=3, useVal="result", outputFilePath)
		
		#makeing EXCEL
		sampleInfo =  data.frame("Sample" = libn, 
			                     "Karyotype" = Karyotype,
    	        			     #"Chr.p" = paste0(which( Chr.p < 0.6), collapse = ","),
								 #FIXME shows the sign of chromosomes cn.p < 0.6
			                     "CNV" = CNV, stringsAsFactors = F)
		syndrome = cnv_surf(sampleInfo, threshold = 0.2)
		CNV.filted = paste0(sapply(syndrome[,"Chr.Site"][syndrome[,"OMIM"]!="-"],
							function(site) strsplit(site, ":")[[1]][1]), collapse=",")
		if(CNV.filted == "") CNV.filted = "-"
		sampleInfo[,"CNV.filted"] = CNV.filted
		# SHEET 1
		SampleInfo[[libn]] = sampleInfo[, c("Sample", "Karyotype", "CNV.filted", "CNV")]
		# SHEET 2
		Syndrome[[libn]] = syndrome
		#copy '.pdf' documents
        graghToReport(syndrome, pathout = outputFilePath, libn )
	}
    # # save data: Results/final2Details.csv
    # # save data: toReport/Details.csv
    # save data: toReport/Samples.csv
	print("Making excel... ... ...")
	SHEET1 = do.call(rbind, SampleInfo)
	SHEET2 = do.call(rbind, Syndrome)
	SHEET3 = libDat[[ 'qcStat' ]]
	write.table(SHEET1, paste0(outputFilePath, "/", reportName, ".table"), row.names = F)

	try(makeAllInOneExcel(SHEET1, SHEET2, SHEET3, reportName, outputFilePath), T)
	
	print("Batch finished")
	warnings()
## model7:
    #qcReport: generate qcReport
    # qcReport(libDat, outputFilePath)
    
    ## model8:
    #plot circos map
    #circos_plot(cytobandPath = "/home/s/Documents/Library/data/cytoband", 
    # reportPath = file.path(outputFilePath, "/toReport/Samples.csv"))
}
