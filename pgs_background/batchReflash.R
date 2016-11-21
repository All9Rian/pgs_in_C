#!/annoroad/data1/bioinfo/PMO/shenqingrui/software/R-3.2.3/bin/Rscript
doList = read.table( "/annoroad/data1/bioinfo/PMO/shenqingrui/PGS/todolist", header=F, stringsAsFactors=F)[[1]]

dataPath = "/annoroad/data1/bioinfo/PROJECT/RD/Medical/PD/PGS/"
dbPath = "/annoroad/data1/bioinfo/PMO/shenqingrui/PGS/PGS_RD/example/dbPot1.1/"
source('/annoroad/data1/bioinfo/PMO/shenqingrui/PGS/PGS_RD/bin/pgs_GLM.R')
#source('/annoroad/data1/bioinfo/PMO/shenqingrui/PGS/PGS_RD/bin/circos.R')

tr = lapply( doList, function( libName) {
	path_100K = paste0( dataPath, libName, "/wangxiaolin")
	batchNames = dir( path_100K)[grepl( "result", dir( path_100K))]

	lapply( batchNames, function( batchName) {
		#loadData: load sample100K documents
		inputFilePath = paste0( path_100K, "/", batchName, "/analysis")
		libDat = loadData( inputFilePath)
		#updateDB: update data base
		#updateDB(dbPath, gather = 40, useVal = "residD", depth = 11438000, family="poisson")
		load( file.path( dirname( dbPath), "/db.Rda"))
		#create output file path
		outputFilePath = paste0( dataPath, libName, "/shenqingrui/", batchName)
		try(system(paste0('rm -r ', outputFilePath)), T)
		dir.create(outputFilePath)
		print("inputFilePath in:"); message(inputFilePath)
		print("dbPath in:"); message(dbPath)
		print("pathout in:"); message(outputFilePath)

		PGS( libDat=libDat, db=db, outputFilePath=outputFilePath, reportName = paste0(libName, '-', batchName), useVal="u_reads")
		#try(PGS( libDat=libDar, db=db, outputFilePath=outputFilePath, reportName = paste0(libName, '-', batchName), useVal="u_reads"), T)
	})
})
