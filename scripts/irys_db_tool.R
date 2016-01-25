#!/usr/bin/env Rscript
library("parallel")
library("RSQLite")

ignore_regexp<-"NONE"
only_regexp<-".*"

argv<-commandArgs(trailingOnly=TRUE)

R_threads<-10
cl<-NULL

table_prefix<-""


p<-function(...)paste(..., sep="")

dbWriteTableOrig<-dbWriteTable

dbWriteTable<-function(...) {
	while(1) {
		X<-try({dbWriteTableOrig(...)}, silent=TRUE)
		if(class(X)!="try-error")return(X)
	
		print(X)
		}
	}

load_map <- function(filename, sep="") {
  #header0 <- read.table(filename, header=FALSE, stringsAsFactors=FALSE, comment.char="", nrows=1000, fill=TRUE)
  header0 <- readLines(filename, n=200)
  header1 <- header0[regexpr("^#h", header0)>=0]
  #print(header0)
  #print(header1)
  #header<-read.table(pipe(paste("grep -E '^#h'", filename)), header=FALSE, stringsAsFactors=FALSE, comment.char="")
  header<-read.table(textConnection(header1[1]), header=FALSE, stringsAsFactors=FALSE, comment.char="")

  A<-NULL
  try({
  A <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=5000, sep=sep)
  names(A) <- unlist(header[1+(1:(dim(A)[2]))])
  })
  if(is.null(A)){
	A<-data.frame(Empty=0)[FALSE,,drop=FALSE]
	attr(A, "header")<-header0[regexpr("^#", header0)>=0]
	return(A)
	}
  #print(A[1:10,])
  
  if(dim(A)[1]>=5000) {
    B <- as.data.frame(scan(filename, what=A, comment.char="#", sep=sep), stringsAsFactors=FALSE)
    A <- B
  }
  attr(A, "header")<-header0[regexpr("^#", header0)>=0]
  return(A)
}

load_cmap<-load_map	



init_db<-function(db) {
	dbGetQuery(db, "CREATE TABLE IF NOT EXISTS manifest(path TEXT , name TEXT, type VARCHAR(255), rows INTEGER, header TEST)")
	dbGetQuery(db, "PRAGMA synchronous=OFF")
	dbGetQuery(db, "PRAGMA count_changes=OFF")
	dbGetQuery(db, "PRAGMA journal_mode=MEMORY")
	#dbGetQuery(db, "PRAGMA threads=4")
	#dbGetQuery(db, "PRAGMA cache_size=100000")
	#dbGetQuery(db, "PRAGMA temp_store=MEMORY")
	}

upload_cmap<-function(db, path, type="cmap", rootname=NULL) {
	if(is.null(rootname)) {
		#rootname<-gsub("^(.*).cmap$", "\\1", path)
		rootname<-path
		rootname<-gsub("/", "__", rootname)
		rootname<-gsub("\\.", "__", rootname)
		rootname<-p(table_prefix, rootname)
		}
	cat("Loading", path, "\n")

	Z<-file.info(path)
	if((dim(Z)[1]<1) || is.na(Z[1, "size"]) || Z[1, "size"]<1) {
		dbWriteTable(db, "manifest", data.frame(path=path, name=rootname, type=type, rows=0, header="*** ZERO LENGTH FILE"), row.names=FALSE, overwrite=FALSE, append=TRUE)
		return(NULL)
		}
	
	A<-load_cmap(path)
	
	cat("Populating", rootname, "\n")
	
	dbWriteTable(db, "manifest", data.frame(path=path, name=rootname, type=type, rows=dim(A)[1], header=p(attr(A, "header"), collapse="\n")), row.names=FALSE, overwrite=FALSE, append=TRUE)
	if(dim(A)[1]>0)
		dbWriteTable(db, rootname, A, row.names=FALSE)
	return(NULL)
	}
	
export_cmap<-function(db, filename, rootname, header=NULL) {
	cat("Exporting", rootname, "into", filename, "\n")
	df<-dbReadTable(db, rootname)
	if(is.null(header)) {
		header<-"#CMAP export"
		}
	dir.create(dirname(filename), recursive=TRUE)
	writeLines(c(header, "\n"), filename, sep="")
	write.table(df, file=filename, append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
	return(NULL)
	}

upload_xmap<-function(db, path, rootname=NULL) {
	if(is.null(rootname)) {
		#rootname<-gsub("^(.*).xmap$", "\\1", path)
		rootname<-path
		rootname<-gsub("/", "__", rootname)
		rootname<-gsub("\\.", "__", rootname)
		rootname<-p(table_prefix, rootname)
		}
	cat("Loading", path, "\n")
	if(file.info(path)[1, "size"]<1) {
		dbWriteTable(db, "manifest", data.frame(path=path, name=rootname, type="xmap", header="*** ZERO LENGTH FILE"), row.names=FALSE, overwrite=FALSE, append=TRUE)
		return(NULL)
		}
	
	A<-load_map(path)
	
	if("Alignment" %in% names(A)) {
	
		alignment<-A[,"Alignment"]
		
		parse_alignment<-function(i) {
			s<-alignment[i]
			X<-read.table(textConnection(gsub("\\)\\(", "\n", substr(s, 2, nchar(s)-1))), sep=",", header=FALSE)
			names(X)<-c("Ref", "Qry")
	
			X[,"XmapEntryID"]<-A[i, "XmapEntryID"]
			X[,"QryContigID"]<-A[i, "QryContigID"]
			X[,"RefContigID"]<-A[i, "RefContigID"]
			return(X)
			}
		
		parse_alignments<-function(idx) {
			return(do.call(rbind, lapply(idx, parse_alignment)))
			}
		
		if(length(alignment)>5000 && !is.null(cl)) {
			cat("ClusterCall\n")
			X<-clusterCall(cl, assign, "alignment", alignment)
			X<-clusterCall(cl, assign, "A", A[,c("XmapEntryID", "QryContigID", "RefContigID"),drop=FALSE])
			X<-clusterCall(cl, assign, "parse_alignment", parse_alignment)
			cat("ClusterApply\n")
			idx<-1:length(alignment)
			L<-parLapply(cl, split(idx, idx %/% 1023), parse_alignments)
			cat("rbind\n")
			df.alignment<-do.call(rbind, L)	
			#print(df.alignment[1:10,])
			} else {
			df.alignment<-do.call(rbind, lapply(1:length(alignment), parse_alignment))	
			}
		} else {
		df.alignment<-data.frame(Ref=-1, Qry=-1, XmapEntryID=-1, QryContigID=-1, RefContigID=-1)
		}
	
	cat("Populating", rootname, "\n")
	
	dbWriteTable(db, "manifest", data.frame(path=path, name=rootname, type="xmap", rows=dim(A)[1], header=p(attr(A, "header"), collapse="\n")), row.names=FALSE, overwrite=FALSE, append=TRUE)
	if(dim(A)[1]>0) {
		dbWriteTable(db, rootname, A[,names(A)[names(A)!="Alignment"],drop=FALSE], row.names=FALSE)
		dbWriteTable(db, p(rootname, "_alignment"), df.alignment, row.names=FALSE)
		}
	return(NULL)
	}
		
export_xmap<-function(db, filename, rootname, header=NULL) {
	cat("Exporting", rootname, "into", filename, "\n")
	
	df<-dbReadTable(db, rootname)
	df_align<-dbReadTable(db, p(rootname, "_alignment"))
	
	if(dim(df_align)[1]>0 && df_align[1, "XmapEntryID"]>=0) {
		L<-split(1:(dim(df_align)[1]), df_align[,"XmapEntryID"])
		
		collapse_alignment<-function(idx) {
			X<-df_align[idx,,drop=FALSE]
			X<-X[order(X[,"Ref"]),,drop=FALSE]
			return(data.frame(XmapEntryID=X[1,"XmapEntryID"], Alignment=paste("(", X[,"Ref"], ",", X[,"Qry"], ")", sep="", collapse="")))
			}
		
		collapse_alignments<-function(L) {
			return(do.call(rbind, lapply(L, collapse_alignment)))
			}
		
		if(length(L)>5000 && !is.null(cl)) {
			X<-clusterCall(cl, assign, "df_align", df_align)
			df_collapsed<-do.call(rbind, parLapply(cl, split(L, (1:length(L)) %/% 1023), collapse_alignments))
			} else {
			df_collapsed<-do.call(rbind, lapply(L, collapse_alignment))
			}
		print(df[1:5,])
		print(df_collapsed[1:5,])
		
		df<-merge(df, df_collapsed, by="XmapEntryID", all.x=TRUE)
		print(df[1:5,])
		}
	
	if(is.null(header)) {
		header<-"#XMAP export"
		}
	dir.create(dirname(filename), recursive=TRUE)
	writeLines(c(header, "\n"), filename, sep="")
	write.table(df, file=filename, append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
	return(NULL)
	}
		
	
upload_all<-function(db, dir) {
	olddir<-getwd()
	setwd(dir)
	cmaps<-list.files(getwd(), "*.cmap", recursive=TRUE)
	cmaps<-cmaps[regexpr(ignore_regexp, cmaps)<0]
	cmaps<-cmaps[regexpr(only_regexp, cmaps)>=0]
	
	cat("Loading", length(cmaps), "cmaps\n")

	if(length(cmaps)>0){
	if(is.null(cl)) {
		for(i in 1:length(cmaps))upload_cmap(db, cmaps[i])
		} else {
		X<-clusterCall(cl, setwd, dir)
		parLapply(cl, cmaps, function(name){upload_cmap(db, name)})
		
# 		fi<-file.info(cmaps)
# 		FLarge<-fi[,"size"]>10e6
# 
# 		if(sum(!FLarge)>0)parLapply(cl, cmaps[!FLarge], function(name){upload_cmap(db, name)})
# 		
# 		cmaps_large<-cmaps[FLarge]
# 		for(i in 1:length(cmaps_large)) {
# 			Y<-clusterCall(cl, function(){dbGetQuery(db, "SELECT COUNT(*) FROM manifest"); return(NULL)})
# 			upload_cmap(db, cmaps_large[i])
# 			X<-clusterCall(cl, function(){dbGetQuery(db, "SELECT COUNT(*) FROM manifest"); return(NULL)})
# 			}
		
		}
	}

	smaps<-list.files(getwd(), "*.smap", recursive=TRUE)
	smaps<-smaps[regexpr(ignore_regexp, smaps)<0]
	smaps<-smaps[regexpr(only_regexp, smaps)>=0]
	
	cat("Loading", length(smaps), "smaps\n")

	if(length(smaps)>0){
	if(is.null(cl)) {
		for(i in 1:length(smaps))upload_cmap(db, smaps[i], "smap")
		} else {
		X<-clusterCall(cl, setwd, dir)
		parLapply(cl, smaps, function(name){upload_cmap(db, name, "smap")})
		
# 		fi<-file.info(smaps)
# 		FLarge<-fi[,"size"]>10e6
# 
# 		if(sum(!FLarge)>0)parLapply(cl, smaps[!FLarge], function(name){upload_smap(db, name)})
# 		
# 		smaps_large<-smaps[FLarge]
# 		for(i in 1:length(smaps_large)) {
# 			Y<-clusterCall(cl, function(){dbGetQuery(db, "SELECT COUNT(*) FROM manifest"); return(NULL)})
# 			upload_smap(db, smaps_large[i])
# 			X<-clusterCall(cl, function(){dbGetQuery(db, "SELECT COUNT(*) FROM manifest"); return(NULL)})
# 			}
		
		}
	}
		
	xmaps<-list.files(getwd(), "*.xmap", recursive=TRUE)
	xmaps<-xmaps[regexpr(ignore_regexp, xmaps)<0]
	xmaps<-xmaps[regexpr(only_regexp, xmaps)>=0]

	cat("Loading", length(xmaps), "xmaps\n")
	if(length(xmaps)>0) {
	
	if(is.null(cl)) {
		for(i in 1:length(xmaps))upload_xmap(db, xmaps[i])
		} else {
		fi<-file.info(xmaps)
		FLarge<-fi[,"size"]>1e6
		if(sum(!FLarge)>0)parLapply(cl, xmaps[!FLarge], function(name){upload_xmap(db, name)})

		if(sum(FLarge)>0) {
			xmaps_large<-xmaps[FLarge]
			for(i in 1:length(xmaps_large)) {
				upload_xmap(db, xmaps_large[i])
				}
			}
		}
	}
	#dbGetQuery(db, "VACUUM")
	X<-clusterCall(cl, setwd, olddir)
	setwd(olddir)
	}

export_all<-function(db, dir) {
	olddir<-getwd()
	setwd(dir)

	cmaps<-dbGetQuery(db, "SELECT path, name, header FROM manifest WHERE type='cmap'")
	cmaps<-cmaps[regexpr(ignore_regexp, cmaps[,"name"])<0,,drop=FALSE]
	cmaps<-cmaps[regexpr(only_regexp, cmaps[,"name"])>=0,,drop=FALSE]
	#cmaps<-cmaps[regexpr(".cmap$", cmaps[,"path"])>=0,,drop=FALSE][1:2,]
	
	if(dim(cmaps)[1]>0)
	for(i in 1:(dim(cmaps)[1])) {
		if(cmaps[i, "rows"]<1) {
			writeLines(cmaps[i, "header"], cmaps[i, "path"])
			next
			}
		export_cmap(db, cmaps[i, "path"], cmaps[i, "name"], cmaps[i, "header"])
		}
	
	xmaps<-dbGetQuery(db, "SELECT path, name, header FROM manifest WHERE type='xmap'")
	xmaps<-xmaps[regexpr(ignore_regexp, xmaps[,"name"])<0,,drop=FALSE]
	xmaps<-xmaps[regexpr(only_regexp, xmaps[,"name"])>=0,,drop=FALSE]
	#xmaps<-xmaps[regexpr(".xmap$", xmaps[,"path"])>=0,,drop=FALSE]
	
	if(dim(xmaps)[1]>0)
	for(i in 1:(dim(xmaps)[1])) {
		if(xmaps[i, "rows"]<1) {
			writeLines(xmaps[i, "header"], xmaps[i, "path"])
			next
			}
		export_xmap(db, xmaps[i, "path"], xmaps[i, "name"], xmaps[i, "header"])
		}
	}
	
merge_db<-function(path1, paths) {
	db<-dbConnect(sqldriver(), dbname=path1)
	
	init_db(db)
	
	manifest1<-dbGetQuery(db, "SELECT * FROM manifest")
	
	if(dim(manifest1)[1]>0)
	for(path2 in paths) {
		cat("Checking", path2, "\n")
		dbGetQuery(db, p("attach \"", path2, "\" as toMerge"))
		

		manifest2<-dbGetQuery(db, "SELECT * FROM toMerge.manifest")
		manifest2<-manifest2[regexpr(ignore_regexp, manifest2[,"name"])<0,,drop=FALSE]
		manifest2<-manifest2[regexpr(only_regexp, manifest2[,"name"])>=0,,drop=FALSE]
		
		Fname<-manifest2[,"name"] %in% manifest1[,"name"]
		Fpath<-manifest2[,"path"] %in% manifest1[,"path"]
		if(any(Fname) || any(Fpath)) {
			cat("Collision while merging", path2, "into", path1, ":\n")
			print(manifest1[Fname | Fpath, c("path", "name"),drop=FALSE])
			print(manifest2[Fname | Fpath, c("path", "name"),drop=FALSE])
			dbDisconnect(db)
			return(NULL)
			}
		dbGetQuery(db, p("detach toMerge"))
		}
		
	for(path2 in paths) {
		cat("Merging", path2, "")
		dbGetQuery(db, p("attach \"", path2, "\" as toMerge"))
		
		manifest2<-dbGetQuery(db, "SELECT * FROM toMerge.manifest")
		manifest2<-manifest2[regexpr(ignore_regexp, manifest2[,"name"])<0,,drop=FALSE]
		manifest2<-manifest2[regexpr(only_regexp, manifest2[,"name"])>=0,,drop=FALSE]
		
		
		dbGetQuery(db, "INSERT INTO manifest SELECT * FROM toMerge.manifest")
		dbGetQuery(db, p("BEGIN"))
		
		if(dim(manifest2)[1]>0)
		for(i in 1:(dim(manifest2)[1])) {
			cat(".")
			#cat("Copying", manifest2[i, "name"], "\n")
			name<-manifest2[i, "name"]
			if(manifest2[i, "rows"]<1)next
			
			dbGetQuery(db, p("CREATE TABLE ", name, " AS SELECT * FROM toMerge.", name, " WHERE 0"))
			dbGetQuery(db, p("INSERT INTO ", name, " SELECT * FROM toMerge.", name))
			
			if(manifest2[i, "type"]=="xmap") {
				name<-p(manifest2[i, "name"], "_alignment")
				dbGetQuery(db, p("CREATE TABLE ", name, " AS SELECT * FROM toMerge.", name, " WHERE 0"))
				dbGetQuery(db, p("INSERT INTO ", name, " SELECT * FROM toMerge.", name))
				}
			}
		dbGetQuery(db, p("COMMIT"))
		dbGetQuery(db, p("detach toMerge"))
		cat("\n")
		}
	#cat("Cleaning up\n")
	#dbGetQuery(db, "VACUUM")
	dbDisconnect(db)
	return(NULL)
	}
	
sqldriver<-function() {
	return(dbDriver("SQLite"))
	}
	
if(!("worker" %in% ls())) {
	
	if(length(argv)<3) {
		cat("\nUsage: \n")
		cat("\n\tirys_db_tool.R [options] load /path/db.sqlite /path/input_directory\n")
		cat("\n\tirys_db_tool.R [options] export /path/db.sqlite /path/\n\n")
		cat("\n\tirys_db_tool.R [options] merge /path/target_db.sqlite /path/to_merge_db.sqlite [/path/db2.sqlite ... ]\n\n")
		cat("\t\toptions: [--ignore-regexp=NONE] [--only-regexp=.*]\n\n")
		cat("\t\tYou can use merge with only one argument to create a subset of existing database.\n\n")
		q(save="no")
		}

	set_vars<-list()

	for(i in 1:length(argv)) {
		if(argv[i]=="--")break
		if(regexpr("^--", argv[i])<0)break
			
		var_name<-gsub("-", "_", gsub("^--([^=]*)=.*", "\\1", argv[i]))
		value<-gsub("^--[^=]*=", "", argv[i])
		
		old_value<-"UNSET"
		try({
			old_value<-get(var_name)
			})
		if(length(old_value)==1 && (old_value=="UNSET")) {
			cat("Unknown option: ", argv[i], "\n")
			q(save="no", status=-1)
			}
		

		if(class(old_value)=="numeric" || class(old_value)=="integer")value<-unlist(strsplit(value, ","))
			
		assign(var_name, as(value, class(old_value)), envir=globalenv())
		set_vars[[length(set_vars)+1]]<-var_name
		}
		
	argv<-argv[i:length(argv)]

	if(tolower(argv[[1]]) %in% c("load", "export")) {

		cat("Starting R cluster\n")
		port<-10187
		cl <- NULL
		while(is.null(cl)) {
			try({
				cl <- makeCluster(getOption("cl.cores", R_threads), useXDR=TRUE, port=port)
				})
			port<-port+1
			Sys.sleep(1)
			}
		X<-clusterCall(cl, setwd, getwd())
		X<-clusterCall(cl, function(){assign("cl", NULL, env=globalenv())})
		X<-clusterCall(cl, function(){library("RSQLite")})
		clusterExport(cl, c("table_prefix", "p", "load_map", "load_cmap", "upload_cmap", "upload_xmap", "sqldriver", "init_db"))
	# 	for(name in c("table_prefix", "p", "load_map", "load_cmap", "upload_cmap", "upload_xmap", "sqldriver")) {
	# 		X<-clusterCall(cl, function(a, b){assign(a, b, env=globalenv())}, name, get(name))
	# 		}
		}

	if(tolower(argv[[1]])=="load") {

		db<-dbConnect(sqldriver(), dbname=argv[[2]])
		init_db(db)
		
		#X<-clusterCall(cl, function(x){assign("db", dbConnect(sqldriver(), dbname=x), env=globalenv())}, argv[[2]])
		X<-parLapply(cl, 1:length(cl), function(k, x){assign("db", dbConnect(sqldriver(), dbname=p(x, "_node", k)), env=globalenv()); init_db(db)}, argv[[2]])


		upload_all(db, argv[[3]])		

		dbDisconnect(db)
		
		print(warnings())
		
		merge_db(argv[[2]], p(argv[[2]], "_node", 1:length(cl)))

		print(warnings())
		
		file.remove(p(argv[[2]], "_node", 1:length(cl)))

		cat("END OF OUTPUT\n")
		
		} else
	if(tolower(argv[[1]])=="export") {
		db<-dbConnect(sqldriver(), dbname=argv[[2]])
		X<-clusterCall(cl, function(x){assign("db", dbConnect(sqldriver(), dbname=x), env=globalenv())}, argv[[2]])
		
		dir.create(argv[[3]], recursive=TRUE)
		
		export_all(db, argv[[3]])
		
		dbDisconnect(db)
		
		print(warnings())

		cat("END OF OUTPUT\n")
		} else
	if(tolower(argv[[1]])=="merge") {
		merge_db(argv[[2]], argv[3:length(argv)])
		print(warnings())
		cat("END OF OUTPUT\n")
		} else {
		cat("Unknown command:", argv[[1]], "\n\n")
		}
	#upload_cmap(db, "contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap")

	#upload_xmap(db, "contigs/exp_refineFinal1/alignref/EXP_REFINEFINAL1.xmap")
	}

