####                                       ####<< class definition >>####
setClass("CountTable",
         representation(
                        expInf="ExperimentInfo",
                        tab="data.frame",
                        tabNorm="data.frame",
                        tabList="list"
               )
         )
####                                       ####<< generics  >>####
setMethod("show", signature="CountTable",
          function(object) {
            cat("CountTable Object \n")
            if (is.null(object@expInf)){
              cat(paste("Experiment not set\n"))
            } else {
              cat(paste("Experiment:",getDescription(object@expInf),"\n"))
            }
          })
####                                       ####<< methods >>####
####                                       ####< init (initCountTable) >####
#### inits @tab from the covBedCountedHits_... files given in expInf
initCountTableFromCovBedCounts<-function(x){
  x@tab=.getCountTable(x@expInf,asDataFrame=TRUE)
  return(x)
}
#### returns the counted number of reads per gene or exon (depending on reference) for every sample
.getCountTable<-function(expInf,asDataFrame=FALSE){
  cat("Processing CoverageBed output files:\n")
  files=getCovBedCountFiles(expInf)
  
  resl=list()
  pb=txtProgressBar(min = 0, max = length(files), style = 3)
  for (i in 1:length(files)){
    setTxtProgressBar(pb, i)
    resl[[i]]=.readCovBedOutput(files[i])
  }
  close(pb)
  
  countsl=list()
  coveredFracl=list()
  for (i in 1:length(resl)){
    countsl[[i]]=resl[[i]]$numReads
    coveredFracl[[i]]=resl[[i]]$coveredFraction
  }

  resmat=do.call("cbind",countsl) ## -> resmat now matrix with numReads 
  resmatFrac=do.call("cbind",coveredFracl)
  ## set row and column names
  row.names(resmat)=resl[[1]]$attributes
  row.names(resmatFrac)=resl[[1]]$attributes
  colnames(resmat)=paste("summedCov_",getSampleDescriptions(expInf),sep="")
  colnames(resmatFrac)=paste("coveredFraction_",getSampleDescriptions(expInf),sep="")

  if (asDataFrame){
    resdf=as.data.frame(resl[[1]])[c("seqname","source","feature","start","end","score","strand","frame","attributes")]
    if(all(rownames(resmat)==resdf$attributes)){
      rownames(resmat)=NULL # prevent warnings of duplicated rownames in cbind
      resdf=cbind(resdf,resmat)
    }
    if(all(rownames(resmatFrac)==resdf$attributes)){
      rownames(resmatFrac)=NULL # prevent warnings of duplicated rownames in cbind
      resdf=cbind(resdf,resmatFrac)
    }
    if (length(unique(resdf$feature))!=1) stop("Your processed reference annotation should only contain one type of feature")
    resmat=resdf
  }
  return(resmat)
}
.readCovBedOutput<-function(filename){
  covBed=read.table(filename,
    colClasses=c(
      "character","character","character",
      "numeric","numeric",
      "character","character","character","character",
      "numeric","numeric","numeric","numeric"),
    header=FALSE,sep="\t")
  
  names(covBed)=c(
             "seqname","source","feature",
             "start","end",
             "score","strand","frame","attributes",
             "numReads","numCoveredBases","ExonLength","coveredFraction")
  return(covBed)
}
####                                       ####< gene_id, unique_id, exonNumber (addStdColsToCountTable,addTabNormToCountTable) >####
#### parses the attributes column and adds a gene_id column for convenience
addGeneIDCol<-function(x,whichTab="tab"){
  cat(paste("Adding gene id column for convenience ..."))
  if(class(x)!="CountTable") stop("addGeneIDCol: object not of class CountTable")
  gff=NULL
  if (whichTab=="tab"){ gff=x@tab }
  if (whichTab=="tabNorm"){ gff=x@tabNorm }

  gff$gene_id=
    .getAttributeField(gff$attributes,
                      getReferenceAnnotationFormat(x@expInf)[["IDFieldName"]],
                      idValSep=getReferenceAnnotationFormat(x@expInf)[["idValSep"]]
                      )
  
  if (whichTab=="tab") x@tab=gff
  if (whichTab=="tabNorm") x@tabNorm=gff
  cat("done\n")
  return(x)
}
#### code adapted from https://stat.ethz.ch/pipermail/bioconductor/2008-October/024669.html
.getAttributeField<-function(x, field, attrsep=NA, idValSep=NA) {
  if (is.na(attrsep)){
    attrsep = "; "
  }
  if (is.na(idValSep)){
    eof=regexpr(field,x)[1]+nchar(field)
    idValSep=substr(x,eof,eof)
    ##cat(paste("using \"",idValSep,"\" as idValSep\n"))
  }
      
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = idValSep, fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}
addUniqueIDCol<-function(x){
  cat("Adding uniqe id column to table for convenience...")
  x@tab$unique_id=paste(x@tab$gene_id,"_",x@tab$start,sep="")
  if(any(duplicated(x@tab$unique_id))) warning("unique id is not unique. Preprocessed annotation might be inappropriate")
  cat("done\n")
  return(x)
}
addExonNumber<-function(x,tab="tab"){
  cat("Adding exon numbers...")
  gff=NULL
  if (tab=="tab"){ gff=x@tab } else {
    if (tab=="tabNorm"){ gff=x@tabNorm } else {
      gff=x@tabList[[tab]]
    }
  }

  gffSplit=split(gff,list(gff$gene_id))
  exonNumberList=lapply(gffSplit,FUN=.getExonNumbers)
  ## do.call("rbind",...) very slow on data.frames 
  ## cmp. http://stackoverflow.com/questions/5980240/performance-of-rbind-data-frame
  tmp=lapply(exonNumberList,function(el){
    m=matrix(el[,"exonNumber"],ncol=1)
    row.names(m)=el[,"unique_id"]
    return(m)
  })
  tmp2=do.call("rbind",tmp) ## much faster !
  tmp3=data.frame(unique_id=row.names(tmp2),exonNumber=tmp2[,1])
  m=merge(gff,tmp3,by.x="unique_id",by.y="unique_id",all.x=TRUE)
  stopifnot(nrow(m)==nrow(gff))
  gff=m

  if (tab=="tab"){ x@tab=gff } else {
    if (tab=="tabNorm"){ x@tabNorm=gff } else {
      x@tabList[[tab]]=gff
    }
  }
  cat("done\n")
  return(x)
}
#### -1 for last exon
.getExonNumbers<-function(gffOfOneGene){
  rnk=rank(gffOfOneGene$start)
  rnk[rnk==max(rnk)]=-1 # -1 denotes last exon
  res=data.frame(unique_id=gffOfOneGene["unique_id"],exonNumber=rnk)
  return(res)
}
####                                       ####< RPKM and expressionInfo (addStdColsToCountTable,addTabNormToCountTable) >####
addRPKMValues<-function(x,deleteExistingRPKMCols=FALSE){
  cat("Adding RPKM values for every exon...")
  tab=x@tab
  if (any(grepl("^RPKM",names(tab))) && !deleteExistingRPKMCols){
    stop("RPKM cols already exist, use deleteExistingRPKMCols=TRUE")
  } else {
    if (any(grepl("^RPKM",names(tab)))){
      tab=tab[,-grep("^RPKM",names(tab))]
    }
  }
  
                                        # RPKM Faktor
  F=10^9
                                        # C - "number of mappable reads that fell onto the gene’s exons"
  ## here: number of reads per exon
  summedCovs=tab[,grepl("^summedCov",names(tab))]

                                        # N - "total number of mappable reads in the experiment"
  ## here: total number of mapped reads in exonic regions
  summedExonicReads=apply(summedCovs,2,sum)

                                        # L - "sum of the exons in base pairs"
  ## here: exon lengths 
  exonLengths=(tab["end"]-tab["start"])+1
  
  for (i in 1:ncol(summedCovs)){
    C=summedCovs[i]
    N=summedExonicReads[i]
    L=exonLengths
    R=F*C/(N*L)
    names(R)=sub("summedCov","RPKM",names(R));
    tab=cbind(tab,R)
  }
  x@tab=tab
  cat("done \n")
  return(x)
}
addExpressionInfo<-function(x,threshold=1){
  cat("Adding expression info...")
  gff=x@tab
  if ("expressionInfo" %in% names(gff)){
    ##cat("droping old expressionInfo column\n")
    gff=subset(gff,select=-expressionInfo)
  }
  gffSplit=split(gff,list(gff$gene_id))
  expressionInfoList=lapply(gffSplit,FUN=.getExpressionInfoRPKM,threshold)
  df=data.frame(gene_id=names(expressionInfoList),expressionInfo=unlist(expressionInfoList))
  rownames(df)=NULL
  x@tab=merge(gff,df,by.x="gene_id",by.y="gene_id",all.x=TRUE)
  cat("done\n")
  return(x)
}
copyExpressionInfoToTabNorm<-function(x){
  cat("Copying expression info to tabNorm...")
  if ("expressionInfo" %in% names(x@tabNorm)){
    ##cat("droping old expressionInfo column\n")
    x@tabNorm=subset(x@tabNorm,select=-expressionInfo)
  }
  geneInfo=x@tab[,c("gene_id","expressionInfo")]
  geneInfo=geneInfo[!duplicated(geneInfo$gene_id),]
  x@tabNorm=merge(x@tabNorm,geneInfo,by.x="gene_id",by.y="gene_id",all.x=TRUE)
  cat("done\n")
  return(x)
}
.getExpressionInfoRPKM<-function(gffOfOneGene,threshold=5){
  meanExonCovsGroup1=gffOfOneGene[,grep("^meanRPKM_group1",names(gffOfOneGene))]
  meanExonCovsGroup2=gffOfOneGene[,grep("^meanRPKM_group2",names(gffOfOneGene))]

  if (is.null(meanExonCovsGroup1) || is.null(meanExonCovsGroup2)) stop()
  if (!all(is.finite(meanExonCovsGroup1)) || !all(is.finite(meanExonCovsGroup2))) {return("infiniteValAfterDivByL")}
  lim=2/3*nrow(gffOfOneGene)
  if (sum(meanExonCovsGroup1>threshold)>lim & sum(meanExonCovsGroup2>threshold)>lim){res="twoThirdExonMeanAboveThreshold"}
  if (all(meanExonCovsGroup1>threshold) & all(meanExonCovsGroup2>threshold)){res="everyExonMeanAboveThreshold"}
  if (all(meanExonCovsGroup1>threshold) & any(meanExonCovsGroup2<threshold)){res="zeroExonMeanInGroup2"}
  if (any(meanExonCovsGroup1<threshold) & all(meanExonCovsGroup2>threshold)){res="zeroExonMeanInGroup1"}
  if (any(meanExonCovsGroup1<threshold) & any(meanExonCovsGroup2<threshold)){res="zeroExonMeanInBoth"}
  if (all(meanExonCovsGroup1<threshold) & all(meanExonCovsGroup2<threshold)){res="notExpressedInBoth"}
  if (all(meanExonCovsGroup1<threshold) & !(all(meanExonCovsGroup2<threshold))){res="notExpressedInGroup1"}
  if (!(all(meanExonCovsGroup1<threshold)) & all(meanExonCovsGroup2<threshold)){res="notExpressedInGroup2"}

  ##cat(".")
  return(res)
}
####                                       ####< flags >####
addZeroMedianInGeneExpressionSum<-function(x){
  cat("addZeroMedianInGeneExpressionSum...")
  gff=x@tab
  if ("flagZeroMedianInGeneExpressionSum" %in% names(gff)){
    ##cat("droping old flag column\n")
    gff=subset(gff,select=-flagZeroMedianInGeneExpressionSum)
  }
  gffSplit=split(gff,list(gff$gene_id))
  flagList=lapply(gffSplit,FUN=.getZeroMedianGeneExpressionSum)
  df=data.frame(gene_id=names(flagList),flagZeroMedianInGeneExpressionSum=unlist(flagList))
  rownames(df)=NULL
  m=merge(gff,df,by.x="gene_id",by.y="gene_id",all.x=TRUE)
  x@tab=m
  cat("done\n")
  return(x)
}
.getZeroMedianGeneExpressionSum<-function(gffOfOneGene){
  covs=gffOfOneGene[,grep("^summedCov_",names(gffOfOneGene))]
  return(sum(apply(covs,2,median)==0))
}
addZeroMedianInGeneExpressionFlag<-function(x){
  cat("addZeroMedianGeneExpressionFlag...")
  gff=x@tab
  if ("flagZeroMedianInGeneExpression" %in% names(gff)){
    ##cat("droping old flag column\n")
    gff=subset(gff,select=-flagZeroMedianInGeneExpression)
  }
  gffSplit=split(gff,list(gff$gene_id))
  flagList=lapply(gffSplit,FUN=.getZeroMedianGeneExpressionFlag)
  df=data.frame(gene_id=names(flagList),flagZeroMedianInGeneExpression=unlist(flagList))
  rownames(df)=NULL
  m=merge(gff,df,by.x="gene_id",by.y="gene_id",all.x=TRUE)
  x@tab=m
  cat("done\n")
  return(x)
}
.getZeroMedianGeneExpressionFlag<-function(gffOfOneGene){
  covs=gffOfOneGene[,grep("^summedCov_",names(gffOfOneGene))]
  return(any(apply(covs,2,median)==0))
}
####                                       ####< library size normalisation >####
#### divides every entry by its DESeq library size parameter to normalize for different library sizes
normaliseTabDESeqSizeFactor<-function(x,normaliseJunctionReadTab=FALSE){
  expInf=x@expInf
  sf=NULL
  
  if (!normaliseJunctionReadTab){
    smplCovs=getPartOfTab(x,select="SummedCov")
    sf=.estimateSizeFactorsForMatrix(smplCovs)
    smplCovsNorm=t(apply(smplCovs,1,function(row){row/sf}))
    x@tab[,c(grep("^summedCov_",names(x@tab)))]=as.data.frame(smplCovsNorm) 
  }
  if (normaliseJunctionReadTab){
    ## grep for desc$ to avoid e.g. finding "desc_10" when searching "desc_1"
    colsToChoose=unlist(sapply(getSampleDescriptions(expInf),function(el){grep(paste(el,"$",sep=""),colnames(x@tabList[["junctionReadTab"]]))}))
    junctionReadCounts=x@tabList[["junctionReadTab"]][,colsToChoose]
    sf=.estimateSizeFactorsForMatrix(junctionReadCounts)
    junctionReadCountsNormed=t(apply(junctionReadCounts,1,function(row){row/sf}))
    x@tabList[["junctionReadTab"]][,colsToChoose]=as.data.frame(junctionReadCountsNormed)
  }
  return(x)
}
#### Estimates library size factors as described in:
#### Anders S, Huber W.
#### Differential expression analysis for sequence count data. Genome Biology. 2010;11(10):R106.
#### Available at: http://genomebiology.com/2010/11/10/R106 [Accessed October 29, 2010].
#### Code of this method is taken from the corresponding open source package DESeq.
.estimateSizeFactorsForMatrix<-function( counts ){
   geomeans <- exp( rowMeans( log(counts) ) )
   apply( counts, 2, function(cnts) 
      median( ( cnts / geomeans )[ geomeans>0 ] ) )
}
####                                       ####< means, meds and log ratios (addGroupColsToCountTable,addTabNormToCountTable) >####
#### adds row mean values for both groups
addRowMeansForConditions<-function(x,whichTab="tab"){
  cat("Adding row means ...")
  
  if (whichTab=="tab") gff=x@tab
  if (whichTab=="tabNorm") gff=x@tabNorm
  
  covsGroup1=gff[,grep("^summedCov",names(gff))][,getGroupInfo(x@expInf)[[1]][["sampleNums"]]]
  covsGroup2=gff[,grep("^summedCov",names(gff))][,getGroupInfo(x@expInf)[[2]][["sampleNums"]]]
  gff$meanExonCovs_group1=apply(covsGroup1,1,mean)
  gff$meanExonCovs_group2=apply(covsGroup2,1,mean)

  if (whichTab=="tab"){
    covsGroup1=gff[,grep("^RPKM",names(gff))][,getGroupInfo(x@expInf)[[1]][["sampleNums"]]]
    covsGroup2=gff[,grep("^RPKM",names(gff))][,getGroupInfo(x@expInf)[[2]][["sampleNums"]]]
    gff$meanRPKM_group1=apply(covsGroup1,1,mean)
    gff$meanRPKM_group2=apply(covsGroup2,1,mean)
  }
  
  if (whichTab=="tab") x@tab=gff
  if (whichTab=="tabNorm") x@tabNorm=gff
  
  cat("done\n")
  return(x)
}
#### adds row median values for both groups
addRowMedsForConditions<-function(x,whichTab="tab"){
  cat("Adding row meds...")
  if (whichTab=="tab") gff=x@tab
  if (whichTab=="tabNorm") gff=x@tabNorm
  

  covsGroup1=gff[,grep("^summedCov",names(gff))][,getGroupInfo(x@expInf)[[1]][["sampleNums"]]]
  covsGroup2=gff[,grep("^summedCov",names(gff))][,getGroupInfo(x@expInf)[[2]][["sampleNums"]]]
  gff$medianExonCovs_group1=apply(covsGroup1,1,median)
  gff$medianExonCovs_group2=apply(covsGroup2,1,median)
  
  if (whichTab=="tab") x@tab=gff
  if (whichTab=="tabNorm") x@tabNorm=gff
  
  cat("done\n")
  return(x)
}
addLogRatioRowMeansForConditions<-function(x){
  cat("Adding log ratios...")
  x@tab$logRatioMeanGroup1_Group2=log2(x@tab$meanExonCovs_group1/x@tab$meanExonCovs_group2)
  cat("done\n")
  return(x)
}
####                                       ####< junction reads  >####
addJunctionReadTab<-function(x,nCores=1,printDotPerGene=TRUE){
  cat("====================  ...constructing junction dimensions (this will take some time)... ==============\n")
  x@tabList[["junctionReadTab"]]=calculateJunctionDimensions(x@expInf,nCores=nCores,printDotPerGene=printDotPerGene)
  cat("====================  done constructing junction dimensions ==============\n")
  return(x)
}
filterJunctionReadTab<-function(x,minSumScore=5){
  cat("Filtering junctions with few read count support...")
  jrt=x@tabList[["junctionReadTab"]]
  scores=jrt[,grep("score",names(jrt))]
  sumScores=apply(scores,1,sum)
  keep=sumScores>minSumScore
  x@tabList[["junctionReadTab"]]=x@tabList[["junctionReadTab"]][keep,]
  cat("done\n")
  return(x)
}
#### adds nExon1,nExon2,consJunction,ExonNested
addSpliceInfoToJunctionReadTab<-function(x){
  cat("addSpliceInfoToJunctionReadTab (nExon1,nExon2,consJunction,ExonNested) ...")
  expInf=x@expInf

                                        # get reference
  gffNames=c("seqname","source","feature","start","end","score","strand","frame","attributes")
  gff=x@tab[,c("gene_id","unique_id",gffNames,"exonNumber")]
  gffGeneSplit=split(gff,list(gff[,"gene_id"]))
  nExons=unlist(lapply(gffGeneSplit,nrow))

                                        # get junctionReadTab
  junctionReadTab=x@tabList[["junctionReadTab"]]
  
                                        #+ add numbers of connected exons (nExon1/nExon2 columns)
  junctionReadTab[,"nExon1"]=apply(junctionReadTab,1,function(junctionReadTabRow){return(as.numeric(strsplit(as.character(junctionReadTabRow["Exon1"]),"_")[[1]][3]))})
  junctionReadTab[,"nExon2"]=apply(junctionReadTab,1,function(junctionReadTabRow){return(as.numeric(strsplit(as.character(junctionReadTabRow["Exon2"]),"_")[[1]][3]))})

                                        # set -1 to actual number
  isLastExon=!is.na(junctionReadTab[,"nExon1"]) & junctionReadTab[,"nExon1"]==-1
  junctionReadTab[isLastExon,"nExon1"]=nExons[junctionReadTab[isLastExon,"gene_id"]]
  isLastExon=!is.na(junctionReadTab[,"nExon2"]) & junctionReadTab[,"nExon2"]==-1
  junctionReadTab[isLastExon,"nExon2"]=nExons[junctionReadTab[isLastExon,"gene_id"]]

                                        #+ add consJunction column
  consJunction=apply(junctionReadTab,1,function(row){return((as.numeric(row["nExon1"])+1)==as.numeric(row["nExon2"]))})
  consJunction[is.na(consJunction)]=FALSE 
  junctionReadTab=cbind(junctionReadTab,consJunction)
 
                                        #+ add nested exons as unique ids (ExonNested column)
  junctionReadTab[,"ExonNested"]=NA
  knownNonConsec=junctionReadTab[,"consJunction"]==FALSE & !is.na(junctionReadTab[,"Exon1"]) & !is.na(junctionReadTab[,"Exon2"])
  junctionReadTabKnownNonConsec=junctionReadTab[knownNonConsec,] 

  if (nrow(junctionReadTabKnownNonConsec) > 0){
    pb=txtProgressBar(min = 0, max = nrow(junctionReadTabKnownNonConsec), style = 3)
    for (i in 1:nrow(junctionReadTabKnownNonConsec)){
      setTxtProgressBar(pb, i)
      qGeneID=junctionReadTabKnownNonConsec[i,"gene_id"]
      nExon1=junctionReadTabKnownNonConsec[i,"nExon1"]
      nExon2=junctionReadTabKnownNonConsec[i,"nExon2"]
      junctionReadTabKnownNonConsec[i,"ExonNested"]=.getNestedExons(gff,qGeneID,nExon1,nExon2)
    }
    close(pb)
  }

  junctionReadTab[knownNonConsec,]=junctionReadTabKnownNonConsec
  
  x@tabList[["junctionReadTab"]]=junctionReadTab
  cat("done \n")
  return(x)
}
.getNestedExons<-function(gff,qGeneID,nExon1,nExon2){
  res=NULL
  if (is.na(nExon1) | is.na(nExon2)) return(NA)
  if(nExon1==nExon2){
    res=subset(gff,gene_id==qGeneID & exonNumber==nExon1,select="unique_id")
  } else {
    res=subset(gff,gene_id==qGeneID & exonNumber > nExon1 & exonNumber < nExon2,select="unique_id")
    if (nrow(res)==0) return(NA)
  }
  return(paste(res[,1],collapse=";"))
}
####                                       ####< tabNorm (addTabNormToCountTable) >####
#### normalizes count table and sets @tabNorm (@tab  is not altered)
#### #r'=asinh(((#r/e_length))/e_median))
setTabNorm<-function(x,transformFun="asinh"){
  cat("Adding tab norm...\n")
  gff=getPartOfTab(x,select="gffPlusAndSummedCov")

                                        #+ divide summedcovs by exon lengths (no multiplication with read length any more)
  gff[,c(grep("^summedCov_",names(gff)))]=gff[,c(grep("^summedCov_",names(gff)))]/(gff$end-gff$start)
  smplCovs=gff[,c(grep("^summedCov_",names(gff)))]
  smplCovs=apply(smplCovs,1:2,function(x) ifelse(is.infinite(x),NA,x)) ## some exons have length 0 !?
                                       
                                        #+ divide every exon by exon-median of its gene
  smplCovsMeds=aggregate(smplCovs,by=list(gff$attributes),FUN=median)
  names(smplCovsMeds)=sub("^summedCov","summedCovGeneMedian",names(smplCovsMeds))
  m=merge(gff,smplCovsMeds,by.x="attributes",by.y="Group.1",all.x=TRUE)

  summedCovColNames=grep("^summedCov_",names(m),value=TRUE)
  summedCovGeneMedianColNames=grep("^summedCovGeneMedian_",names(m),value=TRUE)
  normColNames=sub("summedCovGeneMedian","normedSummedCov",summedCovGeneMedianColNames)

  pb=txtProgressBar(min = 0, max = length(summedCovColNames), style = 3)  
  for (i in 1:length(summedCovColNames)){
    setTxtProgressBar(pb, i)
    numerator=summedCovColNames[i]
    denom=summedCovGeneMedianColNames[i]
    normCol=normColNames[i]
    val=m[,numerator]/m[,denom]
    val[is.infinite(val)]=NA ## median of gene = 0 -> val = Inf
    m[,normCol]=val
  }
  close(pb)

                                        #+ transform values with given function
  if(transformFun=="log2") m[,grepl("^normedSummedCov_",names(m))]=log2(m[,grepl("^normedSummedCov_",names(m))])
  if(transformFun=="asinh") m[,grepl("^normedSummedCov_",names(m))]=asinh(m[,grepl("^normedSummedCov_",names(m))])
  
                                        #+ select gff and normed columns from merge tab and rename normedSummedCov_ to summedCov_ (makes life easier later on ...)
  gffPlusCols=c("seqname","source","feature","start","end","score","strand","frame","attributes","gene_id","unique_id")
  normedCols=grep("^normedSummedCov_",names(m),value=TRUE)
  selCols=c(gffPlusCols,normedCols)
  normedTable=m[,selCols]
  names(normedTable)=sub("^normedS","s",names(normedTable),perl=TRUE)
  
  x@tabNorm=normedTable
  cat("done\n")
  return(x)
}
setMinValInTabNorm<-function(x,thresholdMin=-Inf){
  cat(paste("Setting values in tabNorm to at least",thresholdMin,"..."))
  summedCovs=x@tabNorm[,grep("^summedCov_",names(x@tabNorm))]

  summedCovsNew=list()
  for (col in 1:ncol(summedCovs)){
    summedCovsNew[[col]]=pmax(summedCovs[,col],thresholdMin)
  }
  summedCovsNew=as.data.frame(do.call("cbind",summedCovsNew))
  names(summedCovsNew)=names(summedCovs)

  x@tabNorm[,grep("^summedCov_",names(x@tabNorm))]=summedCovsNew
  cat("done\n")
  return(x)
}
####                                       ####< get  >####
getSummedCovsFromCountTable<-function(x,whichTab="tab",additionalColNames=NULL,reorderCols=TRUE){
  tab=NULL
  if (whichTab=="tab"){ tab=x@tab }
  if (whichTab=="tabNorm"){ tab=x@tabNorm }
  colsToChoose=grep("^summedCov",names(tab),value=TRUE)
  if (!all(is.null(additionalColNames))){
    colsToChoose=c(colsToChoose,additionalColNames)
  }
  if (reorderCols){ 
    colsToChoose=colsToChoose[order(colsToChoose)]
  }
  return(tab[,colsToChoose])
}
getPartOfTab<-function(x,select="gffAndSummedCov",additionalColNames=NULL,whichTab="tab"){
  columnsToSelect=NULL
  
  ## define "logical subsets"
  gffNames=c("seqname","source","feature","start","end","score","strand","frame","attributes")
  gffNames=c(gffNames,"gene_id")
  smplCovNames=grep("^summedCov_",names(x@tab),value=TRUE)

  ## combine "logical subsets"
  if (select=="SummedCov"){ columnsToSelect=smplCovNames }
  if (select=="gffAndSummedCov"){ columnsToSelect=c(gffNames,smplCovNames) }
  if (select=="gffPlusAndSummedCov"){ columnsToSelect=c(gffNames,smplCovNames,"gene_id","unique_id") }
  
  ## additional columns
  if (!is.null(additionalColNames)){ columnsToSelect=c(columnsToSelect,additionalColNames) }
  
  if (whichTab=="tab") gff=x@tab[,columnsToSelect]
  if (whichTab=="tabNorm") gff=x@tabNorm[,columnsToSelect]
  return(gff)
}
getSubtableFromGeneID<-function(x,gene_id,whichTab="tab"){
  stopifnot(is.character(gene_id) & gene_id != "NA")

  if (whichTab=="tab"){ gff=x@tab }
  if (whichTab=="tabNorm"){ gff=x@tabNorm }
  gffOfOneGene=gff[gff$gene_id==gene_id,]
  gffOfOneGene=gffOfOneGene[order(gffOfOneGene$start),]
  if (nrow(gffOfOneGene)==0){ warning(paste(gene_id,"not found in countTable"))}
  return(gffOfOneGene)
}
getExperimentInfo<-function(x){
  return(x@expInf)
}
getStdTab<-function(x){
  return(x@tab)
}
getRPKMValuesFromCountTable<-function(x,additionalColNames="unique_id",reorderCols=TRUE){
  colsToChoose=grep("^RPKM_",names(x@tab))
  colsToChoose=c(colsToChoose,grep(additionalColNames,names(x@tab)))
  if (reorderCols){ 
    colsToChoose=colsToChoose[order(colsToChoose)]
  }
  return(tab[,colsToChoose])
}
#### returns a table with gene expression values defined as the median exon expression of that gene/transcript in RPKM
getGeneExpressionTable<-function(x){
  RPKMTab=getRPKMValuesFromCountTable(x,additionalColNames="gene_id")
  RPKMTabGenes=aggregate(RPKMTab[,-1],by=list(RPKMTab[,"gene_id"]),FUN=median)
  names(RPKMTabGenes)=sub("Group.1","gene_id",names(RPKMTabGenes))
  return(RPKMTabGenes)
}
getRPKMForGene<-function(x,geneID){
  geneTab=subset(x@tab,gene_id==geneID)
  res=apply(geneTab[,grepl("^RPKM_",names(geneTab))],2,median)
  return(res)
}
getGeneExpressionInfo<-function(x){
  return(unique(x@tab[,c("gene_id","expressionInfo")]))
}
####                                       ####< plot  >####
.plotExonProfile<-function(x,geneID,eps=-10,cexLegend=1,numberLastExonAsMinusOne=FALSE,addExonCoordsToPlot=FALSE){
  stopifnot(is.character(geneID))
  if (!(geneID %in% getAllGeneIDs(x))) stop(paste(geneID,"not in countTable"))
  expInf=x@expInf

  subTab=getSubtableFromGeneID(x,geneID,whichTab="tab")
  subTabNorm=getSubtableFromGeneID(x,geneID,whichTab="tabNorm")
  
                                        # means
  meanGroup1Norm=pmax(subTabNorm[,grep("meanExonCovs_group1",names(subTabNorm))],eps)
  meanGroup2Norm=pmax(subTabNorm[,grep("meanExonCovs_group2",names(subTabNorm))],eps)

  m=subTabNorm[,paste("summedCov_",getSampleDescriptions(expInf),sep="")]
  m=apply(m,c(1,2),function(el){return(max(el,eps))})

                                        # create x axis labels
  exonNumbering=subTabNorm[,"exonNumber"]
  if (!numberLastExonAsMinusOne) exonNumbering[exonNumbering==-1]=length(exonNumbering)
  if(addExonCoordsToPlot){
    coords1=apply(subTabNorm[,c("seqname","start")],1,paste,collapse="_")
    coords2=paste("\n",subTabNorm[,c("end")],sep="_")
    coords=paste(coords1,coords2,sep="")
    exonNumbering=paste(coords,"_(",exonNumbering,")",sep="")
  }
  xlabs=exonNumbering

                                        # plot frame
  n=nrow(m)
  rownames(m)=1:nrow(m)
  yl=paste("normalised read coverage",unique(subTab[,"gene_id"],sep=""))
  ylim=c(min(m,na.rm=TRUE),max(m,na.rm=TRUE))
  plot(c(1,n),c(min(m,na.rm=TRUE),max(m,na.rm=TRUE)),cex=0,ylab=yl,main="",ylim=ylim,xaxt="n",xlab="")
  if(!addExonCoordsToPlot) axis(side=1,at=1:n,labels=xlabs,font=1,las="1")
  if(addExonCoordsToPlot) axis(side=1,at=1:n,labels=xlabs,font=1,las="2")
  title(sub="exon")


                                        # plot patients
  cols=rep(hcl(h=0,c=65,l=20,alpha=0.3),ncol(m))
  cols[getGroupInfo(expInf)[[1]]$sampleNums]=hcl(h=240,c=65,l=20,alpha=0.5)    ## alpha changed from 0.2 to 0.5
  cols[getGroupInfo(expInf)[[2]]$sampleNums]=hcl(h=0,c=65,l=20,alpha=0.5)      ## for better visability
  if(exists("i")){
    ibackup<-i
  }
  assign("i",1,envir=.GlobalEnv)
  apply(m,2,function(cl) {
    points(1:n,cl,type="b",col=cols[i],lwd=0.8,ylim=ylim,pch=0)
    assign("i",i+1,envir=.GlobalEnv)
  })
  if(exists("ibackup")){
    assign("i",ibackup,envir=.GlobalEnv)
  }
                                        # plot mean trends
  points(1:n,meanGroup1Norm,type="b",col="blue",lwd=2,ylim=ylim,pch=0)
  points(1:n,meanGroup2Norm,type="b",col="red",lwd=2,ylim=ylim,pch=0)

                                        # add legend
  #legend("top",fill=c("blue","red"),
  legend("top", inset = c(0,-0.27), bty = "n", xpd= TRUE,fill=c("blue","red"),
         legend=c(getGroupInfo(x@expInf)[[1]][["groupName"]],getGroupInfo(x@expInf)[[2]][["groupName"]]),
         cex=cexLegend,
         ncol=2
         )
}
.plotJunctionProfile<-function(x,geneID,numberLastExonAsMinusOne=FALSE,minOverallJunctionReadSupport=5){
  expInf=x@expInf
  
  bigJunctionMatrix=x@tabList[["junctionReadTab"]]
  if (is.null(bigJunctionMatrix)){cat("no junction tab in countTable\n");return()}
  junctionMatrix=bigJunctionMatrix[bigJunctionMatrix$gene_id==geneID,]
  junctionMatrix=junctionMatrix[order(junctionMatrix$splitPos1),]

  scores=junctionMatrix[,paste("score_",getSampleDescriptions(expInf),sep="")] 
  keep=as.numeric(apply(scores,1,sum))>=minOverallJunctionReadSupport
  junctionMatrix=junctionMatrix[keep,]

                                        # build xtickNames
  lastExonNumber=NULL
  if(!numberLastExonAsMinusOne){ 
    subTab=getSubtableFromGeneID(x,geneID,whichTab="tab")
    exonNumbering=subTab[,"exonNumber"]
    lastExonNumber=length(exonNumbering)
  }
  i1=unlist(lapply(strsplit(junctionMatrix$Exon1,"_"),function(el){el[1]}))
  i1=na.omit(unique(i1))[1]
  #TODO: check if el[3] is.na
  i2=unlist(lapply(strsplit(junctionMatrix$Exon1,"_"),function(el){ifelse(!is.na(el[4]),paste(el[3],el[4],sep="_"),el[3])}))
  if(!numberLastExonAsMinusOne) i2=gsub("-1",lastExonNumber,i2)
  i3=unlist(lapply(strsplit(junctionMatrix$Exon2,"_"),function(el){ifelse(!is.na(el[4]),paste(el[3],el[4],sep="_"),el[3])}))
  if(!numberLastExonAsMinusOne) i3=gsub("-1",lastExonNumber,i3)
  xtickNames=paste(i2,"-^-",i3,sep="")

                                        # normalise junction read counts
  m=junctionMatrix[,grepl("score_",names(junctionMatrix))]
  means=apply(m,2,mean)
  for (i in 1:ncol(m)){
    numerator=m[,i]
    denom=means[i]
    m[,i]=numerator/denom
  }

                                        # calculate group means
  mGroup1=m[,getGroupInfo(expInf)[[1]]$sampleNums]
  mGroup2=m[,getGroupInfo(expInf)[[2]]$sampleNums]
  meanGroup1Norm=apply(mGroup1,1,mean,na.rm=TRUE)
  meanGroup2Norm=apply(mGroup2,1,mean,na.rm=TRUE)
    
                                        # plot frame
  n=nrow(m)
  rownames(m)=1:nrow(m)
  yl="normalised number of junction reads"
  ylim=c(min(m,na.rm=TRUE),max(m,na.rm=TRUE))
  
  plot(c(1,n),c(min(m,na.rm=TRUE),max(m,na.rm=TRUE)),cex=0,ylab=yl,main="",ylim=ylim,xaxt="n",xlab="")
  axis(side=1,at=1:n,labels=xtickNames,las="2",font=1,cex.axis=1)
  title(sub="observed junctions")

                                        # plot patients
  cols=rep(hcl(h=0,c=65,l=20,alpha=0.3),ncol(m))
  cols[getGroupInfo(expInf)[[1]]$sampleNums]=hcl(h=240,c=65,l=20,alpha=0.5)  ## alpha changed form 0.2 to 0.5
  cols[getGroupInfo(expInf)[[2]]$sampleNums]=hcl(h=0,c=65,l=20,alpha=0.5)    ## for better visability
  if(exists("i")){
    ibackup<-i
  }
  for (i in 1:ncol(m)){
    points(1:n,m[,i],type="b",col=cols[i],lwd=0.8,ylim=ylim,pch=0)
  }

                                        # plot mean trends
  points(1:n,meanGroup1Norm,type="b",col="blue",lwd=2,ylim=ylim,pch=0)
  points(1:n,meanGroup2Norm,type="b",col="red",lwd=2,ylim=ylim,pch=0)
  

}
.plotJunctionHeatMap<-function(x,geneID,numberLastExonAsMinusOne=FALSE,minOverallJunctionReadSupport=0,comparisonPlot=TRUE){
  expInf=x@expInf
  
  bigJunctionMatrix=x@tabList[["junctionReadTab"]]
  if (is.null(bigJunctionMatrix)){cat("no junction tab in countTable\n");return()}
  junctionMatrix=bigJunctionMatrix[bigJunctionMatrix$gene_id==geneID,]
  junctionMatrix=junctionMatrix[order(junctionMatrix$splitPos1),]
  
  scores=junctionMatrix[,paste("score_",getSampleDescriptions(expInf),sep="")] 
  keep=as.numeric(apply(scores,1,sum))>=minOverallJunctionReadSupport
  junctionMatrix=junctionMatrix[keep,]
  
  lastExonNumber=NULL
  if(!numberLastExonAsMinusOne){ 
    subTab=getSubtableFromGeneID(x,geneID,whichTab="tab")
    exonNumbering=subTab[,"exonNumber"]
    lastExonNumber=length(exonNumbering)
  }
  
  # normalise junction read counts
  m=junctionMatrix[,grepl("score_",names(junctionMatrix))]
  means=apply(m,2,mean)
  for (i in 1:ncol(m)){
    numerator=m[,i]
    denom=means[i]
    m[,i]=numerator/denom
  }
  
  # calculate group means
  mGroup1=m[,getGroupInfo(expInf)[[1]]$sampleNums]
  mGroup2=m[,getGroupInfo(expInf)[[2]]$sampleNums]
  meanGroup1Norm=apply(mGroup1,1,mean,na.rm=TRUE)
  meanGroup2Norm=apply(mGroup2,1,mean,na.rm=TRUE)
  #replace NaN with 0
  meanGroup1Norm[is.nan(meanGroup1Norm)]<-0
  meanGroup2Norm[is.nan(meanGroup2Norm)]<-0
  
  # get exon pairings
  junctionMatrix=junctionMatrix[complete.cases(junctionMatrix$nExon1&junctionMatrix$nExon2),]
  if(nrow(junctionMatrix)<1){
    stop("Error: no junction reads found")
  }
  exonPairList <- mapply( c , junctionMatrix$nExon1, junctionMatrix$nExon2, SIMPLIFY = FALSE )
  exonPairListWODuplicates<-exonPairList
  
  #check for duplicate pairings and mark alternativ exons ends/starts
  if(length(exonPairList[duplicated(exonPairList)]) > 0)
  {
    for(i in 1:length(exonPairList[duplicated(exonPairList)])){
      duplic1 = exonPairList[duplicated(exonPairList)][[i]][1]
      duplic2 = exonPairList[duplicated(exonPairList)][[i]][2]
      if(length(unique(junctionMatrix[junctionMatrix$nExon1==duplic1&junctionMatrix$nExon2==duplic2,]$Exon1)) > 1){
        exonPairListWODuplicates[duplicated(exonPairList)][[i]][1]<-paste(exonPairList[duplicated(exonPairList)][[i]][1],"alt",sep="_")
      }
      else if(length(unique(junctionMatrix[junctionMatrix$nExon1==duplic1&junctionMatrix$nExon2==duplic2,]$Exon2)) > 1){
        exonPairListWODuplicates[duplicated(exonPairList)][[i]][2]<-paste(exonPairList[duplicated(exonPairList)][[i]][2],"alt",sep="_")
      }
    }
  }
  #if there are still duplicates repeat previous step
  i = 1
  while(length(exonPairListWODuplicates[duplicated(exonPairListWODuplicates)]) > 0){
    exonPairList=exonPairListWODuplicates
    for(i in 1:length(exonPairList[duplicated(exonPairList)])){
      duplic1 = exonPairList[duplicated(exonPairList)][[i]][1]
      duplic1 = sub("_.*","",duplic1)
      duplic2 = exonPairList[duplicated(exonPairList)][[i]][2]
      duplic2 = sub("_.*","",duplic2)
      if(length(unique(junctionMatrix[junctionMatrix$nExon1==duplic1&junctionMatrix$nExon2==duplic2,]$Exon1)) > 1){
        exonPairListWODuplicates[duplicated(exonPairList)][[i]][1]<-paste(exonPairList[duplicated(exonPairList)][[i]][1],i,sep="_")
      }
      else if(length(unique(junctionMatrix[junctionMatrix$nExon1==duplic1&junctionMatrix$nExon2==duplic2,]$Exon2)) > 1){
        exonPairListWODuplicates[duplicated(exonPairList)][[i]][2]<-paste(exonPairList[duplicated(exonPairList)][[i]][2],i,sep="_")
      }
    }
    i=i+1
    if(i>4){
      break
    }
  }
  #update exon names to junctionMatrix
  junctionMatrix$nExon1=unlist(exonPairListWODuplicates)[ c(TRUE,FALSE) ]
  junctionMatrix$nExon2=unlist(exonPairListWODuplicates)[ c(FALSE,TRUE) ]
  
  #create heat matrixes for probe 1 and 2, each field represents the coverage of an exon-exon junction
  if (any(grepl("package:gtools",search()))){
    exonNames=mixedsort(unique(unlist(exonPairListWODuplicates)))
  }
  juncMat1=matrix(nrow=length(exonNames),ncol=length(exonNames))
  juncMat2=matrix(nrow=length(exonNames),ncol=length(exonNames))
  
  rownames(juncMat1)=exonNames
  colnames(juncMat1)=exonNames
  rownames(juncMat2)=exonNames
  colnames(juncMat2)=exonNames
  
  for(i in exonNames){
    i=toString(i)
    for(j in exonNames){
      j=toString(j)
      if(i==j){
        juncMat1[i,j]=0
        juncMat2[i,j]=0
      }
      else{
        exonName=rownames(junctionMatrix[junctionMatrix$nExon1==i & junctionMatrix$nExon2==j,])
        if(length(exonName)){
          juncMat1[i,j]<-meanGroup1Norm[[rownames(junctionMatrix[junctionMatrix$nExon1==i&junctionMatrix$nExon2==j,])]]
          juncMat2[i,j]<-meanGroup2Norm[[rownames(junctionMatrix[junctionMatrix$nExon1==i&junctionMatrix$nExon2==j,])]]
        }
        else{
          juncMat1[i,j]=0
          juncMat2[i,j]=0
        }
      }
    }
  }
  
  #plot heat map, either both groups in one layout (comparisonPlot=true) or in seperate layouts
  if(comparisonPlot){
    layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),widths=c(5,1), heights=c(1,1))
    moeImagePlot(juncMat1,countTable@expInf@data$groupInfo[[1]]$groupName,"blue",cexValue=1,marValue=c(3,3,2,1))
    moeImagePlot(juncMat2,countTable@expInf@data$groupInfo[[2]]$groupName,"red",cexValue=1,marValue=c(3,3,2,1))
    layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  }
  else {
    layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
    moeImagePlot(juncMat1,countTable@expInf@data$groupInfo[[1]]$groupName,"blue",cexValue=1,marValue=c(3,3,2,1))
    moeImagePlot(juncMat2,countTable@expInf@data$groupInfo[[2]]$groupName,"red",cexValue=1,marValue=c(3,3,2,1))
  }
}
####                                       ####<< construction >>####
####                                       ####< steps  >####
#### gene_id, unique_id, exonNumber, RPKM_
addStdColsToCountTable<-function(x){
  x=addGeneIDCol(x,whichTab="tab") ## extract gene_id from attributes
  x=addUniqueIDCol(x) ## ccdsID_start
  x=addExonNumber(x,tab="tab")
  x=addRPKMValues(x)
  return(x)
}
addGroupColsToCountTable<-function(x){
  x=addRowMeansForConditions(x)
  x=addRowMedsForConditions(x)
  x=addLogRatioRowMeansForConditions(x) 
  return(x)
}
addTabNormToCountTable<-function(x){
                                        #+ set tabNorm
  x=setTabNorm(x,transformFun="asinh") ## #r'=asinh(((#r/e_length)*r_length)/e_median)
  x=setMinValInTabNorm(x,-4)
  x=addExonNumber(x,tab="tabNorm")
  x=copyExpressionInfoToTabNorm(x)
    
  x=addGeneIDCol(x,whichTab="tabNorm")
  x=addRowMeansForConditions(x,"tabNorm")
  x=addRowMedsForConditions(x,"tabNorm")
  return(x)
}
####                                       ####< checks >####
#### checks that @tab was constructed from parsed coverageBed files
checkPoint_tab<-function(x){
  expInf=x@expInf
  t1=!is.null(x@tab)
  expectedColNames=paste("summedCov_",getSampleDescriptions(expInf),sep="")
  t2=all(expectedColNames %in% colnames(x@tab))
  return(t1&t2)
}
####                                       ####<< INTERFACE >>####
##' Constructs an object of class CountTable
##'
##' @title constructCountTable
##' @param x object of class CountTable
##' @return object of class CountTable
##' @author Moritz Aschoff
##' @export
constructCountTable<-function(x,nCores=1,printDotPerGene=TRUE){

  if (nCores>1){ 
    if (any(grepl("package:parallel",search()))) {  ## check for package:parallel
      cat(paste("using",nCores,"cores "))
    } else if (any(grepl("package:multicore",search()))) { ## multicore depreciated, not recommended
      cat(paste("using",nCores,"cores "))
    } else {
      stop("Please attach parallel package when choosing nCores > 1 \n")
    }
  }
  
  x=initCountTableFromCovBedCounts(x) 
  if(!checkPoint_tab(x)) stop("Processing of coveragebed files went wrong, please check files and provided experimental information")
  
                                        #+ gene_id, unique_id, exonNumber, RPKM_
  x=addStdColsToCountTable(x) 
                                        #+ library size normalisation 
  x=normaliseTabDESeqSizeFactor(x,normaliseJunctionReadTab=FALSE) 
                                        #+ means, meds and log ratios
  x=addGroupColsToCountTable(x) 
                                         #+ juncReadTab
  x=addJunctionReadTab(x,nCores=nCores,printDotPerGene)
  x=normaliseTabDESeqSizeFactor(x,normaliseJunctionReadTab=TRUE) 
  x=addSpliceInfoToJunctionReadTab(x)
 
                                        #+ expression info for filtering
  x=addExpressionInfo(x,threshold=5)
  x=addZeroMedianInGeneExpressionFlag(x)
  x=addZeroMedianInGeneExpressionSum(x)
                                        #+ normalised reads
  x=addTabNormToCountTable(x)
  return(x)
}
####                                       ####< get >####
##' Sets the experiment information so that a CountTable object can be constructed
##'
##' @title setExperimentInfo
##' @param x object of class CountTable
##' @param expInf object of class ExperimentInfo
##' @return object of class CountTable
##' @author Moritz Aschoff
##' @export
setExperimentInfo<-function(x,expInf){
  x@expInf=expInf
  return(x)
}
##' returns a vector of gene symbols represented in CountTable object
##'
##' @title getAllGeneIDs
##' @param x object of class CountTable 
##' @return vector of gene symbols represented in CountTable object
##' @author Moritz Aschoff
##' @export
getAllGeneIDs<-function(x){     
  return(unique(x@tab[,"gene_id"]))
}
####                                       ####< plot >####
##' Plots a read count profile for an arbitrary gene represented in the CountTable object
##'
##' The read count profile is a normalised representation of the read coverage for all exons from the union transcript of that gene
##' @title plotGeneProfile
##' @param x object of class CountTable
##' @param geneID gene symbol of gene to plot
##' @param minOverallJunctionReadSupport minimal number of junction reads (summed over all samples) required to support a junction that should be plotted (defaults to 5)
##' @return no return value - plots gene profile as a side effect
##' @author Moritz Aschoff
##' @export
plotGeneProfile<-function(x,geneID,minOverallJunctionReadSupport=5){
  par(mfrow=c(2,1))
  par(mar=c(5.1, 4.1, 4.1, 2))
  .plotExonProfile(x,geneID,addExonCoordsToPlot=FALSE) ## TODO improve and add as parameter for plotGeneProfile
  par(mar=c(8, 4.1, 4.1, 2))
  .plotJunctionProfile(x,geneID,minOverallJunctionReadSupport=minOverallJunctionReadSupport)
}

