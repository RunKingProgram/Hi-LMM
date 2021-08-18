QQ_plot <- function(methodname,pValue){
    print("drawing Q-Q plot")
    jpeg(file = paste(methodname,".QQ.Plot.jpeg",sep = ""),
     width = 600 * 5,height = 600 * 5,res = 72 * 4)
    par(mfrow = c(1,1),mar = c(3.1,3.1,0.7,1),oma = c(2.5,2.5,2.5,0),
     tcl=-0.5,mgp=c(2.8,1.3,0),lwd = 1.5)
    P.values <- as.matrix(pValue)
    N = nrow(P.values)
    P.values <- as.matrix((P.values)[order(P.values)])
    log.P.values <- as.matrix(rev(P.values))
    p_value_quantiles <- (1:N)/(N+1)
    log.Quantiles <- -log10(p_value_quantiles)
    N1 = length(log.Quantiles)
    c95 <- rep(NA,N1)
    c05 <- rep(NA,N1)
    for(j in 1:N1){
        k=ceiling((10^-log.Quantiles[j])*N)
        if(k==0)k=1
        c95[j] <- qbeta(0.95,k,N-k+1)
        c05[j] <- qbeta(0.05,k,N-k+1)
    }
    plot(NULL, xlim = c(0,max(log.Quantiles)),
     ylim = c(0,max(c(log.P.values,-log10(c05)))),
     cex.axis = 3.0, cex.lab = 2.2, type = "l",lty = 1, lwd = 5,
     axes = TRUE, xlab = "", ylab = "",col = "gray",yaxt = "n",xaxt = "n")
    index = length(c95):1
    polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),
    col='gray',border = NA)
    abline(a = 0, b = 1, col = "black",lwd=2)
    color <- c("firebrick2")
    points(log.Quantiles, log.P.values ,col=color,cex=1.1)
    x.lim <- max(log.Quantiles)
    y.lim <- max(c(log.P.values,-log10(c05)))
    axis(1,at = 1:ceiling(x.lim),cex.axis=2.3,
     labels = c(1:ceiling(max(x.lim))),tick=TRUE,lwd.ticks=3)
    axis(2,at = seq(1,ceiling(y.lim),3),cex.axis=2.3,
     labels = seq(1,ceiling(y.lim),3),tick=TRUE,lwd.ticks=3)
    box()
    palette("default")
    mtext(expression(Observed~~-log[10](italic(p))),side=2,cex=1.9,outer=TRUE,line=-0.4)
    mtext(expression(Expected~~-log[10](italic(p))),side=1,cex=1.9,outer=TRUE)
    dev.off()
}
Manh_plot <- function(methodname,manh){
    print("drawing manhattan plot")
    jpeg(file = paste(methodname,".Manhattan.Plot.Genomewise.jpeg",sep = ""),
     width=900*5,height=500*5,res=72*4)
    par(mfrow=c(1,1),mar=c(3.1,3.1,0.7,1),oma=c(2.5,2.5,2.5,0),
     tcl=-0.5,mgp=c(2.8,1.3,0),lwd = 1.5)
    manh0 <- manh
    cutOff <- 0.05
    manh <- matrix(as.numeric(as.matrix(manh0)),nrow(manh0),ncol(manh0))
    manh <- manh[manh[,1]!=0,]
    numMarker <- nrow(manh)
    bonferroniCutOff <- -log10(cutOff/numMarker)
    y.lim <- ceiling(max(manh[,3]))
    chmtoanalyze <- unique(manh[,1])
    nchr <- length(chmtoanalyze)
    chrcolor <- c("forestgreen","firebrick2","royalblue2","orangered1",
    "gold2","violetred2","black")
    plotcolor <- rep(chrcolor,ceiling(nchr/5))
    mypch=20
    manh <- manh[order(manh[,2]),]
    manh <- manh[order(manh[,1]),]
    ticks=NULL
    lastbase=0
    for (i in chmtoanalyze){
        index=(manh[,1]==i)
        ticks <- c(ticks, lastbase+mean(manh[index,2]))
        manh[index,2]=manh[index,2]+lastbase
        lastbase=max(manh[index,2])
    }
    x <- as.numeric(manh[,2])
    y <- as.numeric(manh[,3])
    z <- as.numeric(manh[,1])
    size = 1;ratio = 10;base = 1
    themax = ceiling(max(y))
    themin = floor(min(y))
    wd = ((y - themin + base) / (themax - themin + base)) * size * ratio
    s = size - wd / ratio / 2
    plot(y ~ x,xlab = "",ylab = "" ,ylim = c(0,y.lim),cex.axis = 2.1, cex.lab = 2.2,col = plotcolor[z],
     axes = FALSE,type = "p",pch = mypch,lwd = wd,cex = s+.5,main = "",cex.main = 2)
    abline(h = bonferroniCutOff,col = "forestgreen")
    axis(1, at = ticks,cex.axis=2.3,labels = chmtoanalyze,tick = TRUE,lwd.ticks = 3)
    axis(2, at = seq(1,floor(y.lim),3),cex.axis = 2.3,
     labels = seq(1,floor(y.lim),3),tick = TRUE,lwd.ticks = 3)
    box()
    palette("default")
    mtext(expression(Observed~~-log[10](italic(p))),side = 2,cex = 1.9,outer = TRUE,line = -0.4)
    mtext(expression(Chromosome),side = 1,cex = 1.9,outer = TRUE)
    dev.off()
}

threadHiLMM<-function(i,Data_HiLMM,thread,nmar){
    nobs=length(Data_HiLMM$ebv)
    bim.filename=paste(Data_HiLMM$filename,".bim",sep="")
    bed.filename=paste(Data_HiLMM$filename,".bed",sep="")
    outfilename="HiLMM"
    Hilmm_multithreads(as.matrix(Data_HiLMM$ebv),as.matrix(Data_HiLMM$gnewb0),
      bed.filename,paste(outfilename,i,sep="")
    	,thread,nmar,i)
  }

Data_HiLMM <- function(filename,h2 = NULL,ebv = NULL){
    bed.filename <- paste(filename, ".bed", sep = "")
    fam.filename <- paste(filename, ".fam", sep = "")
    bim.filename <- paste(filename, ".bim", sep = "")
    y <- as.numeric(read.table(fam.filename)[, 6])
    nobs=length(y)
    write.table(y,"pheno",col.names=F,row.names=F)
    cmdline <- paste("gemma -bfile",filename, " -gk 2 -o kinship")
   system(cmdline)
    tamsposeXy <-a
    if(is.null(h2) & is.null(ebv)){
        ipmin <- optimize(Loglike,c(0,1), tol = 0.05,sg=tamsposeXy$sg,
        ynew = tamsposeXy$tugy,gnew = tamsposeXy$tugx)
        h2 <- ipmin$minimum
    }
    if(is.null(ebv)){
        ebv<-gblup(y,h2)
    }
    return(list(ebv=ebv,h2=h2,filename=filename,gnewb0=tamsposeXy$gnewb0))
}

HiLMM <- function(Data_Hilmm, Test = c("Jiont","Separate") ,thread=NULL,QQ=F,Manh=F){
    if (is.null(thread)) thread <- detectCores()
    cl <- makeCluster(thread)
    y = Data_Hilmm$ebv
    bimfile=fread(paste(Data_Hilmm$filename,".bim",sep=""))
    nmar=nrow(bimfile)
    parLapply(cl,c(1:thread),threadHiLMM,Data_Hilmm,thread,nmar)
    stopCluster(cl)
    result<-c()
    for(i in 1:thread){
        result=rbind(fread(paste("HiLMM",i,sep=""),head=F),result)
    }
    return(cbind(bim,result))
    pval = pchisq(as.matrix(result[,2]), 1, lower.tail = F)
    logp <- -log10(pval)
    print("Caculation Finished")
    position <- which(pval < 0.05 / nmar)
    fwrite(cbind(bimfile[position,], logp[position]),
    "QTNs", sep = " ", quote = F)
    if(QQ) QQ_plot("HiLMM", logp)
    if(Manh) Manh_plot("HiLMM", cbind(bimfile[, c(1,4)], logp))
    Test = Test[1]
    if(!any(Test == c("Jiont","Separate"))) {
        Test = c("Jiont")
        warning("You give wrong Test name. use Jiont")
    }
    if (Test == "Jiont"){
        maxmar = length(y)
        if (maxmar > 500) maxmar = 500
        index = which( pval < 0.01)
        p0 = pval[index]
        index <- Cposi_Choice(cbind(index,p0),nmar/maxmar)[,1]
        thrd = qchisq(0.0001,1,lower.tail=FALSE)
        jiontany=fast_jiont(as.matrix(y),as.matrix(Data_Hilmm$gnewb0)
        ,as.matrix(index),thrd,paste(Data_Hilmm$filename,".bed",sep=""))
        pv=pchisq(jiontany$chi[-1],1,lower.tail=FALSE)
        position=which(pv < 0.05 / nmar)
        pval[jiontany$position[position]]<-pv[position]
        position <- which(pval < 0.05 / nmar)
        logp <- -log10(pval)
        fwrite(cbind(bimfile[position,], logp[position]),
        "jiontQTNs", sep = " ", quote = F)
        if(QQ) QQ_plot("jiontGrammar_Lambda", logp)
        if(Manh) Manh_plot("jiontGrammar_Lambda", cbind(bimfile[, c(1,4)],logp))
    }
}

