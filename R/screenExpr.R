screenExpr <-
function(Yexpr, sdCutoff) {
    # Data set is presumed to contain non-unique geneIDs, so cannot read geneIDs into rownames
    data.all <- Yexpr
    
    print(paste("    Read expression data for ",dim(data.all)[1]," genes.",sep=""));
    
    # Remove Affy probes that do not map to geneIDs
    data.mapped <- data.all[toupper(data.all[,1]) != "UNMAPPED" & !is.na(data.all[,1]),];
    print(paste("    Mapped data to ",dim(data.mapped)[1]," genes IDs.",sep=""));
    
    # Restrict analysis to genes with sample var above threshold
    nLastCol <- dim(data.mapped)[2];
    print("Removing genes with low sample variance...");
    sampleVar <- apply(data.mapped[,2:nLastCol],1,var);	# compute sample variance
    idxHighVar <- sampleVar > sdCutoff^2;
    data.highVar <- data.mapped[idxHighVar,]
    sampleVar.highVar <- sampleVar[idxHighVar];
    print(paste("    Found ",dim(data.highVar)[1]," genes exceeding SD threshold of ",sdCutoff,".",sep=""));
    
    nGenesHighVar <- dim(data.highVar)[1];
    # Cannot have replicate gene names for clustering; choose gene with largest sample variance
    genesUnique <- as.vector(unique(data.highVar[,1]));
    nGenesUnique <- length(genesUnique);
    nSamples <- dim(data.highVar)[2]-1;		# first col is still gene name
    data <- array(dim=c(nGenesUnique,nLastCol));
    data[,1] <- genesUnique;
    colnames(data) <- colnames(data.highVar);
    print("Removing duplicate genes (selecting for max standard deviation)...");
    
    for (gene in genesUnique) {
        # index/indices of all genes matching gene (detect duplicates)
        idxGenes <- seq(along=1:nGenesHighVar)[data.highVar[,1]==gene];
        
        data.slice <- data.highVar[idxGenes,2:nLastCol];
        if (length(idxGenes)>1) {
            idxMaxVar <- which.max(sampleVar.highVar[idxGenes]);	# find dupls with max var
            data.slice <- data.slice[idxMaxVar,];
        }
        data[data[,1]==gene,2:nLastCol] <- as.matrix(data.slice);
    }
    print(paste("    ",nGenesUnique," unique genes IDs remain.",sep=""));
    print("Screened data is the output:... ");
    print("");
    rm(data.all,data.mapped,data.highVar);
    sdata<-data
    sdata
}
