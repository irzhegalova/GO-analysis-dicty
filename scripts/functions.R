library(dplyr)
library(stringr)
library(tximport)
library(ggplot2)
library(tibble)
library(DESeq2)
library("RColorBrewer")
library("pheatmap")
library(tidyr)
library(GO.db)
library(GOSemSim)

# Functions

#create GO db merged with GO-term annotations
create_GO.dicty <- function() {

    #BP = Biological process; CC = Cellular component; MF = Molecular function
    GO.dicty <- read.delim("/home/izhegalova/projects/dicty_hic_loop_study/data/gene_association.dictyBase",
                         header=FALSE, comment.char="!")
    colnames(GO.dicty) <- c('DB', 'Gene', 'Gene_symbol', 'Qualifier', 'GO_id', 'DB:Reference',
                       'Evidence Code', 'With_from', 'Aspect', 'Gene_name', 'Gene_synonym',
                       'DBobject_type', 'taxon', 'Date', 'Assigned_by', 'Annp_extension', 'Gene_Product_Form_ID')
    GO.dicty.filtered <- GO.dicty %>%
        dplyr::select(Gene, GO_id, Gene_name)
    #load GO descruptions
    goterms <- data.frame(GO = GOID(GOTERM),
                        TERM = Term(GOTERM),
                        # Synonym(GOTERM),
                        # Secondary(GOTERM),
                        Def = Definition(GOTERM),
                        ONTOLOGY = Ontology(GOTERM))
    # create an equivalent to goAnno
    GO.dicty.filtered.withTerms <- merge(GO.dicty.filtered, goterms, by.x = 'GO_id', by.y='GO')
    colnames(GO.dicty.filtered.withTerms)[1] <- 'GO'
    return(GO.dicty.filtered.withTerms)
}

# compute hypergeometric test for GO terms
hypergeomGO <- function(GO.hypergeom, m) {
    d <- data.frame(numWdrawn=GO.hypergeom$x[m]-1, ## White balls drawn
                          numW=GO.hypergeom$X[m],        ## White balls
                          numB=GO.hypergeom$N[m]-GO.hypergeom$X[m],      ## Black balls
                          numDrawn=GO.hypergeom$n[m]-GO.hypergeom$x[m])    ## balls drawn

    ## calcute pvalues based on hypergeometric model
    pvalues <- apply(d, 1, function(n)
                     phyper(n[1], n[2], n[3], n[4], lower.tail=FALSE))
    return(pvalues)

  }

# adapted GO.analysis.Forclusters for gene lists
GO.analysis.ForGeneLists <- function(res.forGO = NULL, mode = 'DE', ONT = 'BP', GOfilter = NULL,
                                     computeGOfilter = FALSE) {
    #res.forGO has three columns: id - GO-term ID
#                                 cl - number of group
                        #         GENE = gene name

      # simplify GOterms
    if (computeGOfilter)
        {
        GO.dicty <- create_GO.dicty()
        GOfilter <- simplifyGOterms.dicty(levels(GO.dicty.G0SemSim@geneAnno$GO),
                        cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", ONT,
                        GO.dicty.G0SemSim)
    }

  ##------------Gene onthology-------------
    GO.dicty.filtered <- create_GO.dicty() %>%
        dplyr::filter(ONTOLOGY %in% ONT)
    res.GO <- merge(res.forGO, GO.dicty.filtered, by.x='GENE', by.y="Gene") %>%
        dplyr::filter(ONTOLOGY == ONT)

    # hypergeometric test

    #X
    GO.dicty.sum <- GO.dicty.filtered %>%
        dplyr::count(GO)
    colnames(GO.dicty.sum)[2] <- 'X'
    #x
    GO.sum <- res.GO %>%
        dplyr::count(GO)
    colnames(GO.sum)[2] <- 'x'

    #N
    if (mode == 'DE') N <- 8987
    if (mode == 'all') N <- 10350
    #n
    n <- res.forGO %>%
        nrow()
    GO.dicty.unique <- unique(GO.dicty.filtered[,c(1,4,5)])
    GO.hypergeom <- merge(GO.sum, GO.dicty.sum, by="GO")
    GO.hypergeom.annotated <- merge(GO.hypergeom, GO.dicty.unique, by="GO")

    GO.hypergeom.annotated$N <- N
    GO.hypergeom.annotated$n <- n

    # throw away too small or too large GO-terms
    GO.hypergeom.annotated <- GO.hypergeom.annotated[GO.hypergeom.annotated$X > 10, ]
    GO.hypergeom.annotated <- GO.hypergeom.annotated[GO.hypergeom.annotated$X < 200, ]
    GO.hypergeom.annotated$p.value <- sapply(1:nrow(GO.hypergeom.annotated),
                                             function(l) hypergeomGO(GO.hypergeom.annotated, l))
    # create p.adjusted and filter
    GO.hypergeom.annotated$p.adjust <- p.adjust(GO.hypergeom.annotated$p.value, method = "BH")

    # find bootstraping values
    qobj <- tryCatch(qvalue(p=GO.hypergeom.annotated$p.value, lambda=0.05,
                            pi0.method="bootstrap"), error=function(e) NULL)
    if (class(qobj) == "qvalue") {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }
    GO.hypergeom.annotated$qvalues <- qvalues
    GO.hypergeom.significant <- GO.hypergeom.annotated[GO.hypergeom.annotated$p.adjust < 0.05, ]
    GO.hypergeom.significant <- GO.hypergeom.significant[order(GO.hypergeom.significant$p.adjust), ]
    return(GO.hypergeom.significant)
}

# function to be called inside computeIC.dicty to prevent error
countGOsum <- function(i) {
    out <- tryCatch(
        {
            sum(gocount[Offsprings[[i]]], na.rm=TRUE)
        },
        error=function(cond) {
            message(paste("GO term has problems:", i))
            # Choose a return value in case of error
            return(0)
        }
        )
        return(out)
}


computeIC.dicty <- function(goAnno, ont) {
    godata <- goAnno
    goids <- unique(godata[godata$ONTOLOGY == ont, "GO"])
    ## all GO terms appearing in an given ontology
    goterms=goAnno$GO
    gocount <- table(goterms)
    ## goid of specific organism and selected category.
    goname  <- names(gocount)

    ## ensure goterms not appearing in the specific annotation have 0 frequency..
    go.diff        <- setdiff(goids, goname)
    m              <- double(length(go.diff))
    names(m)       <- go.diff
    gocount        <- as.vector(gocount)
    names(gocount) <- goname
    gocount        <- c(gocount, m)

    Offsprings <- switch(ont,
                         MF = AnnotationDbi::as.list(GOMFOFFSPRING),
                         BP = AnnotationDbi::as.list(GOBPOFFSPRING),
                         CC = AnnotationDbi::as.list(GOCCOFFSPRING))
    cnt <- gocount[goids] + sapply(goids, function(i) countGOsum(i))
    names(cnt) <- goids

    ## the probabilities of occurrence of GO terms in a specific corpus.
    p <- cnt/sum(gocount)
    ## IC of GO terms was quantified as the negative log likelihood.
    IC <- -log(p)
    return(IC)
}

godata.dicty <- function(goAnno, ont, computeIC = TRUE) {
    ont <- toupper(ont)
    ont <- match.arg(ont, c("BP", "CC", "MF"))

    goAnno <- goAnno[!is.na(goAnno$GO), ]
    goAnno <- goAnno[goAnno$ONTOLOGY == ont,]
    kk = as.character(goAnno$Gene)

    metadata <- data.frame(name=c('Db type', 'ORGANISM'),
                           value=c('Manual', 'Dictyostelium discoideum')
                           )

    res <- new("GOSemSimDATA",
               keys = kk,
               ont = ont,
               geneAnno = goAnno,
               metadata = metadata)
    if (computeIC) {
        message('preparing IC data...')
        IC <- computeIC.dicty(goAnno, ont)
        res@IC <- IC
    }

    return(res)
}

mgoSim.dicty <- function(GO1, semData, measure="Wang", combine="BMA"){
    scores <- termSim(GO1, GO1, semData, method=measure)
    # we need table, not one results
#     res <- combineScores(scores, combine)
    return(round(scores, digits=3))
}


setClass("GOSemSimDATA",
         representation = representation(
             keys = "character",
             ont = "character",
             IC = "numeric",
             geneAnno = "data.frame",
             metadata = "data.frame"
         ))

simplifyGOterms.dicty <- function(GO, cutoff=0.7, by="p.adjust", select_fun=min,
                                  measure="Rel", ontology, semData) {
     if (missing(semData) || is.null(semData)) {
        if (measure == "Wang") {
            semData <- godata(ont = ontology)
        } else {
            stop("godata should be provided for IC-based methods...")
        }
    } else {
        if (ontology != semData@ont) {
            msg <- paste("semData is for", semData@ont, "ontology, while enrichment result is for", ontology)
            stop(msg)
        }
    }
    sim <- mgoSim.dicty(GO,
                  semData = semData,
                  measure=measure)
    diag(sim) <- 0
    sim[lower.tri(sim)] <- 0
    sim[sim > cutoff] <- NA
    res <- sim[, colSums(is.na(sim)) == 0]
    return(colnames(res))
    }
