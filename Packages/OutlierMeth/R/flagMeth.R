flagMeth <-
function(beta, reference, p = 0.01)  {

options(scipen = 999)

binDat=matrix(rep(NA,prod(dim(beta))),ncol=ncol(beta))

rownames(binDat)=rownames(beta)
        colnames(binDat)=colnames(beta)

   
        #### get common CpGs
        cpgs=intersect(rownames(beta),rownames(reference))

binDat[cpgs,]=sweep(beta[cpgs,],1, reference[cpgs,paste0("P",p)],">")
binDat[cpgs,][which(sweep(beta[cpgs,],1,reference[cpgs,paste0("P",p)],">"))]=1
binDat[cpgs,][which(sweep(beta[cpgs,],1,reference[cpgs,paste0("N",p)],"<"))]=-1

return(binDat)
    }

