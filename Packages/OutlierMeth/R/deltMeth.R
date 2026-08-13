deltMeth <-
function(beta,reference=c(all,tcga,geo,blood), p=(0.01)) {
options(scipen = 999)

binDat=matrix(rep(NA,prod(dim(beta))),ncol=ncol(beta))

rownames(binDat)=rownames(beta)
        colnames(binDat)=colnames(beta)

   
        #### get common CpGs
        cpgs=intersect(rownames(beta),rownames(reference))

binDat[cpgs,]=sweep(beta[cpgs,],1, reference[cpgs,paste0("P",p)],">")
binDat[cpgs,][which(sweep(beta[cpgs,],1,reference[cpgs,paste0("P",p)],">"))]=1
binDat[cpgs,][which(sweep(beta[cpgs,],1,reference[cpgs,paste0("N",p)],"<"))]=-1


binDat_Pos=as.data.frame(binDat==1)*1
Beta_Map = (beta * binDat_Pos)
Beta_Map[Beta_Map == 0] <- NA

Beta_Map_ID= cbind(row.names(Beta_Map),Beta_Map)
all_ID=as.data.frame(cbind(row.names(reference),(reference[,paste0("P",p)])))

names(Beta_Map_ID) [1] <- 'ID'
names(all_ID) [1] <- 'ID'
all_ID$V2=as.numeric(all_ID$V2)

Beta_Map1_thr=merge(Beta_Map_ID,all_ID,all=T)
row.names(Beta_Map1_thr) <- Beta_Map1_thr$'ID'

Beta_Map_delt_Pos= ((Beta_Map1_thr[,2:(length(Beta_Map_ID))] - (Beta_Map1_thr[,(length(Beta_Map1_thr))])))

Beta_Map_delt_Pos [is.na(Beta_Map_delt_Pos)]=0


#---------Neg

binDat_Neg=as.data.frame(binDat==-1)*1
Beta_Map = (beta * binDat_Neg)
Beta_Map[Beta_Map == 0] <- NA

Beta_Map_ID= cbind(row.names(Beta_Map),Beta_Map)
all_ID=as.data.frame(cbind(row.names(reference),(reference[,paste0("N",p)])))

names(Beta_Map_ID) [1] <- 'ID'
names(all_ID) [1] <- 'ID'
all_ID$V2=as.numeric(all_ID$V2)

Beta_Map1_thr=merge(Beta_Map_ID,all_ID,all=T)
row.names(Beta_Map1_thr) <- Beta_Map1_thr$'ID'

Beta_Map_delt_Neg= (((Beta_Map1_thr[,(length(Beta_Map1_thr))])-(Beta_Map1_thr[,2:(length(Beta_Map_ID))])))

Beta_Map_delt_Neg [is.na(Beta_Map_delt_Neg)]=0

#combine---------

Beta_Map_delt = (Beta_Map_delt_Pos-Beta_Map_delt_Neg)
Beta_Map_delt[Beta_Map_delt == 0] = NA
                
return(Beta_Map_delt)
    }


