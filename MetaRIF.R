# Load packages
# If the package is not installed, install it
# If the package is not loaded, load it
    packages <- c("KEGGREST","MetaDCN","MetaPath","MetaDE","MetaQC","preproc","plyr","magrittr","caret", "glmnet", "lightgbm", "nnet", "randomForest", "e1071")
    lapply(packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
    })

# a function named MSig_pred_factor that predicts molecular subtypes of Recurrent Implantation Failure (RIF) 
# using gene expression data.
MSig_pred_factor=function(data){
        # load("./MetaSig_classifier_gene_plus.rda")
        data01=data[match(gene1$ENTREZID,rownames(data)),]
        data02=data[match(gene2$ENTREZID,rownames(data)),]
        rownames(data01)=paste0("gene_",rownames(data01))
        rownames(data02)=paste0("gene_",rownames(data02))
        data01[is.na(data01)] <- 0
        data02[is.na(data02)] <- 0
        data01=as.data.frame(t(data01));data02=as.data.frame(t(data02))
        if(class(best_model_RIF11)[1]!="multnet"){
                prob01 <- predict(best_model_RIF11, newdata = as.matrix(data01),type="prob")
            }else {
                prob01 <- predict(best_model_RIF11, newx = as.matrix(data01),type = "response")
                prob01<-prob01[,,1]}
            
        if(class(best_model_RIF22)[1]!="multnet"){
                prob02 <- predict(best_model_RIF22, newdata = as.matrix(data02),type="prob")
            }else {
                prob02 <- predict(best_model_RIF22, newx = as.matrix(data02),type = "response")
                prob02<-prob02[,,1]}
        prob01_out1=prob01[,1:2];colnames(prob01_out1)=c("Prob_Normal1","Prob_RIF1")
        prob02_out2=prob02[,1:2];colnames(prob02_out2)=c("Prob_Normal2","Prob_RIF2")
        prob_out=cbind(prob01_out1,prob02_out2);prob_out=as.data.frame(prob_out)
        # return(prob_out)
        pred <- ifelse(prob_out$Prob_RIF1 > prob_out$Prob_Normal1 & prob_out$Prob_RIF1 > prob_out$Prob_RIF2 + 0.1, "RIF1", 
                        ifelse(prob_out$Prob_RIF2 > prob_out$Prob_Normal2 & prob_out$Prob_RIF2 > prob_out$Prob_RIF1 + 0.1, "RIF2",
                            ifelse(prob_out$Prob_RIF1 < prob_out$Prob_Normal1 & prob_out$Prob_RIF2 < prob_out$Prob_Normal2, "Normal", "Uncertain")))
        return(pred)
    }
# a function named MSig_pred_prob that predicts molecular subtypes probability of Recurrent Implantation Failure (RIF)
# using gene expression data.
MSig_pred_prob=function(data){
        # load("./MetaSig_classifier_gene_plus.rda")
        data01=data[match(gene1$ENTREZID,rownames(data)),]
        data02=data[match(gene2$ENTREZID,rownames(data)),]
        rownames(data01)=paste0("gene_",rownames(data01))
        rownames(data02)=paste0("gene_",rownames(data02))
        data01[is.na(data01)] <- 0
        data02[is.na(data02)] <- 0
        data01=as.data.frame(t(data01));data02=as.data.frame(t(data02))
        if(class(best_model_RIF11)[1]!="multnet"){
            prob01 <- predict(best_model_RIF11, newdata = as.matrix(data01),type="prob")
        }else {
            prob01 <- predict(best_model_RIF11, newx = as.matrix(data01),type = "response")
            prob01=prob01[,,1]}
        
        if(class(best_model_RIF22)[1]!="multnet"){
            prob02 <- predict(best_model_RIF22, newdata = as.matrix(data02),type="prob")
        }else {
            prob02 <- predict(best_model_RIF22, newx = as.matrix(data02),type = "response")
            prob02=prob02[,,1]}
        prob01_out1=prob01[,1:2];colnames(prob01_out1)=c("Prob_Normal1","Prob_RIF1")
        prob02_out2=prob02[,1:2];colnames(prob02_out2)=c("Prob_Normal2","Prob_RIF2")
        prob=cbind(prob01_out1,prob02_out2);prob=as.data.frame(prob)

        pred <- ifelse(prob$Prob_RIF1 > prob$Prob_Normal1 & prob$Prob_RIF1 > prob$Prob_RIF2 + 0.1, "RIF1", 
                    ifelse(prob$Prob_RIF2 > prob$Prob_Normal2 & prob$Prob_RIF2 > prob$Prob_RIF1 + 0.1, "RIF2",
                        ifelse(prob$Prob_RIF1 < prob$Prob_Normal1 & prob$Prob_RIF2 < prob$Prob_Normal2, "Normal", "Uncertain")))

        prob_2 <- data.frame(
            Normal_prob = ifelse(pred == "Normal" | pred == "RIF1", prob$Prob_Normal1, prob$Prob_Normal2),
            RIF_prob = ifelse(pred == "Normal" |pred == "RIF1", prob$Prob_RIF1, prob$Prob_RIF2))

        return(prob_2)}
##################Data preprocessing
# For microarray data, we recommend using quartile normalized data. For RNA sequencing data, we recommend using FPKM data.
# The input data should be a matrix with rows as genes and columns as samples.
# The gene ENTREZID should be the row names of the matrix.
# The sample names should be the column names of the matrix.
MetaRIF=function(data){
    load("MetaRIF_classifier_gene_plus.rda")
    sdat=t(scale(t(data)))
    sdat<-t(scale(t(sdat)))
    pred=MSig_pred_factor(sdat)
    prob=MSig_pred_prob(sdat)
    return(list(pred=pred,prob=prob))
}


    

