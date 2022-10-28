expectedVsObserved <- function(data,predictor){

  data<-as.data.frame(data)
  predicted_positions <- paste(predictor,"_disordered",sep = "")
  
  phosphoDiso_ST <- list()
  phosphoDiso_obs <- numeric()
  phosphoDiso_expct <- numeric()
  phosphoDiso_expct_prob <- numeric()
  for (i in 1:nrow(data)) {
    TStotalIndexes <- as.numeric(gregexpr("S|T", data[i,"sequence"])[[1]]) 
    TSinDiso_count <- sum(TStotalIndexes %in% data[i,predicted_positions][[1]])
    TSinDiso_fraction <- TSinDiso_count/length(TStotalIndexes)
    phosphoDiso_ST[[i]] <- TStotalIndexes
    phosphoDiso_expct_prob[i] <- TSinDiso_fraction
    phosphoDiso_expct[i] <- data[i,"psites_count"][[1]]*TSinDiso_fraction
    phosphoDiso_obs[i] <- sum(data[i,"psites"][[1]] %in% data[i,predicted_positions][[1]])
  }
  return(list("phosphoDiso_ST"=phosphoDiso_ST,"phosphoDiso_obs"=phosphoDiso_obs,"phosphoDiso_expct"=phosphoDiso_expct,"phosphoDiso_expct_prob"=phosphoDiso_expct_prob))
}