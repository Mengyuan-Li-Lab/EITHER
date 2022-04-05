#' @title Evaluation of the tumor heterogeneity level.
#'
#' @description
#' \code{EITHER} EITHER evaluates the tumor heterogeneity level of each tumor sample based on Shannon Entropy Heterogeneity score vignette.
#'
#' @details
#' EITHER evaluates the tumor heterogeneity level of each tumor sample based on gene expression profiles.
#' 
#' @param exp gene expression profiles in both tumor and normal samples (microarray or RNA-Seq data, log2 transformed is recommended).
#' EITHER function support 2 kind of datasets: 
#' 1) if the "exp" dataset has normal samples and tumor samples, we can caculate the heterogeneity score with both normal samples and tumor samples.
#' 2) else if the "exp" dataset only has tumor samples, we can caculate the heterogeneity score by with tumor samples.
#' @param match samples in "exp" denoted by "Tumor" or "Normal".
#' @import
#' @export
#' @return A matrix with 2 columns:
#' \item{Sample}{Tumor samples to be predicted.}
#' \item{heterogeneity score}{The heterogeneity score of each sample.}
#' @author Mengyuan Li <300562@@njucm.edu.cn>
#' @examples
#' EITHER(exp, match)
EITHER <- function(exp, match)
{#Two files need to be input into this function.
  exp <- as.matrix(exp); 
  gene <- exp[-1, 1]; 
  samp <- exp[1, -1]; 
  exp <- as.matrix(exp[-1, -1]); 
  colnames(exp) <- samp; 
  rownames(exp) <- gene; 
  match <- as.matrix(match); 
  #Pick up normal samples.
  nor_pos <- which(samp %in% match[which(match[, 2] == "Normal"), 1]) 
  #Pick up tumor samples.
  tum_pos <- which(samp %in% match[which(match[, 2] == "Tumor"), 1]) 
  exp_tum <- exp[, tum_pos]; samp_tum <- samp[tum_pos]; 
  score <- matrix(0, nrow <- dim(exp_tum)[1], ncol <- dim(exp_tum)[2])

  if(length(nor_pos) > 0){
    nor <- c(); 
	for(j in 1:dim(exp)[1]){
	  nor[j] <- mean(as.numeric(exp[j, nor_pos]))
	  }#Caculate the average values of each gene in normal sample.
    #Caculate the heterogeneity score of each gene.
    for(s in 1 : dim(exp_tum)[1]){
	  for(u in 1 : dim(exp_tum)[2]){
	    score[s, u] <- (as.numeric(exp_tum[s, u]) - as.numeric(nor[s]))^2
	   }
	}
  }else if(length(nor_pos)==0){
	for(s in 1 : dim(exp_tum)[1]){ 
	  for(u in 1 : dim(exp_tum)[2]){
	    score[s, u] <- (as.numeric(exp_tum[s, u]) - mean(as.numeric(exp_tum[s, ])))^2
	  }
	}
  }
  colnames(score) <- samp_tum; rownames(score) <- gene;
    for(s in 1:dim(score)[1]){
      for(u in 1:dim(score)[2]){
        score[s, u]=as.numeric(score[s, u])%/%(5)+1
		}
	}
  heterogeneity_score <- c();
  for(z in 1:length(samp_tum)){
    p_data <- unique(score[,z]); 
	p <- c();
	for(st in 1:length(p_data)){
	  p[st]=length(which(score[,z]==p_data[st]))/dim(score)[1]
	  }
	  for(v in 1:length(p)){
	     if(p[v]==0){sm <- 0}else{sm <- p[v]*log2(p[v])}
         heterogeneity_score[z] <- heterogeneity_score[z]-sm
		 }
    }
  heterogeneity_score <- cbind(samp_tum, heterogeneity_score); 
  #caculate the heterogeneity score of each sample.
  colnames(heterogeneity_score) <- c("sample", "heterogeneity score")
  #EITHER function will output the heterogeneity score of each tumor sample.
  return(heterogeneity_score)
}





