#
#
#
######## the following function returns quantile plot object in ggplot2  #########
######## as well as statistics of the quantiles  #########
# 
#
#

quantile_plot_and_stats <- function(prs = prs1, 
                                    pheno.data.all = PTSD_all_pheno_data,
                                    pcs = pcs1,
                                    binary.pheno = T,
                                    pheno.name = "M3_PCL5",
                                    num.quantiles = 10, 
                                    quant.ref = 1, 
                                    covary = T, 
                                    covariates = c("C1", "C2"))
{
  
  require(ggplot2); require(plyr)
  
  if (binary.pheno == F)
  {
    names(prs) <- c("ID", "SCORE");
    pheno.data <- pheno.data.all[,c("IID", pheno.name)];
    names(pheno.data) <- c("ID", "PHEN");
    
    for.quantiles <- merge(x = pheno.data, by.x = "ID", y = prs, by.y = "ID")
    
    for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
    for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))
    
    if(covary == F){
      coef.quantiles <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1]
      ci.quantiles.u <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
      ci.quantiles.l <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
    }
    if(covary == T){
      for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
      coef.quantiles <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1]
      ci.quantiles.u <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
      ci.quantiles.l <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
    }
    
    coef.quantiles[1] <- 0 
    ci.quantiles.u[1] <- 0
    ci.quantiles.l[1] <- 0
    quant.list <- seq(1, num.quantiles, 1)
    quant.list <- quant.list[quant.list != quant.ref]     
    quantiles.for.table <- c(quant.ref, quant.list)
    quantiles.df <- data.frame(coef.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
    names(quantiles.df) <- c("Coef", "CI.U", "CI.L", "DEC")
    quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
    quantiles.plot <- ggplot(quantiles.df) + 
      geom_point(aes(x = DEC, y = Coef), colour = "royalblue", size=4) + 
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.line = element_line(colour = "black",size=0.5)) + 
      geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = Coef, x = DEC), colour = "royalblue2", size = 0.9) + 
      ylab("Change in PCL") + 
      xlab("Quantile of PTSD-PRS") #+ 
    #ylim(0, 200) +
    scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
    
  }
  
  
  ##########################################################################################
  ######################### Binary phenotype
  
  if (binary.pheno == T)
  {
    names(prs) <- c("ID", "SCORE");
    pheno.data <- pheno.data.all[,c("IID", pheno.name)]; 
    names(pheno.data) <- c("ID", "PHEN");
    pheno.data$PHEN <- factor(pheno.data$PHEN)
    
    for.quantiles <- merge(x = pheno.data, by.x = "ID", y = prs, by.y = "ID")
    
    for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
    for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))
    
    if(covary == F){
      or.quantiles <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1])
      ci.quantiles.u <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
      ci.quantiles.l <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
    }
    if(covary == T){
      for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
      or.quantiles <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1])
      ci.quantiles.u <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
      ci.quantiles.l <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
    }
    or.quantiles[1] <- 1 
    ci.quantiles.u[1] <- 1
    ci.quantiles.l[1] <- 1
    quant.list <- seq(1, num.quantiles, 1)
    quant.list <- quant.list[quant.list != quant.ref]     
    quantiles.for.table <- c(quant.ref, quant.list)
    quantiles.df <- data.frame(or.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
    names(quantiles.df) <- c("OR", "CI.U", "CI.L", "DEC")
    quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
    quantiles.plot <- ggplot(quantiles.df) + 
      geom_point(aes(x = DEC, y = OR), colour = "royalblue2", size=4) + 
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.line = element_line(colour = "black",size=0.5)) + 
      geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = OR, x = DEC), colour = "royalblue2", size = 0.9) + 
      ylab("Odds Ratio") + 
      xlab("Quantile of PTSD-PRS") + 
      scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
  }
  
  SE1 = (quantiles.df[,"CI.U"] - quantiles.df[,1])/1.96;
  quantiles.df <- cbind(quantiles.df,SE1)
  print(SE1)
  print(quantiles.df)
  return(quantiles.plot)
}

