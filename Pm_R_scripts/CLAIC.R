CLAIC <- function(J, H, LL){
  #specify J = variability matrix; 
  # H = Hessian matrix; 
  # LL = Log-likelihood 
  require(matlib)
  if(det(H) == 0){
    stop("Hessian matrix is singular!") #verifies AIC is calculable
  }
  else{
    JH <- J %*% inv(H); #calculates the matrix product of J, H^{-1}
    CLAIC <- ((2*tr(JH)) - (2*LL));
    return(CLAIC)
  }
}

setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/1D_demographics")

#had to modify Godambe.py to get it to print out the H and J matrices

#dadi-cli StatDM --fs ../dadi_1pop_wsaf.fs --model three_epoch --demo-popt three_epoch.demog.params.InferDM.bestfits --grids 80 90 100 --bootstrapping-dir bootstraps/ --output three_epoch.godambe.ci --nomisid > three_epoch_matrices.txt

H <- data.table::fread("test_H.txt") |> as.matrix()

J <- data.table::fread("test_J.txt") |> as.matrix()

LL <- -556

test_CLAIC <- CLAIC(J,H,LL)

#to get LL values, as well as parameter estimates, pulling the converged optimization results from the Godambe.ci file and then taking the median values

bottlegrowth_converged <- data.table::fread("bottlegrowth_converged.txt")

bottlegrowth_medians <- bottlegrowth_converged |> dplyr::summarize(LL = median(`Log(likelihood)`), nuB = median(nuB), nuF = median(nuF), T = median(T), theta = median(theta))

#pulling middle H and J matrices, which should correspond to the 95% CI with a 0.01 step size - have to reformat to gget them to play nice with R

bottlegrowth_H <- data.table::fread("bottlegrowth_H.txt") |> as.matrix()

bottlegrowth_J <- data.table::fread("bottlegrowth_J.txt") |> as.matrix()

bottlegrowth_LL <- bottlegrowth_medians$LL

bottlegrowth_CLAIC <- CLAIC(bottlegrowth_H, bottlegrowth_J, bottlegrowth_LL)

growth_converged <- data.table::fread("growth_converged.txt")

growth_medians <- growth_converged |> dplyr::summarize(LL = median(`Log(likelihood)`), nu = median(nu), T = median(T), theta = median(theta))

growth_H <- data.table::fread("growth_H.txt") |> as.matrix()

growth_J <- data.table::fread("growth_J.txt") |> as.matrix()

growth_LL <- growth_medians$LL

growth_CLAIC <- CLAIC(growth_H, growth_J, growth_LL)

snm <- data.table::fread("snm_1d.demog.params.InferDM.bestfits", skip = 1)

names(snm) <- c("LL", "theta", "V3")

snm_H <- 0.00022408 |> as.matrix()

snm_J <- 23.32883465 |> as.matrix()

snm_LL <- snm$LL

snm_CLAIC <- CLAIC(snm_H, snm_J, snm_LL)

two_epoch_converged <- data.table::fread("two_epoch_converged.txt")

two_epoch_medians <- two_epoch_converged |> dplyr::summarize(LL = median(`Log(likelihood)`), nu = median(nu), T = median(T), theta = median(theta))

two_epoch_H <- data.table::fread("two_epoch_H.txt") |> as.matrix()

two_epoch_J <- data.table::fread("two_epoch_J.txt") |> as.matrix()

two_epoch_LL <- two_epoch_medians$LL

two_epoch_CLAIC <- CLAIC(two_epoch_H, two_epoch_J, two_epoch_LL)

#three epoch model did not converge after 1000 optimizations but did after 10000

three_epoch_converged <- data.table::fread("three_epoch_converged.txt")

three_epoch_medians <- three_epoch_converged |> dplyr::summarize(LL = median(`Log(likelihood)`), nuB = median(nuB), nuF = median(nuF), TB = median(TB), TF = median(TF), theta = median(theta))

three_epoch_H <- data.table::fread("three_epoch_H.txt") |> as.matrix()

three_epoch_J <- data.table::fread("three_epoch_J.txt") |> as.matrix()

three_epoch_LL <- three_epoch_medians$LL

three_epoch_CLAIC <- CLAIC(three_epoch_H, three_epoch_J, three_epoch_LL)
