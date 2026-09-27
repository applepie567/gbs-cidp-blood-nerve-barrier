# Independent base-R check of the conditional linear models and correlations.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else "."
d <- read.csv(file.path(root,"results/external_donor_scores.csv"))
y <- d$Fc_score
group <- as.numeric(d$diagnosis == "CIDP")
index <- as.numeric(scale(d$macrophage_index))
model <- lm(y ~ group + index)
X <- model.matrix(model)
u <- residuals(model)/(1-hatvalues(model))
B <- solve(crossprod(X))
V <- B %*% crossprod(X, X * as.numeric(u^2)) %*% B
se <- sqrt(V[2,2])
estimate <- coef(model)[2]
df <- df.residual(model)
CI <- estimate + c(-1,1)*qt(.975,df)*se
p <- 2*pt(-abs(estimate/se),df)
truth <- read.csv(file.path(root,"results/external_conditional_models.csv"))
t <- truth[truth$model == "Macrophage index adjusted",]
stopifnot(abs(t$difference-estimate)<1e-10, max(abs(c(t$ci_low,t$ci_high)-CI))<1e-9, abs(t$model_p-p)<1e-10)
cat("Conditional model coefficient",estimate,"CI",CI,"P",p,"df",df,"\n")
for(scheme in c("donor","centre")) {
 x <- read.csv(file.path(root,paste0("results/",scheme,"_holdout_predictions.csv")))
 cat(scheme,"nucleus fraction Spearman",cor(x$index,x$nucleus_fraction,method="spearman"),"\n")
 cat(scheme,"captured count fraction Spearman",cor(x$index,x$captured_count_fraction,method="spearman"),"\n")
}
cat("All conditional coefficient, interval and P-value checks passed.\n")
