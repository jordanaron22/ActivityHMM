##
## Make function nestedloop available to generate data for
## nested loop plot
##
source("C:/Users/jorda/Documents/Grad School/Research/ActivityHMM/NestedLoop/nestedloop.r")


##
## Load dataset 'res' with simulation results
##
load("C:/Users/jorda/Documents/Grad School/Research/ActivityHMM/NestedLoop/res.rda")
##
res$theta <- round(res$theta, 2)
res$rho2  <- round(res$rho^2, 2)
##
## Reorder dataset with simulation results according to
## simulation parameters theta, rho2, p.c, tau2, k.
##
## Simulation parameters theta, p.c and k are sorted in
## decreasing order (see argument sign).
##
nldata <- nestedloop(res,
                     varnames=c("theta", "rho2", "p.c", "tau2", "k"),
                     varlabels=
                     c("Odds ratio", "Selection",
                       "Control event proportion", "Heterogeneity",
                       "Number of studies"),
                     sign=c(-1, 1, -1, 1, -1))


##
##
## Generate Figure 2
##
##

##
## R object pd2 is used in nested-loop plot (Figure 2)
##
pd2 <- nldata
##
## Use labels instead of numeric values for rho2 and p.c
##
pd2$rho2 <- factor(pd2$rho2,
                   levels=c(0, 0.36, 0.64, 1),
                   labels=c("no", "weak", "moderate", "strong"))
pd2$p.c <- factor(pd2$p.c,
                  levels=c(0.05, 0.1, 0.2, 0.3),
                  labels=c("5%", "10%", "20%", "30%"))


pdf("Figure2.pdf", paper="a4r", width=18, height=15)
##
par(pty="m")
##
## Create skeleton of nested-loop plot using standard R plot function
##
plot(pd2$theta,
     log=c("y"), type="n",
     ylim=c(0.1, 1.1), bty="n",
     xlab="4 x 4 x 4 x 4 x 3 = 768 ordered scenarios",
     ylab="Odds ratio",
     las=1, xaxt="n")
##
## Add vertical lines (using R function lines.nestedloop)
##
lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
##
## Add reference lines (using R function lines.nestedloop)
##
lines(pd2, which="r",
      ymin.refline=0.09, ymax.refline=0.20,
      cex.ref=0.7)
##
## Estimates and legend (using standard R functions lines and legend)
##
lines(pd2$theta, col="black", lwd=2, type="s")   # True theta
lines(pd2$exp.Peto, col="darkgray", type="s")    # Peto
lines(pd2$exp.Trimfill, col="green", type="s")   # Trimfill
lines(pd2$exp.Peters, col="blue", type="s")      # Peters
lines(pd2$exp.G2, col="red", type="s")           # Limit meta-analysis, method 1
lines(pd2$exp.LimF, col="violet", type="s")      # Limit meta-analysis, method 2
lines(pd2$exp.expect, col="gold", type="s")      # Limit meta-analysis, method 3
##
legend(1, 0.55,
       lwd=c(2, rep(1, 6)),
       col=c("black", "darkgray", "green",
         "blue", "red", "violet", "gold"),
       cex=0.9,
       bty="n",
       c("True theta", "Peto method",
         "Trim and Fill method", "Peters method",
         "Limit meta-analysis, method 1",
         "Limit meta-analysis, method 2",
         "Limit meta-analysis, method 3"))
##
dev.off()


##
##
## Generate Figure 3
##
##

##
## R object pd3 is used in hybrid plot (Figure 3)
##
pd3 <- nldata
##
## Use labels instead of numeric values for theta and rho2
##
pd3$theta <- factor(pd3$theta,
                    levels=c(0.5, 0.67, 0.75, 1),
                    labels=c("OR=0.5", "OR=0.67", "OR=0.75", "OR=1"))
pd3$rho2 <- factor(pd3$rho2,
                   levels=c(0, 0.36, 0.64, 1),
                   labels=c("No selection", "Weak selection",
                     "Moderate selection", "Strong selection"))
##
## Generate new variable 'cond' with information for subpanels.
##
## First two criteria are theta (4 categories) and rho (4 categories),
## like in the lattice plot, where both these were used to
## characterise the 4x4 panels.
##
pd3$cond <- interaction(pd3[,attr(pd3, "varnames")[2]],
                        pd3[,attr(pd3, "varnames")[1]],
                        sep=", ")
layout1 <- length(unique(pd3[,attr(pd3,
                                   "varnames")[1]]))
layout2 <- length(unique(pd3[,attr(pd3,
                                   "varnames")[2]]))
##
## Generate new variable 'index' indicating different scenarios within
## subpanel.
##
## For each subpanel, i.e. combination of theta and rho, we have
## 4*4*3 = 48 scenarios (p.c, tau2, k).
##
pd3$index <- seq(1:unique(table(pd3$cond)))
##
varcommon <- c("index", "cond", "theta", "rho2", "p.c", "tau2", "k")
##
method <- c("Peto method",
            "Trim and Fill method",
            "Peters method",
            "Limit meta-analysis, method 1",
            "Limit meta-analysis, method 2",
            "Limit meta-analysis, method 3")
##
pd3.Peto     <- pd3[, c(varcommon, "bias.Peto")]
pd3.Trimfill <- pd3[, c(varcommon, "bias.Trimfill")]
pd3.Peters   <- pd3[, c(varcommon, "bias.Peters")]
pd3.G2       <- pd3[, c(varcommon, "bias.G2")]
pd3.LimF     <- pd3[, c(varcommon, "bias.LimF")]
pd3.expect   <- pd3[, c(varcommon, "bias.expect")]
##
names(pd3.Peto) <- names(pd3.Trimfill) <-
  names(pd3.Peters) <- names(pd3.G2) <-
  names(pd3.LimF) <- names(pd3.expect) <-
  c(varcommon, "bias")
##
pd3.Peto$method     <- method[1]
pd3.Trimfill$method <- method[2]
pd3.Peters$method   <- method[3]
pd3.G2$method       <- method[4]  # beta-lim = LMA method 1
pd3.LimF$method     <- method[5]  # mu-lim   = LMA method 2
pd3.expect$method   <- method[6]  # beta-0   = LMA method 3
##
pd3 <- rbind(pd3.Peto,
             pd3.Trimfill,
             pd3.Peters,
             pd3.G2,
             pd3.LimF,
             pd3.expect)
##
pd3$method <- factor(pd3$method, levels=method)


## Plot figure (hybrid plot)
##
library(lattice)
##

pdf("Figure3.pdf", paper="a4r", width=18, height=15)
##
col.method <- c("gray", "green", "blue", "red", "violet", "gold")
lty.method <- rep(1,6)
##
xyplot(bias ~ index | cond,
       data=pd3,
       groups=method,
       type="s",
       ylim=c(-1.25, 0.25),
       xlab="Simulation scenarios",
#       xlab=expression(paste("4 x 4 x 3 = 48 scenarios per layer: event proportion (5%, 10%, 20%, 30%); heterogeneity (",tau^2," = 0, 0.05, 0.1, 0.2); number of studies (5, 10, 20)")),
       ylab="Bias",
       col=col.method,
       lty=lty.method,
       scales=list(y=list(at=c(0.2, 0, -0.2, -0.4))),
       panel=function(x,y,...){
         ##panel.grid(h=-1, v=-1)
         panel.abline(v=1:48, col="lightgrey")
         panel.abline(h=0)
         ##panel.lines(x=1:48, y=nldata$p.c[1:48]-1, type="s")
         panel.nestedloop(pd2[1:48,],
                          varnames=c("p.c", "tau2", "k"),
                          varlabels=
                          c("Control event proportion", "Heterogeneity",
                            "Number of studies"),
                          which="r",
                          ymin.refline=-1.2, ymax.refline=-0.60,
                          log=FALSE,
                          cex.ref=0.35)
         panel.superpose(x,y, ...)
       },
       par.strip.text=list(cex=0.7),
       key=list(lines=list(col=col.method,lty=lty.method),
                text=list(method), space="right"),
       layout=c(layout2, layout1))
##
dev.off()


pdf("Figure3-old.pdf", width=18, height=15)
##
col.method <- c("gray", "green", "blue", "red", "violet", "gold")
lty.method <- rep(1,6)
##
xyplot(bias ~ index | cond,
       data=pd3,
       groups=method,
       type="s",
       xlab=expression(paste("4 x 4 x 3 = 48 scenarios per layer: event proportion (5%, 10%, 20%, 30%); heterogeneity (",tau^2," = 0, 0.05, 0.1, 0.2); number of studies (5, 10, 20)")),
       ylab="Bias",
       col=col.method,
       lty=lty.method,
       panel=function(x,y,...){
         ##panel.grid(h=-1, v=-1)
         ##panel.abline(v=1:48, col="lightgrey")
         panel.abline(h=0)
         panel.superpose(x,y, ...)
       },
       key=list(lines=list(col=col.method,lty=lty.method),
                text=list(method)),
       layout=c(layout2, layout1))
##
dev.off()
