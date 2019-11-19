#' fit_model
#'
#' an internal function for HarrellPlot
#' @import emmeans
#' @import lme4
#' @import lmerTest
#' @import data.table
#' @export
fit_model <- function(
  x,
  y,
  g,
  covcols=NULL,
  rintcols=NULL,
  rslopecols=NULL,
  dt,
  fit.model='lm', # lm, lmm, glm
  glm_family='gaussian', # family for glm()
  REML=TRUE, # if FALSE then fit by ML
  add_interaction=TRUE,
  interaction.group=TRUE,
  interaction.treatment=TRUE,
  mean_intervals.method='raw', # model for CI of mean
  conf.mean=0.95, # confidence level for CI of mean
  contrasts.method='revpairwise', # which contrasts to show
  contrasts.scaling='raw',
  conf.contrast=0.95,
  adjust=FALSE
){
  
  if(g=='dummy_g'){
    xcols <- x
    grouping <- FALSE
    add_interaction <- FALSE
  }else{
    xcols <- c(x,g)
    grouping <- TRUE
  }

  if(add_interaction==TRUE){
    icols <- c(x,g)
  }else{
    icols <- NULL
  }

  # make model formula and fit model
  model_formula <- formula(make_formula_str(y, xcols, rintcols, rslopecols, icols, covcols))
  if(fit.model=='lm'){
    fit <- lm(model_formula, data=dt)
  }
  if(fit.model=='lmm'){
    fit <- lmerTest::lmer(model_formula, data=dt, REML=REML)
  }
  if(fit.model=='glm'){
    if(glm_family=="nb") # negative binomial
      fit <- glm.nb(model_formula, data=dt)
    if(glm_family=="poisson")
      fit <- glm(model_formula, family="poisson", data=dt) # log link
  }

  emm <- emmeans(fit, specs=xcols)
  if(fit.model=="glm"){
    emm.link <- copy(emm)
    emm <- emmeans(fit, specs=xcols, type="response")
  }
  
  # save global
  tables <- list(NULL)
  tables$fit <- fit
  tables$form_str <- model_formula
  tables$coeffs <- coefficients(summary(fit))
  tables$summary <- broom::glance(fit)
  tables$summary.raw <- summary(fit)
  tables$contrasts.scaled <- data.table(NULL)
  # grudgingly include anova tables
  if(fit.model=='lm'){
    tables$anova.1 <- anova(fit)
    tables$anova.2 <- car::Anova(fit, type='II')
    options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))
    tables$anova.3 <- car::Anova(lm(model_formula, data=dt), type='III')
    options(contrasts=c(unordered="contr.treatment", ordered="contr.poly"))
  }
  if(fit.model=='lmm'){
    # tables$anova.3 <- Anova(lmer(model_formula, data=dt), type='III')
    tables$anova.1 <- anova(fit, type=1)
    tables$anova.2 <- anova(fit, type=2)
    tables$anova.3 <- anova(fit, type=3)
  }
  if(fit.model=='glm'){
    tables$anova.1 <- anova(fit)
    tables$anova.2 <- car::Anova(fit, type='II')
    options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))
    if(glm_family=="nb") # negative binomial
      refit <- glm.nb(model_formula, data=dt)
    if(glm_family=="poisson")
      refit <- glm(model_formula, family="poisson", data=dt) # log link
    tables$anova.3 <- car::Anova(refit, type='III')
    options(contrasts=c(unordered="contr.treatment", ordered="contr.poly"))
  }
  

  #     Bayes model
  #  mad <- median(abs(dt[,y] - mean(dt[,y])))
  #  #dt[,y.mad:=y/mad]
  #  y.mad <- dt[,y]/mad
  #  fit.mcmc <- MCMCregress(y.mad ~ x, data=dt, b0 = 0, B0 = 0.1, c0=2, d0=0.11)
  #  post.lsm <- emmeans(fit.mcmc, specs='x')
  # # dt <- dt[, .SD, .SDcols=c('x','y')] #drop y.mad because need to bind ci later

  #### Compute the table of estimated marginal means ci_means ##########
  if(mean_intervals.method=='lm'){
    tables$means <- confint(emm, level=conf.mean)
    ci_means <- data.table(tables$means) # mean intervals not adjusted
    if(fit.model=="lm"|fit.model=="lmm"){
      ci_means <- ci_means[, .SD, .SDcols=c(xcols,'emmean','lower.CL','upper.CL')]
    }
    if(fit.model=="glm"){
      ci_means <- ci_means[, .SD, .SDcols=c(xcols,'response','asymp.LCL','asymp.UCL')]
    }
  }
  if(mean_intervals.method=='raw'){
    conf.tail <- conf.mean + (1-conf.mean)/2
    tables$means <- dt[, .(
      mean=mean(get(y)),
      sem=sd(get(y))/sqrt(.N),
      lower=mean(get(y))-sd(get(y))/sqrt(.N)*qt(conf.tail,(.N-1)),
      upper=mean(get(y))+sd(get(y))/sqrt(.N)*qt(conf.tail,(.N-1))),
      by=xcols]
    ci_means <- tables$means[, .SD, .SDcols=c(xcols,'mean', 'lower','upper')]
  }
  tables$means.raw <- tables$means # these will be the same
  if(mean_intervals.method=='boot'){
    dt_boot <- data.table(dt[, Hmisc::smean.cl.boot(get(y),conf.int=conf.mean), by=xcols])
    dt_boot[, tile:=c('a','lower','upper')]
    form <- formula(paste(paste(xcols,collapse='+'),'tile',sep='~'))
    ci_means <- dcast(dt_boot, form, value.var='V1') #**** change x+g to formula
  }
  # if(mean_intervals.method=='bayes'){
  #   conf.tail <- conf.mean + (1-conf.mean)/2
  #   res <- summary(as.mcmc(post.lsm), quantiles = c(0.5, (1-conf.tail), conf.tail))$quantiles*mad
  #   ci_means <- data.table(x=row.names(res),res)
  #   ci_means[, (x):=factor(substr(x,3,nchar(x)))]
  # }
  if(grouping==FALSE){
    ci_means[, (g):='dummy']
    setnames(ci_means, old=colnames(ci_means), new=c(x, y, 'lower', 'upper', g))
    ci_means <- ci_means[,.SD, .SDcols=c(x, g, y,'lower','upper')]
  }else{
    setnames(ci_means, old=colnames(ci_means), new=c(xcols, y, 'lower','upper'))
  }

  #######################################################
  #### Compute the table of contrasts: ci_diffs ##########
  #######################################################
  
  x_levels <- levels(dt[, get(x)]) # levels for means plot
  g_levels <- levels(dt[, get(g)])
  n_levels <- length(x_levels)
  n_groups <- length(g_levels)

  if(contrasts.method=='coefficients'){
    ci_diffs <- coefficients(summary(fit))

    # get rid of df column from lmer fit
    if(fit.model=='lmm'){
      ci_diffs <- ci_diffs[,-which(colnames(ci_diffs)=='df')]
    }

    x_names <- c(x_levels[-1], g_levels[-1])
    if(add_interaction==TRUE){
      temp <- expand.grid(x_levels[-1], g_levels[-1])
      x_names <- c(x_names, paste(temp$Var1, temp$Var2, sep=':'))
    }
    # confint
    ci_ci <- confint(fit, level=conf.contrast)[row.names(ci_diffs),]
    ci_diffs <- cbind(ci_diffs, ci_ci)
    tables$contrasts <- data.table(contrast=x_names, ci_diffs[-1,])
    ci_diffs <- data.table(contrast=x_names, g='dummy', ci_diffs[-1,])
    setnames(ci_diffs, old=colnames(ci_diffs), new = c('contrast', 'g', 'estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'lower', 'upper'))
    ci_diffs <- ci_diffs[, .SD, .SDcols=c('contrast','g','estimate','lower','upper')]
    ci_diffs[, contrast:=factor(contrast, levels=x_names)]
    ci_diffs<-ci_diffs[nrow(ci_diffs):1,] # reorder to match table output
  }
  if(contrasts.method!='coefficients'){
    ci.adjust <- 'none'
    if(adjust==TRUE){
      ci.adjust <- ifelse(contrasts.method=='trt.vs.ctrl1', 'dunnettx','tukey')
    }
    if(grouping==FALSE){
      ci_diffs <- summary(contrast(emm,
                          method=contrasts.method,
                          adjust=ci.adjust,
                          level=conf.contrast),
                          infer=c(TRUE,TRUE))
      tables$contrasts.raw <- ci_diffs
      ci_diffs <- data.table(ci_diffs)
    }
    if(grouping==TRUE){
      if(contrasts.method=="trt.vs.ctrl1"){ # similar to coefficients but interaction is a contrast
        ci_diffs <- summary(contrast(emm,
                                     method=contrasts.method,
                                     adjust=ci.adjust,
                                     level=conf.contrast),
                            infer=c(TRUE,TRUE))
      }
      if(contrasts.method=="revpairwise" & add_interaction==FALSE){
        ci_diffs.p1 <- summary(contrast(
          emmeans(fit, specs=xcols[1], type="response"),
          method=contrasts.method,
          adjust=ci.adjust,
          level=conf.contrast), 
          infer=c(TRUE,TRUE))
        ci_diffs.p2 <- summary(contrast(
          emmeans(fit, specs=xcols[2], type="response"),
          method=contrasts.method,
          adjust=ci.adjust,
          level=conf.contrast), 
          infer=c(TRUE,TRUE))
        ci_diffs <- rbind(ci_diffs.p1, ci_diffs.p2)
      }
      if(contrasts.method=="revpairwise" & add_interaction==TRUE){ # subset into pairwise within each group
        ci_diffs <- summary(contrast(emm,
                                     method=contrasts.method,
                                     adjust=ci.adjust,
                                     level=conf.contrast,
                                     simple = "each",
                                     combine=TRUE),
                            infer=c(TRUE,TRUE))
      }
      tables$contrasts.raw <- ci_diffs
      ci_diffs <- data.table(ci_diffs)
      # if(fit.model=="glm"){
      #   setnames(ci_diffs,
      #            old=c("ratio", "asymp.LCL", "asymp.UCL"),
      #            new=c("estimate", "lower.CL", "upper.CL")
      #            )
      # }
      
      if(add_interaction==FALSE){ # delete redundant rows
        # this is inelegant
        # ci_diffs <- ci_diffs[, .(
        #   estimate=mean(estimate),
        #   SE=mean(SE),
        #   df=mean(df),
        #   lower.CL=mean(lower.CL),
        #   upper.CL=mean(upper.CL),
        #   t.ratio=mean(t.ratio),
        #   p.value=mean(p.value)
        # ), by=contrast]
      }
      if(add_interaction==TRUE & contrasts.method=='revpairwise'){
        ci_diffs[, contrast:=ifelse(ci_diffs[, get(g)]!=".",
                      paste(ci_diffs[, get(g)], contrast, sep=":"),
                      paste(ci_diffs[, get(x)], contrast, sep=":"))]
        ci_diffs.g <- data.table(NULL)
        ci_diffs.x <- data.table(NULL)
        if(interaction.group==TRUE){
          ci_diffs.g <- ci_diffs[get(g)!=".",]
        }
        if(interaction.treatment==TRUE){
          ci_diffs.x <- ci_diffs[get(x)!=".",]
        }
        ci_diffs <- rbind(ci_diffs.g, ci_diffs.x)
      }
    }
    ycols <- setdiff(colnames(ci_diffs), c(x,g))
    tables$contrasts <- ci_diffs[, .SD, .SDcols=ycols]
    if(fit.model=="glm"){
      setnames(ci_diffs,
               old=c("ratio", "asymp.LCL", "asymp.UCL"),
               new=c("estimate", "lower.CL", "upper.CL")
      )
    }
    ci_diffs[, g:='dummy'] # why do I have this again?
    ci_diffs <- ci_diffs[, .SD, .SDcols=c('contrast','g','estimate','lower.CL','upper.CL')]
  }
  
  if(fit.model=='bayes'){
    conf.tail <- conf.contrast + (1-conf.contrast)/2
    res <- summary(as.mcmc(contrast(post.lsm, method=contrasts.method)), quantiles = c(0.5, (1-conf.tail), conf.tail))$quantiles*mad
    ci_diffs <- data.table(x=row.names(res),res)
    ci_diffs[, x:=factor(substr(x,10,nchar(x)))]
  }

  # make sure factor order of ci_diffs is in order they appear in the table
  ci_diffs[, contrast:=factor(contrast, ci_diffs$contrast)]

  ci_diffs.raw <- copy(ci_diffs)
  
  yscale <- 1 # default
  if(contrasts.scaling=='standardized'){
    yscale <- summary(fit)$sigma
    ci_diffs[, estimate:=estimate/yscale]
    ci_diffs[, lower.CL:=lower.CL/yscale]
    ci_diffs[, upper.CL:=upper.CL/yscale]
    ycols <- setdiff(colnames(ci_diffs), "g")
    tables$contrasts.scaled <- data.table(ci_diffs[, .SD, .SDcols=ycols])
  }
  if(contrasts.scaling=='percent' & fit.model!="glm"){
    # SE and CI and t and p are not simple transformations like mean
    # SE is computed using the delta method following
    # https://www2.census.gov/programs-surveys/acs/tech_docs/accuracy/percchg.pdf
    # which has excellent coverage given my small simulation, which results are here
    # https://rdoodles.rbind.io/notebooks/relative_standard_errors.nb.html
    scale.o <- ci_diffs[1, estimate]
    if(grouping==FALSE){
      working <- ci_diffs
      split1 <- data.frame(t(do.call("cbind", strsplit(as.character(working$contrast)," - "))))
      ref <- as.character(split1[, "X2"])
      nonref <- as.character(split1[, "X1"])
      cell_id <- ci_means[, get(x)]
      match_it_ref <- match(ref, cell_id)
      denom <- ci_means[match_it_ref, get(y)]
      se.denom <- summary(emm)[match_it_ref, "SE"]
      
      match_it_nonref <- match(nonref, cell_id)
      num <- ci_means[match_it_nonref, get(y)]
      se.num <- summary(emm)[match_it_nonref, "SE"]
    }
    if(grouping==TRUE & add_interaction==TRUE & contrasts.method=='revpairwise'){
      working <- data.table(tables$contrasts.raw)
      split1 <- data.frame(t(do.call("cbind", strsplit(as.character(working$contrast)," - "))))
      ref.p1 <- ifelse(working[, get(g)]!=".", working[, get(g)], as.character(split1[, "X2"]))
      ref.p2 <- ifelse(working[, get(x)]!=".", working[, get(x)], as.character(split1[, "X2"]))
      ref <- paste(ref.p1, ref.p2, sep="-")
      cell_id <- paste(ci_means[, get(g)], ci_means[, get(x)], sep="-")
      match_it_ref <- match(ref, cell_id)
      denom <- ci_means[match_it_ref, get(y)]
      se.denom <- summary(emm)[match_it_ref, "SE"]
      
      nonref.p1 <- ifelse(working[, get(g)]!=".", working[, get(g)], as.character(split1[, "X1"]))
      nonref.p2 <- ifelse(working[, get(x)]!=".", working[, get(x)], as.character(split1[, "X1"]))
      nonref <- paste(nonref.p1, nonref.p2, sep="-")
      match_it_nonref <- match(nonref, cell_id)
      num <- ci_means[match_it_nonref, get(y)]
      se.num <- summary(emm)[match_it_nonref, "SE"]
    }
    if(grouping==TRUE & add_interaction==TRUE & contrasts.method=='trt.vs.ctrl1'){
      working <- data.table(tables$contrasts.raw)
      split1 <- data.frame(t(do.call("cbind", strsplit(as.character(working$contrast)," - "))))
      split2 <- data.frame(t(do.call("cbind", strsplit(as.character(split1$X2),","))))
      ref <- paste(split2[, "X2"], split2[, "X1"], sep="-")
      cell_id <- paste(ci_means[, get(g)], ci_means[, get(x)], sep="-")
      match_it_ref <- match(ref, cell_id)
      denom <- ci_means[match_it_ref, get(y)]
      se.denom <- summary(emm)[match_it_ref, "SE"]
      
      split3 <- data.frame(t(do.call("cbind", strsplit(as.character(split1$X1),","))))
      nonref <- paste(split3[, "X2"], split3[, "X1"], sep="-")
      match_it_nonref <- match(nonref, cell_id)
      denom <- ci_means[match_it_nonref, get(y)]
      se.num <- summary(emm)[match_it_nonref, "SE"]
    }
    if(grouping==TRUE & add_interaction==FALSE){
      marginal_means <- rbind(emmeans(fit, specs=g, adjust=ci.adjust),
                              emmeans(fit, specs=x, adjust=ci.adjust))
      marginal_means <- data.table(summary(marginal_means))
      marginal_means[, cell_id:=ifelse(get(g)!=".", get(g), get(x))]
      # the denominator will be marginal means
      working <- ci_diffs
      split1 <- data.frame(t(do.call("cbind", strsplit(as.character(working$contrast)," - "))))
      ref <- as.character(split1[, "X2"])
      cell_id <- ifelse(marginal_means[, get(g)!="."], marginal_means[, get(g)], marginal_means[, get(x)])
      match_it_ref <- match(ref, cell_id)
      denom <- marginal_means[match_it_ref, emmean]
      se.denom <- marginal_means[match_it_ref, "SE"]
      
      nonref <- as.character(split1[, "X1"])
      match_it_nonref <- match(nonref, cell_id)
      num <- marginal_means[match_it_nonref, emmean]
      se.num <- marginal_means[match_it_nonref, "SE"]
    }
    ci_diffs[, estimate:=100*estimate/denom]
    df <- tables$contrasts$df
    prob <- conf.contrast + (1-conf.contrast)/2
    tcrit <- qt(prob, df)
    # se of percent difference using delta method
    se <- abs(num/denom)*sqrt((se.num^2/num^2 + se.denom^2/denom^2))
    # # check using deltamethod function for trt vs cntrl
    #   se2[1] <- msm::deltamethod(~x2/x1, mean=coef(fit), cov=vcov(fit))
    #   se2[2] <- msm::deltamethod(~x3/x1, mean=coef(fit), cov=vcov(fit))
    #   se2[3] <- msm::deltamethod(~x4/x1, mean=coef(fit), cov=vcov(fit))
    # }
    ci_diffs[, lower.CL:=estimate - 100*tcrit*se]
    ci_diffs[, upper.CL:=estimate + 100*tcrit*se]

    # # rescale back to scale.o
    # yscale <- scale.o/ci_diffs[1, estimate]
    # ci_diffs[, estimate:=estimate*yscale]
    # ci_diffs[, lower.CL:=lower.CL*yscale]
    # ci_diffs[, upper.CL:=upper.CL*yscale]
    ycols <- setdiff(colnames(ci_diffs), "g")
    tables$contrasts.scaled <- data.table(ci_diffs[, .SD, .SDcols=ycols])
  }

  if(contrasts.scaling=='percent' & fit.model=="glm"){
    ci_diffs[, estimate:=ifelse(estimate>1, 100*(estimate-1), -100*(1-estimate))]
    ci_diffs[, lower.CL:=ifelse(lower.CL>1, 100*(lower.CL-1), -100*(1-lower.CL))]
    ci_diffs[, upper.CL:=ifelse(upper.CL>1, 100*(upper.CL-1), -100*(1-upper.CL))]
    ycols <- setdiff(colnames(ci_diffs), "g")
    tables$contrasts.scaled <- data.table(ci_diffs[, .SD, .SDcols=ycols])
  }
  setnames(ci_diffs, old=colnames(ci_diffs), new=c(x, g, y,'lower','upper'))
  return(list(fit=fit, ci_means=ci_means, ci_diffs=ci_diffs, tables=tables, yscale=yscale))
}
