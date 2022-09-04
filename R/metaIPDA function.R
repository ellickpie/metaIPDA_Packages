#' meta-analysis by IPDA
#'
#' add later
#' @param dir.location The directory location that contains individual studies data in csv format
#' @param DV the name of DV, must be in character format, e.g., "DV01"
#' @return kable outputs of IPDA meta analysis results
#' @examples 
#' output <- metaIPDA 
#' @import lmerTest rvest
#' @export
metaIPDA = function(
    dir.location,
    DV,
    DV.nature,
    IV,
    IV.nature,
    Moderator = "",
    Moderator.nature = "",
    Control = ""#must be numeric
)
  
{
  # preloading####
  
  {
    if (!require(devtools)) {install.packages('rvest')}
    if (!require(devtools)) {install.packages('MASS')}
    if (!require(devtools)) {install.packages('dplyr')}
    if (!require(devtools)) {install.packages('stringr')}
    if (!require(devtools)) {install.packages('lme4')}
    if (!require(devtools)) {install.packages('lmerTest')}
    if (!require(devtools)) {install.packages('nlme')}
    if (!require(devtools)) {install.packages('jtools')}
    if (!require(devtools)) {install.packages('DescTools')}
    if (!require(devtools)) {install.packages('knitr')}
    if (!require(devtools)) {install.packages('kableExtra')}
    
    options("scipen" = 100, "digits" = 4)
    
    
  }
  
  #combining data ready for IPDA analysis below: data.analysis ####
  #source("data combining.R")
  {
    study = list()
    for (i in 1:length(list.files(dir.location)))
    {
      filename = paste(dir.location, "/", list.files(dir.location)[i], sep = "")
      study[[i]] =
        read.csv(filename)
      study[[i]]$study.cluster =
        paste("study.cluster", i, sep = "") %>% as.factor()
      
      if (DV.nature == "nominal")
      {
        study[[i]]$DV = as.factor(study[[i]][DV][, 1])
      }
      if (DV.nature == "continuous")
      {
        study[[i]]$DV = as.numeric(scale(as.numeric(study[[i]][DV][, 1])))
      }
      
      
      if (IV.nature == "nominal")
      {
        study[[i]]$IV = as.factor(study[[i]][IV][, 1])
      }
      if (IV.nature == "continuous")
      {
        study[[i]]$IV = as.numeric(scale(as.numeric(study[[i]][IV][, 1])))
      }
      
      if (Moderator != "")
      {
        if (Moderator.nature == "nominal")
        {
          study[[i]]$Moderator = as.factor(study[[i]][Moderator][, 1])
        }
        if (Moderator.nature == "continuous")
        {
          study[[i]]$Moderator = as.numeric(scale(as.numeric(study[[i]][Moderator][, 1])))
        }
      }
      
      if (Control[1] != "")
        
      {
        for (ii in 1:length(Control))
        {
          temptext = paste("study[[i]]$Control",
                           ii,
                           " = study[[i]][Control[ii]][,1]",
                           sep = "")
          eval(parse(text = temptext))
          
          if (Control.nature[ii] == "continuous")
          {
            temptext = paste(
              "study[[i]]$Control",
              ii,
              " = as.numeric(scale(as.numeric(study[[i]]$Control",
              ii,
              ")))",
              sep = ""
            )
            eval(parse(text = temptext))
          }
          
          if (Control.nature[ii] == "nominal")
          {
            temptext =
              paste("study[[i]]$Control",
                    ii,
                    " = as.factor(study[[i]]$Control",
                    ii,
                    ")",
                    sep = "")
            eval(parse(text = temptext))
          }
          
          
        }
        
        
      }
      
      
    }
    
    data.analysis = do.call(rbind, study)
  }
  
  #set contrasts for all factors below ####
  #source("contrasts setting.R")
  {
    k = nlevels(data.analysis$study.cluster)
    contrasts(data.analysis$study.cluster) = contr.treatment(k) - (1 / k)
    
    if (IV.nature == "nominal")
    {
      k = nlevels(data.analysis$IV)
      contrasts(data.analysis$IV) = contr.treatment(k) - (1 / k)
    }
    
    if (Moderator.nature == "nominal")
    {
      k = nlevels(data.analysis$Moderator)
      contrasts(data.analysis$Moderator) = contr.treatment(k) - (1 / k)
    }
    
    
    if (Control[1] != "")
      
    {
      for (ii in 1:length(Control))
      {
        if (Control.nature[ii] == "nominal")
        {
          temptext = paste("Control", ii, sep = "")
          k = nlevels((data.analysis[temptext])[, 1])
          contrast_text = paste("contrasts(data.analysis$",
                                temptext,
                                ") = contr.treatment(k) - (1 / k)",
                                sep = "")
          eval(parse(text = contrast_text))
          
        }
        
        
        
      }
    }
    
  }
  
  ##creating glm text ####
  #source("glm text.R")
  {
    if (Moderator == "") {
      glmtext_Moderator = ""
    }
    if (Moderator != "")
    {
      glmtext_Moderator = paste("Moderator", collapse = " * ")
      glmtext_Moderator = paste(" * ", glmtext_Moderator, sep = "")
      
      glmtext_Moderator.plus1SD = paste("(Moderator - 1)", collapse = " * ")
      glmtext_Moderator.plus1SD = paste(" * ", glmtext_Moderator.plus1SD, sep = "")
      
      glmtext_Moderator.minus1SD = paste("(Moderator + 1)", collapse = " * ")
      glmtext_Moderator.minus1SD = paste(" * ", glmtext_Moderator.minus1SD, sep = "")
      
    }
    
    
    
    if (Control[1] == "")
    {
      glmtext_Control = ""
    }
    
    if (Control[1] != "")
    {
      tempcontrol = list()
      for (ii in 1:length(Control))
      {
        tempcontrol[[ii]] = paste("Control", ii, sep = "")
      }
      
      glmtext_Control = paste(as.character(tempcontrol), collapse = " + ")
      glmtext_Control = paste(" + ", glmtext_Control, sep = "")
    }
    if (DV.nature == "continuous")
    {
      glmtext.end = ", data = data.analysis))"
    }
    
    if (DV.nature == "nominal")
    {
      glmtext.end = ", data = data.analysis, family = binomial(link = 'logit')))"
    }
    
    
    
    glmtext.head.fixed = "summary(glm(DV ~ IV"
    if (DV.nature == "continuous")
    {
      glmtext.head.random = "summary(lmer(DV ~ IV"
    }
    if (DV.nature == "nominal")
    {
      glmtext.head.random = "summary(glmer(DV ~ IV"
    }
    
  }
  
  ##main effect of IV clean Model####
  #source("analysis main effect clean.R")
  {
    glmtext.clean.fixed =
      paste(glmtext.head.fixed,
            " * study.cluster",
            glmtext.end,
            sep = "")
    
    glmtext.clean.random =
      paste(glmtext.head.random,
            " + (1 + IV | study.cluster)",
            glmtext.end,
            sep = "")
    
    IPDA.clean.fixed =
      eval(parse(text = glmtext.clean.fixed))
    
    IPDA.clean.random =
      eval(parse(text = glmtext.clean.random))
  }
  
  ##main effect of IV full Model ####
  #source("analysis main effect full.R")
  {
    glmtext.full.fixed =
      paste(
        glmtext.head.fixed,
        glmtext_Moderator,
        " * study.cluster",
        glmtext_Control,
        glmtext.end,
        sep = ""
      )
    
    glmtext.full.random =
      paste(
        glmtext.head.random,
        glmtext_Moderator,
        glmtext_Control,
        " + (1 + IV | study.cluster)",
        glmtext.end,
        sep = ""
      )
    
    
    IPDA.full.fixed =
      eval(parse(text = glmtext.full.fixed))
    
    IPDA.full.random =
      eval(parse(text = glmtext.full.random))
  }
  
  ##Simple effects of IV ####
  #source ("analysis simple effect.R")
  {
    if (Moderator != "")
      
    {
      if (Moderator.nature == "nominal")
      {
        mod_levels = nlevels(data.analysis$Moderator)
        
        for (iii in 1:mod_levels)
        {
          mod_contrast = vector()
          mod_contrast[1:mod_levels] = 0
          mod_contrast[iii] = 1
          contrasts(data.analysis$Moderator) = mod_contrast
          var_text = paste("IPDA.full.fixed.", "contrast.", iii, sep = "")
          assign(var_text, eval(parse(text = glmtext.full.fixed)))
          var_text = paste("IPDA.full.random.", "contrast.", iii, sep = "")
          assign(var_text, eval(parse(text = glmtext.full.random)))
        }
        
      }
      
      
      if (Moderator.nature == "continuous")
        
      {
        data.analysis$Moderator.original = data.analysis$Moderator
        data.analysis$Moderator = data.analysis$Moderator.original - 1
        glmtext.full.plus1SD.fixed =
          paste(
            glmtext.head.fixed,
            glmtext_Moderator,
            " * study.cluster",
            glmtext_Control,
            glmtext.end,
            sep = ""
          )
        IPDA.full.plus1SD.fixed = eval(parse(text = glmtext.full.plus1SD.fixed))
        
        data.analysis$Moderator = data.analysis$Moderator.original + 1
        glmtext.full.minus1SD.fixed =
          paste(
            glmtext.head.fixed,
            glmtext_Moderator,
            " * study.cluster",
            glmtext_Control,
            glmtext.end,
            sep = ""
          )
        IPDA.full.minus1SD.fixed = eval(parse(text = glmtext.full.minus1SD.fixed))
        
        data.analysis$Moderator = data.analysis$Moderator.original - 1
        glmtext.full.plus1SD.random =
          paste(
            glmtext.head.random,
            glmtext_Moderator,
            glmtext_Control,
            " + (1 + IV | study.cluster)",
            glmtext.end,
            sep = ""
          )
        IPDA.full.plus1SD.random = eval(parse(text = glmtext.full.plus1SD.random))
        
        data.analysis$Moderator = data.analysis$Moderator.original + 1
        glmtext.full.minus1SD.random =
          paste(
            glmtext.head.random,
            glmtext_Moderator,
            glmtext_Control,
            " + (1 + IV | study.cluster)",
            glmtext.end,
            sep = ""
          )
        IPDA.full.minus1SD.random = eval(parse(text = glmtext.full.minus1SD.random))
        data.analysis$Moderator = data.analysis$Moderator.original
        
      }
    }
  }
  
  ##summary of results  ####
  output_table = data.frame(
    dumnumber = 1:100,
    label = NA,
    es.fixed = NA,
    es.stat.fixed = NA,
    es.df.fixed = NA,
    es.p.fixed = NA,
    es.random = NA,
    es.stat.random = NA,
    es.df.random = NA,
    es.p.random = NA
  )
  
  #text adjustment
  {
    if (IV.nature == "nominal")
    {
      IV_text = "IV2"
    }
    if (IV.nature == "continuous")
    {
      IV_text = "IV"
    }
    
    if (DV.nature == "nominal")
    {
      stat_text = "z value"
      p_text = "Pr(>|z|)"
    }
    if (DV.nature == "continuous")
    {
      stat_text = "t value"
      p_text = "Pr(>|t|)"
    }
  }
  
  #the first row: Clean Model
  {
    i = 1
    output_table$label[i] = "Clean Model"
    output_table$es.fixed[i] =
      round(IPDA.clean.fixed$coefficients[IV_text, "Estimate"], 4)
    output_table$es.stat.fixed[i] =
      paste("t = ",
            round(IPDA.clean.fixed$coefficients[IV_text, stat_text], 4), sep = "")
    output_table$es.df.fixed[i] = IPDA.clean.fixed$df.residual
    output_table$es.p.fixed[i] =
      round(IPDA.clean.fixed$coefficients[IV_text, p_text], 4)
    
    output_table$es.random[i] =
      round(IPDA.clean.random$coefficients[IV_text, "Estimate"], 4)
    
    if (DV.nature == "continuous")
    {
      output_table$es.stat.random[i] =
        paste("t = ",
              round(IPDA.clean.random$coefficients[IV_text, stat_text], 4), sep = "")
      
      output_table$es.df.random[i] =
        round(IPDA.clean.random$coefficients[IV_text, "df"], 4)
    }
    if (DV.nature == "nominal")
    {
      output_table$es.stat.random[i] =
        paste("z = ",
              round(IPDA.clean.random$coefficients[IV_text, stat_text], 4), sep = "")
      output_table$es.df.random[i] = NA
    }
    output_table$es.p.random[i] =
      round(IPDA.clean.random$coefficients[IV_text, p_text], 4)
  }
  
  #@the second row: Full Model
  {
    i = 2
    output_table$label[i] = "Full Model"
    output_table$es.fixed[i] =
      round(IPDA.full.fixed$coefficients[IV_text, "Estimate"], 4)
    output_table$es.stat.fixed[i] =
      paste("t = ",
            round(IPDA.full.fixed$coefficients[IV_text, stat_text], 4), sep = "")
    output_table$es.df.fixed[i] = IPDA.full.fixed$df.residual
    output_table$es.p.fixed[i] =
      round(IPDA.full.fixed$coefficients[IV_text, p_text], 4)
    
    output_table$es.random[i] =
      round(IPDA.full.random$coefficients[IV_text, "Estimate"], 4)
    
    if (DV.nature == "continuous")
    {
      output_table$es.stat.random[i] =
        paste("t = ",
              round(IPDA.full.random$coefficients[IV_text, stat_text], 4), sep = "")
      
      output_table$es.df.random[i] =
        round(IPDA.full.random$coefficients[IV_text, "df"], 4)
    }
    if (DV.nature == "nominal")
    {
      output_table$es.stat.random[i] =
        paste("z = ",
              round(IPDA.full.random$coefficients[IV_text, stat_text], 4), sep = "")
      output_table$es.df.random[i] = NA
    }
    output_table$es.p.random[i] =
      round(IPDA.full.random$coefficients[IV_text, p_text], 4)
    
  }
  
  
  #@the third row: Moderation
  
  output_table$label[3] = "Test of Moderation"
  
  if (Moderator != "")
  {
    if (Moderator.nature == "nominal" & DV.nature == "continuous")
    {
      #source("output moderation nominal DV.con.R")
      {
        ##### Moderation output table when Moderator is nominal
        
        aa = i + 1
        bb = i + nlevels(data.analysis$Moderator) - 1
        
        #### beta of moderators####
        predictor_list = names(IPDA.full.fixed$coefficients[, 1])
        predictOr_filter = (str_detect(predictor_list, "IV"))
        tempfilter = as.data.frame(cbind(predictor_list, predictOr_filter))
        tempfilter = subset(tempfilter, tempfilter$predictOr_filter == "TRUE")
        moderation_list = tempfilter$predictor_list
        moderation_filter =  str_count(moderation_list, pattern = ":")
        tempfilter = as.data.frame(cbind(moderation_list, moderation_filter))
        tempfilter = subset(tempfilter, tempfilter$moderation_filter == "1")
        moderation_list = tempfilter$moderation_list
        moderation_list = moderation_list [1:nlevels(data.analysis$Moderator) -
                                             1]
        
        
        output_table$label[aa:bb] =
          moderation_list [1:nlevels(data.analysis$Moderator) - 1]
        output_table$es.fixed[aa:bb] =
          round(IPDA.full.fixed$coefficients[, "Estimate"][output_table$label[aa:bb]], 4)
        
        output_table$es.stat.fixed[aa:bb] =
          paste("t = ",
                round(IPDA.full.fixed$coefficients[, "t value"][output_table$label[aa:bb]], 4))
        output_table$es.df.fixed[aa:bb] = IPDA.full.fixed$df[2]
        output_table$es.p.fixed[aa:bb] =
          round(IPDA.full.fixed$coefficients[, "Pr(>|t|)"][output_table$label[aa:bb]], 4)
        
        output_table$label[aa:bb] =
          moderation_list [1:nlevels(data.analysis$Moderator) - 1]
        output_table$es.random[aa:bb] =
          round(IPDA.full.random$coefficients[, "Estimate"][output_table$label[aa:bb]], 4)
        output_table$es.stat.random[aa:bb] =
          paste("t = ",
                round(IPDA.full.random$coefficients[, "t value"][output_table$label[aa:bb]], 4))
        output_table$es.df.random[aa:bb] =
          round(IPDA.full.random$coefficients[, "df"][output_table$label[aa:bb]], 4)
        output_table$es.p.random[aa:bb] =
          round(IPDA.full.random$coefficients[, "Pr(>|t|)"][output_table$label[aa:bb]], 4)
        
        
        
        output_table$label[bb + 1] = "Test of Simple Effects"
        
        #### simple effects
        for (i in 1:nlevels(data.analysis$Moderator))
        {
          entry = bb + 1 + i
          output_table$label[entry] = levels(data.analysis$Moderator)[i]
          
          contrast_text =
            paste("round(IPDA.full.fixed.contrast.",
                  i,
                  "$coefficients[2], 4)",
                  sep = "")
          output_table$es.fixed[entry] =
            eval(parse(text = contrast_text))
          
          contrast_text =
            paste('round(IPDA.full.fixed.contrast.',
                  i,
                  '$coefficients[, "t value"][2], 4)',
                  sep = "")
          output_table$es.stat.fixed[entry] =
            paste("t = ", eval(parse(text = contrast_text)), sep = "")
          
          contrast_text =
            paste('round(IPDA.full.fixed.contrast.', i, '$df[2], 4)', sep = "")
          output_table$es.df.fixed[entry] =
            eval(parse(text = contrast_text))
          
          contrast_text =
            paste('round(IPDA.full.fixed.contrast.',
                  i,
                  '$coefficients[, "Pr(>|t|)"][2], 4)',
                  sep = "")
          output_table$es.p.fixed[entry] =
            eval(parse(text = contrast_text))
          
          ####
          contrast_text =
            paste("round(IPDA.full.random.contrast.",
                  i,
                  "$coefficients[2], 4)",
                  sep = "")
          output_table$es.random[entry] =
            eval(parse(text = contrast_text))
          
          contrast_text =
            paste('round(IPDA.full.random.contrast.',
                  i,
                  '$coefficients[, "t value"][2], 4)',
                  sep = "")
          output_table$es.stat.random[entry] =
            paste("t = ", eval(parse(text = contrast_text)), sep = "")
          
          contrast_text =
            paste('round(IPDA.full.random.contrast.',
                  i,
                  '$coefficients[, "df"][2], 4)',
                  sep = "")
          output_table$es.df.random[entry] =
            eval(parse(text = contrast_text))
          
          contrast_text =
            paste('round(IPDA.full.random.contrast.',
                  i,
                  '$coefficients[, "Pr(>|t|)"][2], 4)',
                  sep = "")
          output_table$es.p.random[entry] =
            eval(parse(text = contrast_text))
          
          
          
        }
      }
    }
    
    if (Moderator.nature == "nominal" & DV.nature == "nominal")
    {
      #source("output moderation nominal DV.nom.R")
      {##### Moderation output table when Moderator is nominal
        
        aa = i + 1
        bb = i + nlevels(data.analysis$Moderator) - 1
        
        #### beta of moderators####
        predictor_list = names(IPDA.full.fixed$coefficients[,1])
        predictOr_filter = (str_detect(predictor_list, "IV"))
        tempfilter = as.data.frame(cbind(predictor_list, predictOr_filter))
        tempfilter = subset(tempfilter, tempfilter$predictOr_filter == "TRUE" )
        moderation_list = tempfilter$predictor_list
        moderation_filter =  str_count(moderation_list, pattern = ":")
        tempfilter = as.data.frame(cbind(moderation_list, moderation_filter))
        tempfilter = subset(tempfilter, tempfilter$moderation_filter == "1" )
        moderation_list = tempfilter$moderation_list
        moderation_list = moderation_list [1: nlevels(data.analysis$Moderator)-1]
        
        
        output_table$label[aa:bb] = 
          moderation_list [1: nlevels(data.analysis$Moderator)-1]
        output_table$es.fixed[aa:bb] =
          round(IPDA.full.fixed$coefficients[, "Estimate"][output_table$label[aa:bb]], 4)
        
        output_table$es.stat.fixed[aa:bb] =
          paste("z = ", 
                round(IPDA.full.fixed$coefficients[, "z value"][output_table$label[aa:bb]], 4)
          )
        output_table$es.df.fixed[aa:bb] = IPDA.full.fixed$df[2]
        output_table$es.p.fixed[aa:bb] =                    
          round(IPDA.full.fixed$coefficients[, "Pr(>|z|)"][output_table$label[aa:bb]], 4)
        
        output_table$label[aa:bb] = 
          moderation_list [1: nlevels(data.analysis$Moderator)-1]
        output_table$es.random[aa:bb] =
          round(IPDA.full.random$coefficients[, "Estimate"][output_table$label[aa:bb]], 4)
        output_table$es.stat.random[aa:bb] =
          paste("z = ", 
                round(IPDA.full.random$coefficients[, "z value"][output_table$label[aa:bb]], 4)
          )
        output_table$es.df.random[aa:bb] = NA
        
        output_table$es.p.random[aa:bb] =                    
          round(IPDA.full.random$coefficients[, "Pr(>|z|)"][output_table$label[aa:bb]], 4)
        
        
        
        output_table$label[bb+1] = "Test of Simple Effects"
        
        #### simple effects 
        for (i in 1: nlevels(data.analysis$Moderator))
        {
          entry = bb+ 1 + i
          output_table$label[entry] = levels(data.analysis$Moderator)[i]
          
          contrast_text = 
            paste("round(IPDA.full.fixed.contrast.", i, "$coefficients[2], 4)", sep = "")
          output_table$es.fixed[entry] = 
            eval(parse(text = contrast_text))
          
          contrast_text = 
            paste('round(IPDA.full.fixed.contrast.', i, '$coefficients[, "z value"][2], 4)', sep = "")
          output_table$es.stat.fixed[entry] = 
            paste("z = ", eval(parse(text = contrast_text)), sep = "")
          
          contrast_text = 
            paste('round(IPDA.full.fixed.contrast.', i, '$df[2], 4)', sep = "")
          output_table$es.df.fixed[entry] = NA
          
          contrast_text = 
            paste('round(IPDA.full.fixed.contrast.', i, '$coefficients[, "Pr(>|z|)"][2], 4)', sep = "")
          output_table$es.p.fixed[entry] = 
            eval(parse(text = contrast_text))
          
          ####
          contrast_text = 
            paste("round(IPDA.full.random.contrast.", i, "$coefficients[2], 4)", sep = "")
          output_table$es.random[entry] = 
            eval(parse(text = contrast_text))
          
          contrast_text = 
            paste('round(IPDA.full.random.contrast.', i, '$coefficients[, "z value"][2], 4)', sep = "")
          output_table$es.stat.random[entry] = 
            paste("z = ", eval(parse(text = contrast_text)), sep = "")
          
          contrast_text = 
            paste('round(IPDA.full.random.contrast.', i, '$coefficients[, "df"][2], 4)', sep = "")
          output_table$es.df.random[entry] = NA
          
          contrast_text = 
            paste('round(IPDA.full.random.contrast.', i, '$coefficients[, "Pr(>|z|)"][2], 4)', sep = "")
          output_table$es.p.random[entry] = 
            eval(parse(text = contrast_text))
        }
        
      }
    }
    
    if (Moderator.nature == "continuous" & DV.nature == "continuous")
    {
      #source("output moderation continuous DV.con.R")
      {##### Moderation output table when Moderator is nominal
        
        #### beta of moderators####
        predictor_list = names(IPDA.full.fixed$coefficients[,1])
        moderation_list = paste(predictor_list[2], ":", predictor_list[3], sep = "")
        
        output_table$label[4] = moderation_list
        
        output_table$es.fixed[4] =   
          IPDA.full.fixed$coefficients[, 1][moderation_list]
        
        output_table$es.stat.fixed[4] =
          paste("t = ", 
                round(IPDA.full.fixed$coefficients[, "t value"][moderation_list], 4)
          )
        output_table$es.df.fixed[4] = IPDA.full.fixed$df[2]
        
        output_table$es.p.fixed[4] =                    
          round(IPDA.full.fixed$coefficients[, "Pr(>|t|)"][moderation_list], 4)
        
        output_table$es.random[4] =   
          IPDA.full.random$coefficients[, 1][moderation_list]
        
        output_table$es.stat.random[4] =
          paste("t = ", 
                round(IPDA.full.random$coefficients[, "t value"][moderation_list], 4)
          )
        output_table$es.df.random[4] = 
          round(IPDA.full.random$coefficients[, "df"][moderation_list], 4)
        
        output_table$es.p.random[4] =                    
          round(IPDA.full.random$coefficients[, "Pr(>|t|)"][moderation_list], 4)
        
        
        
        #### simple effects 
        output_table$label[5:7] = 
          c("Test of Simple Effects",
            "-1 SD",
            "+1 SD"
          )
        
        output_table$es.fixed[6] =   
          IPDA.full.minus1SD.fixed$coefficients[, 1][predictor_list[2]]
        
        output_table$es.stat.fixed[6] =
          paste("t = ", 
                round(IPDA.full.minus1SD.fixed$coefficients[, "t value"][predictor_list[2]], 4)
          )
        output_table$es.df.fixed[6] = IPDA.full.minus1SD.fixed$df[2]
        
        output_table$es.p.fixed[6] =                    
          round(IPDA.full.minus1SD.fixed$coefficients[, "Pr(>|t|)"][predictor_list[2]], 4)
        
        output_table$es.random[6] =   
          IPDA.full.minus1SD.random$coefficients[, 1][predictor_list[2]]
        
        output_table$es.stat.random[6] =
          paste("t = ", 
                round(IPDA.full.minus1SD.random$coefficients[, "t value"][predictor_list[2]], 4)
          )
        output_table$es.df.random[6] = 
          round(IPDA.full.minus1SD.random$coefficients[, "df"][predictor_list[2]], 4)
        
        output_table$es.p.random[6] =                    
          round(IPDA.full.minus1SD.random$coefficients[, "Pr(>|t|)"][predictor_list[2]], 4)
        
        
        output_table$es.fixed[7] =   
          IPDA.full.plus1SD.fixed$coefficients[, 1][predictor_list[2]]
        
        output_table$es.stat.fixed[7] =
          paste("t = ", 
                round(IPDA.full.plus1SD.fixed$coefficients[, "t value"][predictor_list[2]], 4)
          )
        output_table$es.df.fixed[7] = IPDA.full.plus1SD.fixed$df[2]
        
        output_table$es.p.fixed[7] =                    
          round(IPDA.full.plus1SD.fixed$coefficients[, "Pr(>|t|)"][predictor_list[2]], 4)
        
        output_table$es.random[7] =   
          IPDA.full.plus1SD.random$coefficients[, 1][predictor_list[2]]
        
        output_table$es.stat.random[7] =
          paste("t = ", 
                round(IPDA.full.plus1SD.random$coefficients[, "t value"][predictor_list[2]], 4)
          )
        output_table$es.df.random[7] = 
          round(IPDA.full.plus1SD.random$coefficients[, "df"][predictor_list[2]], 4)
        
        output_table$es.p.random[7] =                    
          round(IPDA.full.plus1SD.random$coefficients[, "Pr(>|t|)"][predictor_list[2]], 4)
        
        
        
      }
    }
    
    if (Moderator.nature == "continuous" & DV.nature == "nominal")
    {
      #source("output moderation continuous DV.nom.R")
      {##### Moderation output table when Moderator is nominal
        
        #### beta of moderators####
        predictor_list = names(IPDA.full.fixed$coefficients[,1])
        moderation_list = paste(predictor_list[2], ":", predictor_list[3], sep = "")
        
        output_table$label[4] = moderation_list
        
        output_table$es.fixed[4] =   
          IPDA.full.fixed$coefficients[, 1][moderation_list]
        
        output_table$es.stat.fixed[4] =
          paste("z = ", 
                round(IPDA.full.fixed$coefficients[, "z value"][moderation_list], 4)
          )
        output_table$es.df.fixed[4] = NA
        
        output_table$es.p.fixed[4] =                    
          round(IPDA.full.fixed$coefficients[, "Pr(>|z|)"][moderation_list], 4)
        
        output_table$es.random[4] =   
          IPDA.full.random$coefficients[, 1][moderation_list]
        
        output_table$es.stat.random[4] =
          paste("z = ", 
                round(IPDA.full.random$coefficients[, "z value"][moderation_list], 4)
          )
        output_table$es.df.random[4] = NA
        
        output_table$es.p.random[4] =                    
          round(IPDA.full.random$coefficients[, "Pr(>|z|)"][moderation_list], 4)
        
        
        
        #### simple effects 
        output_table$label[5:7] = 
          c("Test of Simple Effects",
            "-1 SD",
            "+1 SD"
          )
        
        output_table$es.fixed[6] =   
          IPDA.full.minus1SD.fixed$coefficients[, 1][predictor_list[2]]
        
        output_table$es.stat.fixed[6] =
          paste("z = ", 
                round(IPDA.full.minus1SD.fixed$coefficients[, "z value"][predictor_list[2]], 4)
          )
        output_table$es.df.fixed[6] = NA
        
        output_table$es.p.fixed[6] =                    
          round(IPDA.full.minus1SD.fixed$coefficients[, "Pr(>|z|)"][predictor_list[2]], 4)
        
        output_table$es.random[6] =   
          IPDA.full.minus1SD.random$coefficients[, 1][predictor_list[2]]
        
        output_table$es.stat.random[6] =
          paste("z = ", 
                round(IPDA.full.minus1SD.random$coefficients[, "z value"][predictor_list[2]], 4)
          )
        output_table$es.df.random[6] = NA
        
        output_table$es.p.random[6] =                    
          round(IPDA.full.minus1SD.random$coefficients[, "Pr(>|z|)"][predictor_list[2]], 4)
        
        
        output_table$es.fixed[7] =   
          IPDA.full.plus1SD.fixed$coefficients[, 1][predictor_list[2]]
        
        output_table$es.stat.fixed[7] =
          paste("z = ", 
                round(IPDA.full.plus1SD.fixed$coefficients[, "z value"][predictor_list[2]], 4)
          )
        output_table$es.df.fixed[7] = NA
        
        output_table$es.p.fixed[7] =                    
          round(IPDA.full.plus1SD.fixed$coefficients[, "Pr(>|z|)"][predictor_list[2]], 4)
        
        output_table$es.random[7] =   
          IPDA.full.plus1SD.random$coefficients[, 1][predictor_list[2]]
        
        output_table$es.stat.random[7] =
          paste("z = ", 
                round(IPDA.full.plus1SD.random$coefficients[, "z value"][predictor_list[2]], 4)
          )
        output_table$es.df.random[7] = NA
        
        output_table$es.p.random[7] =                    
          round(IPDA.full.plus1SD.random$coefficients[, "Pr(>|z|)"][predictor_list[2]], 4)
        
        
        
      }
      
    }
  }
  
  
  
  #@the fourth row: test of heterogeneity
  for (i in 1:nrow(output_table))
  {
    if (is.na(output_table$label[i]) == T)
    {
      break
    }
  }
  output_table$label[i] = "Test of Heterogeneity"
  
  
  if (DV.nature == "continuous")
  {
    m0.fixed = glm(DV ~ IV + study.cluster, data = data.analysis)
    m1.fixed = glm(DV ~ IV * study.cluster, data = data.analysis)
    
  }
  
  if (DV.nature == "nominal")
  {
    m0.fixed = glm(DV ~ IV + study.cluster,
                   data = data.analysis,
                   family = binomial(link = 'logit'))
    m1.fixed = glm(DV ~ IV * study.cluster,
                   data = data.analysis,
                   family = binomial(link = 'logit'))
  }
  
  output_table$es.fixed[i] =
    paste("pesudo r^2 change = ",
          round(PseudoR2(m1.fixed) - PseudoR2(m0.fixed), 4),
          sep = "")
  
  m0vsm1 =
    anova(m0.fixed, m1.fixed)
  
  output_table$es.stat.fixed[i] =
    paste("chi-squared change = ",
          round(m0vsm1$Deviance[2], 4), sep = "")
  output_table$es.df.fixed[i] =
    round(m0vsm1$Df[2], 4)
  m0vsm1.p =
    pchisq(m0vsm1$Deviance[2],
           m0vsm1$Df[2], lower.tail = F)
  output_table$es.p.fixed[i] = round(m0vsm1.p, 4)
  
  for (i in 1:nrow(output_table))
  {
    if (is.na(output_table$label[i]) == T)
    {
      break
    }
  }
  
  output_table = output_table[-c(1)]
  output_table = output_table[1:i - 1,]
  output_table[is.na(output_table)] = ""
  
  colnames(output_table) = 
    c("",
      "estimated fixed",
      "test fixed",
      "df fixed",
      "p fixed", 
      "estimated random",
      "test random",
      "df random",
      "p random" 
    )
  
  
  a = kable(output_table, 
            col.names = c("",
                          "Estimate",
                          "Test",
                          "df",
                          "p",
                          "Estimate",
                          "Test",
                          "df",
                          "p"), 
                          align = "lrrrrrrrr"
  ) %>% kableExtra::kable_styling() %>% 
    kable_classic(full_width = F, html_font = "Arial", font_size = 14) %>% 
    column_spec(., 1, width = "4cm") %>% column_spec(., 2:9, width = "2cm") %>% 
    add_header_above(c(
      " " = 1,
      "Fixed Effects" = 4,
      "Random Effects" = 4))
      
  return(a) 
}