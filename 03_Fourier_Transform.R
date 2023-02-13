# written by Liangying, 1/31/2023

#install.packages("abind")
#install.packages("ggsignif")
#install.packages("ggpubr")
library(abind)
library(ggsignif)
library(rstatix)

df_DA = cbind(DA, df_trials)
df_5HT = cbind(HT, df_trials)
df_NE = cbind(NE, df_trials)
df_pH = cbind(pH, df_trials)
#df_DA = data.frame(DA, df_trials)

#-------------------FFT----------------------
Fs = 10 
sub_level = 1
col_start = 1
col_end = 57
#sub = "014"

#--------------single subject level---------------------
sub_level = 1
power = 0
t_p_list_allsubs = list()
df_Mag_half_valence_allsubs = list()
sub_list = unique(df_trials$sub_id_num)

for(sub in sub_list)
{
  fft_list = Fourier(df_DA, Fs, sub_level, col_start, col_end, sub, power)
  
  t_p_sig = fft_list$t_p_sig
  l = list(t_p_sig)
  names(l) = str_c("sub", sub)
  t_p_list_allsubs = append(t_p_list_allsubs, l)
  
  Mag_tmp = list(fft_list$df_Mag_half_valence)
  names(Mag_tmp) =  str_c("sub", sub)
  df_Mag_half_valence_allsubs = append(df_Mag_half_valence_allsubs, Mag_tmp)
}

df_Mag_half_valence_allsubs
t_p_list_allsubs

#------------------group level------------------------------
sub_level = 0
power = 0
group_DA = Fourier(df_DA, Fs, sub_level, col_start, col_end, sub, power)
group_5HT = Fourier(df_5HT, Fs, sub_level, col_start, col_end, sub, power)
group_NE = Fourier(df_NE, Fs, sub_level, col_start, col_end, sub, power)
group_pH = Fourier(df_pH, Fs, sub_level, col_start, col_end, sub, power)

power = 1
group_power = Fourier(df_pH, Fs, sub_level, col_start, col_end, sub, power)

#--------Fourier transform function---------
Fourier = function(df, Fs, sub_level, col_start, col_end, sub, power)
{
  if(sub_level == 1)
  {
    d = filter(df, sub_id_num == sub) # attention, the formal parameter cannot be the same with the actual parameter
    df_sub = d[, col_start:col_end]
  } else {
    df_sub = df[, col_start:col_end]
  }  # notice the curly brackets the location of else, otherwise R will think it's finished
  
  df_freq = apply(df_sub, 1, fft)
  df_freq = t(df_freq)
  df_Mag = Mod(df_freq)   # Magnitude

  if(power == 1)
  {
    df_Mag = log10(Mod(df_freq)^2)
  }
  
  N = col_end - col_start + 1
  dF = Fs / N
  #f = dF * seq(0, (N/2), 1)
  #f = dF * seq(0, (N-1), 1)
  f_half = dF * seq(0, floor(N/2), 1)   # for both even and odd numbers
  
  df_Mag_half = df_Mag[, 1: (floor(N/2) + 1)]
  if(sub_level == 1)
  {
    df_Mag_half_valence = data.frame(Mag = df_Mag_half, valence = d$valence)
  } else {
    df_Mag_half_valence = data.frame(Mag = df_Mag_half, valence = df_DA$valence)
  }
  
  t_half = sapply(df_Mag_half_valence[, 1: (floor(N/2) + 1)], function(cl) {t.test(cl ~ valence, df_Mag_half_valence) })
  
  my_list = list("df_Mag_half_valence" = df_Mag_half_valence, "f_half" = f_half, "t_half" = t_half)   #R cannot directly return multiple values like python
  #return(my_list)
  
  df_Mag_half_valence = my_list$df_Mag_half_valence
  f_half = my_list$f_half
  f_half = round(f_half, 2)
  t_half = my_list$t_half
  
  f_half_str = as.character(f_half)
  colnames(t_half) = f_half_str     # make frequency as the colnames
  
  p_half = apply(t_half, 2, function(cl){cl$p.value})
  t_half_tmp = apply(t_half, 2, function(cl){cl$statistic})
  p_half = data.frame(as.list(p_half))
  colnames(p_half) = f_half_str
  
  freedom_half = apply(t_half, 2, function(cl){cl$parameter})
  
  p_half_corrected = p.adjust(p_half, method = "fdr", n = length(p_half))
  p_half_corrected = data.frame(as.list(p_half_corrected))
  colnames(p_half_corrected) = f_half_str
  
  p_half_sig = p_half[which(p_half < 0.05)] # find which frequency's p value < 0.05
  p_half_corrected_sig = p_half_corrected[which(p_half < 0.05)]
  t_half_sig = t_half_tmp[which(p_half < 0.05)]
  
  freedom_half_sig = freedom_half[which(p_half < 0.05)]
  
  t_half_sig = data.frame(as.list(t_half_sig))
  colnames(t_half_sig) = colnames(p_half_sig)
  
  t_p_sig = rbind(t_half_sig, p_half_sig, p_half_corrected_sig, freedom_half_sig)
  
  fft_list = list("df_Mag_half_valence" = df_Mag_half_valence, "t_p_sig" = t_p_sig)
  return(fft_list)
}


# the above two independent t-test is incorrect from statistics view, so I use LME 

#-----------------Linear Mixed Effect model for Repeated measures two sample t-test-----------------------
colnames(group_DA$df_Mag_half_valence)

data = group_DA
data = group_5HT 
data = group_NE 
data = group_pH

sub = df_trials$sub_id_num
f_list_str

result_lme_2T = lme_2T(data, sub, f_list_str)
result_lme_2T= result_lme_2T$lme_2T

cols = which(result_lme_2T['p',] < 0.05)
tmp = result_lme_2T[cols]    # [], dataframe
format(tmp, scientific = F)   # get rid of scientific notation

p.adjust(result_lme_2T['p',], method = "fdr", n = length(result_lme_2T['p',]))


lme_2T = function(data, sub, f_list_str)
{
  Mag = data$df_Mag_half_valence
  colnames(Mag) = c(f_str, "valence")
  Mag$sub = sub
  #fre = "3.68"
  
  i = 1
  for(fre in f_list_str)
  {
    df = data.frame(Mag_f = Mag[[fre]], valence = factor(Mag$valence), sub = factor(Mag$sub))
    
    model = lmer(Mag_f ~ valence + (1 | sub), df)   
    
    result = summary(model)$coefficients     # fixed effect
    
    # based on neg reference
    beta = result[2,1]
    freedom = result[2,3]
    t = result[2,4]
    p = result[2,5]
    
    # based on pos reference, -t
    beta = -1 * beta
    t = -1 * t
    
    
    if(i == 1)
    {
      df_result = data.frame(c(beta, freedom, t, p))
    }
    else
    {
      df_result = cbind(df_result, c(beta, freedom, t, p))
    }
    i = i + 1
  }
  
  names(df_result) = f_list_str
  rownames(df_result) = c("beta", "df","t","p")
  
  lme_list = list("lme_2T" = df_result)
  
  return(lme_list)
}




#---------Average magnitude of each frequency across condition for each subject--------
i = 1
for(sub in sub_list)
{
  sub = str_c("sub", sub)
  Mag_sub = df_Mag_half_valence_allsubs[[sub]]      # [[var]] is a dataframe, [var] is a list
  Mag_sub_mean = Mag_sub %>% group_by(valence) %>% dplyr::summarise(across(everything(), mean)) %>% as.data.frame()
  
  if(i == 1)
  {
    Mag_mean_allsubs = Mag_sub_mean
  }
  else
  {            
    Mag_mean_allsubs = abind(Mag_mean_allsubs, Mag_sub_mean, along = 3)
  }
  i = i + 1 
}

col = dim(Mag_mean_allsubs)[2]
Mag_mean_allsubs_del = Mag_mean_allsubs[, 2:col,]   #delete valence col, otherwise error for t-test
t_mean_allsubs = apply(Mag_mean_allsubs_del, 2, function(d) {t.test(d[1,], d[2,]) })
names(t_mean_allsubs) = f_half_str

t_p_mean_sig = lapply(t_mean_allsubs, 
             function(l)
               {
               if(l$p.value < 0.1)
               {
                 t = l$statistic
                 p = l$p.value
                 c(t,p)   # return multiple values
               }
               })
t_p_mean_sig[sapply(t_p_mean_sig, is.null)] = NULL


#------------Plot Mean Magnitude and Power of each frequency--------------
i = 1
k = 1
for(sub in sub_list)
{
  #sub = str_c("sub", sub)
  #Mag_sub = Mag_mean_allsubs[[sub]]
  Mag_sub = Mag_mean_allsubs[,,k]
  Mag_sub = as.data.frame(Mag_sub)
  
  new_names = c("valence", f_half_str)
  names(Mag_sub) = new_names
  
  Mag_sub_arr = melt(Mag_sub, id.vars = "valence", measure.vars = f_half_str, variable_name = "freq", value.name = "Mag")
  #Mag_sub_arr$sub = str_sub(sub, start = 4, end = -1)
  Mag_sub_arr$sub = sub
  
  if(i == 1)
  {
    Mag_allsubs_arr = Mag_sub_arr
  }
  else
  {
    Mag_allsubs_arr = rbind(Mag_allsubs_arr, Mag_sub_arr)
  }
  i = i + 1
  k = k + 1
}

colnames(Mag_allsubs_arr)[which(names(Mag_allsubs_arr) == "value")] = "Mag"
df_Mag_mean_allsubs = Mag_allsubs_arr %>% group_by(sub, freq, valence) %>% dplyr::summarise(Mag = mean(Mag)) %>% as.data.frame()

data = df_Mag_mean_allsubs
data2 = summarySE(df_Mag_mean_allsubs, measurevar = "Mag", groupvars = c("freq","valence"))

#--------Power = log10Mag---------------------
Mag_allsubs_arr$Power = log10(Mag_allsubs_arr$Mag^2)
#Mag_allsubs_arr$Power = Mag_allsubs_arr$Mag^2
df_Power_mean_allsubs = Mag_allsubs_arr %>% group_by(sub, freq, valence) %>% dplyr::summarise(Power = mean(Power)) %>% as.data.frame()

tmp = Mag_allsubs_arr %>% filter(freq == 3.68)
t.test(Power~valence, tmp) # not mean, based on all trials

tmp = df_Power_mean_allsubs %>% filter(freq == 3.68)
t.test(Power~valence, tmp)

data = df_Power_mean_allsubs %>% filter(freq != 0)
data2 = summarySE(data, measurevar = "Power", groupvars = c("freq","valence"))

t.test(df_Mag_half_valence[,21] ~ valence, df_Mag_half_valence)

#---------plot Magnitude of each frequency, no mean------------
df_Mag_half_valence = group$df_Mag_half_valence
df_Mag_half_valence = group_5HT$df_Mag_half_valence
df_Mag_half_valence = group_NE$df_Mag_half_valence
df_Mag_half_valence = group_pH$df_Mag_half_valence

new_names = c(f_half_str, "valence")
names(df_Mag_half_valence) = new_names

df_Mag_half_valence_arr = melt(df_Mag_half_valence, id.vars = "valence", measure.vars = f_half_str, variable_name = "freq", value.name = "Mag")
colnames(df_Mag_half_valence_arr)[which(names(df_Mag_half_valence_arr) == "value")] = "Mag"

data = df_Mag_half_valence_arr
data2 = summarySE(df_Mag_half_valence_arr, measurevar = "Mag", groupvars = c("freq","valence"))

max(data$Mag)

Magplot(data,data2) 

Magplot = function(data, data2)
{
plot = ggplot(data2, mapping = aes(x = freq, y = Mag, fill = factor(valence)))+
    geom_bar(position = position_dodge(), stat='identity', color='black', size=.75, width = 0.75)+
    geom_point(data, mapping = aes(x = freq, y = Mag, color = factor(valence)), position = position_jitterdodge(), size = 1)+ # make dots on the top of errorbar
    geom_errorbar(aes(ymin = Mag-se, ymax = Mag+se), size=.9,width=.2,position=position_dodge(.75))+ 
    geom_signif(stat = "identity",
                data = data.frame(x = c(3.6, 4.6, 15.6, 21.6, 25.5, 26.6, 27.6), xend = c(4.3, 5.3, 16.4, 22.4, 26.3, 27.4,28.5),
                                  y = c(14.6, 14.3, 3.1, 1.8, 1.7, 1.75, 1.9), annotation = c(" * ", "***","*", "**", " ** ", "* ", " *")),
              # data = data.frame(x = c(4.6, 5.7, 28.7), xend = c(5.3, 6.5, 29.3),
              #                  y = c(14.2, 16.4, 2), annotation = c("*"," *", "* ")),
              #data = data.frame(x = c(4.6, 21.6, 25.5), xend = c(5.3, 22.4, 26.3),
              #                  y = c(15, 2.2, 2), annotation = c("*","p = 0.056", "* ")),
              #data = data.frame(x = c(18.6), xend = c(19.4),
              #                  y = c(4.2), annotation = c("p = 0.052")),
                aes(x = x, xend = xend, y = y, yend = y, annotation = annotation,textsize = 6))+
    scale_color_manual(values=c("#7FAEFF","#808A87"),name = "Valence",labels = c("Neg", "Pos"))+
    scale_fill_manual(values=c("#4A82FF","#808069"), name = "Valence",labels = c("Neg", "Pos"))+
   # scale_fill_manual(values=c("#4A82FF","#808069"),guide = FALSE)+
    labs(title = "pH")+
    xlab('Frequency')+
    ylab("Magnitude")+
    #scale_x_discrete(expand = c(0,0))+
    #scale_y_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0), limits = c(0, ceiling(max(data$Mag))))+ 
    #scale_y_continuous(expand = c(0,0), limits = c(-0.5, 2.5))+ 
    theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.line = element_line(color='black',size = 1),
          axis.text.x = element_text(size = 14, face = "bold", angle = -30, vjust = 1, hjust = 0),
          axis.text.y = element_text(size = 16,face = "bold"),
          axis.ticks.y = element_line(size = 1), 
          axis.ticks.x = element_blank(), #删去刻度线
          axis.title.x = element_text(size = 16, face = "bold",margin = margin(t = 10, r = 0, b = 0, l = 0)),   # change the distance between the axis title and the numbers. 
          axis.title.y = element_text(size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),   # change the distance between the axis title and the numbers. 
          axis.ticks.length.y = unit(0.15, 'cm'),
          legend.text = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 16, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
}

plot = ggplot(data2, mapping = aes(x = freq, y = Power, fill = factor(valence)))+
  geom_bar(position = position_dodge(), stat='identity', color='black', size=.75, width = 0.75)+
  geom_point(data, mapping = aes(x = freq, y = Power, color = factor(valence)), position = position_jitterdodge(), size = 1)+ # make dots on the top of errorbar
  geom_errorbar(aes(ymin = Power-se, ymax = Power+se), size=.9,width=.2,position=position_dodge(.75))+ 
  geom_signif(stat = "identity",
              data = data.frame(x = 20.7, xend = 21.5,
                                y = 0.4, annotation = "p = 0.076"),
              aes(x = x, xend = xend, y = y, yend = y, annotation = annotation,textsize = 5))+
  scale_color_manual(values=c("#7FAEFF","#808A87"),name = "Valence",labels = c("Neg", "Pos"))+
  scale_fill_manual(values=c("#4A82FF","#808069"), name = "Valence",labels = c("Neg", "Pos"))+
  # scale_fill_manual(values=c("#4A82FF","#808069"),guide = FALSE)+
  labs(title = "DA")+
  xlab('Frequency')+
  ylab("Power")+
  #scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0), limits = c(-0.7, ceiling(max(data$Power))))+ 
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(color='black',size = 1),
        axis.text.x = element_text(size = 14, face = "bold", angle = -30, vjust = 1, hjust = 0),
        axis.text.y = element_text(size = 16,face = "bold"),
        axis.ticks.y = element_line(size = 1), 
        axis.ticks.x = element_blank(), #删去刻度线
        axis.title.x = element_text(size = 16, face = "bold",margin = margin(t = 10, r = 0, b = 0, l = 0)),   # change the distance between the axis title and the numbers. 
        axis.title.y = element_text(size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),   # change the distance between the axis title and the numbers. 
        axis.ticks.length.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))


root_path = "/Users/liangyingliu/Desktop/Results/Fourier/5-HT"
path = file.path(root_path, "Power_allsubs")
ggsave(plot = plot, filename = str_c("Power_allsubs", ".PNG"), path = path)

root_path = "/Users/liangyingliu/Desktop/Results/Fourier/pH"
path = file.path(root_path, "Mag_allsubs_mean")
ggsave(plot = plot, filename = str_c("Mag_allsubs_mean", ".PNG"), path = path)



#----------plot with single trials------------------
plot_valence_all <- ggplot(df_d %>% filter(sub_id_num == sub, between(trials, 16,36)), mapping = aes(x = Timepoint, y = DA, color = factor(valence)))+
  geom_line()+
  labs(title = sub_name)+
  #facet_wrap(~ valence + session_id + blocks + trials, scales = "fixed", 
  #           labeller = labeller(valence = valence_labs, session_id = session_labs, blocks = blocks_labs))+
  facet_grid(trials ~ valence + session_id + blocks,
             labeller = labeller(valence = valence_labs, session_id = session_labs, blocks = blocks_labs))+
  geom_vline(aes(xintercept = 26),
             color = "black",
             alpha = 0.7)+
  theme_light()+
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 7)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.position = "none")
plot_valence_all




#-------------------------------A function for mean, count, sd, se, and ci-------------------------------------------------
##   Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
