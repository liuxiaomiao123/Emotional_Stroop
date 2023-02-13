# written by Liangying, 2/6/2023

#install.packages("scatterplot3d")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("emmeans")
library(scatterplot3d)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)


#------------------group level FFT, see 03_Fourier_Transform.R code------------------------------
sub_level = 0
power = 0
group_DA = Fourier(df_DA, Fs, sub_level, col_start, col_end, sub, power)
group_5HT = Fourier(df_5HT, Fs, sub_level, col_start, col_end, sub, power)
group_NE = Fourier(df_NE, Fs, sub_level, col_start, col_end, sub, power)
group_pH = Fourier(df_pH, Fs, sub_level, col_start, col_end, sub, power)

Mag_DA = group_DA$df_Mag_half_valence
Mag_5HT = group_5HT$df_Mag_half_valence
Mag_NE = group_NE$df_Mag_half_valence
Mag_pH = group_pH$df_Mag_half_valence

new_names = c(f_half_str, "valence")
names(Mag_DA) = new_names
names(Mag_5HT) = new_names
names(Mag_NE) = new_names
names(Mag_pH) = new_names

Mag_DA_lme = cbind(Mag_DA, sub_id_num = df_trials$sub_id_num)
Mag_5HT_lme = cbind(Mag_5HT, sub_id_num = df_trials$sub_id_num)
Mag_NE_lme = cbind(Mag_NE, sub_id_num = df_trials$sub_id_num)
Mag_pH_lme = cbind(Mag_pH, sub_id_num = df_trials$sub_id_num)


#--------------------------test DA, 5HT, NE, pH colinearity--------------------------
DA_5HT_lme = Colinearity_lme(Mag_DA_lme, Mag_5HT_lme)
HT_DA_lme = Colinearity_lme(Mag_5HT_lme, Mag_DA_lme)

DA_NE_lme = Colinearity_lme(Mag_DA_lme, Mag_NE_lme)
NE_DA_lme = Colinearity_lme(Mag_NE_lme, Mag_DA_lme)

DA_pH_lme = Colinearity_lme(Mag_DA_lme, Mag_pH_lme)
pH_DA_lme = Colinearity_lme(Mag_pH_lme, Mag_DA_lme)

HT_NE_lme = Colinearity_lme(Mag_5HT_lme, Mag_NE_lme)
NE_5HT_lme = Colinearity_lme(Mag_NE_lme, Mag_5HT_lme)

HT_pH_lme = Colinearity_lme(Mag_5HT_lme, Mag_pH_lme)
pH_5HT_lme = Colinearity_lme(Mag_pH_lme, Mag_5HT_lme)

NE_pH_lme = Colinearity_lme(Mag_NE_lme, Mag_pH_lme)
pH_NE_lme = Colinearity_lme(Mag_pH_lme, Mag_NE_lme)

pair_lme_list = list(DA_5HT_lme, HT_DA_lme, DA_NE_lme, NE_DA_lme, DA_pH_lme, pH_DA_lme, HT_NE_lme, NE_5HT_lme, HT_pH_lme, pH_5HT_lme, NE_pH_lme, pH_NE_lme)
pair_name = c("DA_5HT_lme", "HT_DA_lme", "DA_NE_lme", "NE_DA_lme", "DA_pH_lme", "pH_DA_lme", "HT_NE_lme", "NE_5HT_lme", "HT_pH_lme", "pH_5HT_lme", "NE_pH_lme", "pH_NE_lme")
names(pair_lme_list) = pair_name


idx = "interaction"
i = 1

for(pair in pair_name)
{
  ls = Heatmap_row(pair_lme_list[[pair]], idx)
  t = ls[[1]]
  p = ls[[2]]
  
  if(i == 1)
  {
    t_heatmap = t
    p_heatmap = p
    
  }else{
    t_heatmap = rbind(t_heatmap, t)
    p_heatmap = rbind(p_heatmap, p)
  }
  
  i = i + 1
}
rownames(t_heatmap) = pair_name
rownames(p_heatmap) = pair_name

if(idx == "neg")
{
  t_heatmap_neg = t_heatmap
  p_heatmap_neg = p_heatmap
}else if(idx == "pos"){
  t_heatmap_pos = t_heatmap
  p_heatmap_pos = p_heatmap
}else{
  t_heatmap_interaction = t_heatmap
  p_heatmap_interaction = p_heatmap
}


#-------------------------------------------------------------------
Heatmap_row = function(lme_result, idx)
{
  if(idx == "neg")
  {
     t_neg = lme_result$lme_neg['t',]
     p_neg = lme_result$lme_neg['p',]
     neg_list = list("t_neg" = t_neg, "p_neg" = p_neg)
     return(neg_list)
  } else if(idx == "pos") {
    t_pos = lme_result$lme_pos['t',]
    p_pos = lme_result$lme_pos['p',]
    pos_list = list("t_pos" = t_pos, "p_pos" = p_pos)
    return(pos_list)
  }else{
    t_interaction = lme_result$lme_interaction['t',]
    p_interaction = lme_result$lme_interaction['p',]
    interaction_list = list("t_interaction" = t_interaction, "p_interaction" = p_interaction)
    return(interaction_list)
  }
}

#----------------pair of neurotransmitter colinearity test, LME with interaction------------------------
f_list_str = as.character(f_list)

Mag1 = Mag_NE_lme
Mag2 = Mag_5HT_lme
fre = "3.16"

Colinearity_lme = function(Mag1, Mag2)   
{
  i = 1
  for(fre in f_list_str)
  {
    df = data.frame(Mag1_f = Mag1[[fre]], Mag2_f = Mag2[[fre]], valence = factor(Mag1$valence), sub = factor(Mag1$sub_id_num))

    #model = lmer(Mag1_f ~ Mag2_f + valence-1 + Mag2_f:valence + (1 |sub), df)   -1: from reference level coding to level means coding, but level means coding can still not get the regression t and p value for each group
    #coef(model)  #check individual slope and intercept
    
    # Relevel can get each group's regression coefficient, df, t and p values 
    model_negRef = lmer(Mag1_f ~ Mag2_f + valence + Mag2_f:valence + (1 |sub), df)    # neg0 is the reference level
    result_negRef = summary(model_negRef)$coefficients     # fixed effect
    
    df$valence = relevel(df$valence, ref = 2)     #neg0, level1; pos1, level2
    model_posRef = lmer(Mag1_f ~ Mag2_f + valence + Mag2_f:valence + (1 |sub), df)    #pos1 is the reference level
    result_posRef = summary(model_posRef)$coefficients     # fixed effect
    
    beta_neg = result_negRef[2,1]
    freedom_neg = result_negRef[2,3]
    t_neg = result_negRef[2,4]
    p_neg = result_negRef[2,5]
    
    beta_pos = result_posRef[2,1]
    freedom_pos = result_posRef[2,3]
    t_pos = result_posRef[2,4]
    p_pos = result_posRef[2,5]
    
    # based on pos reference
    beta_interaction = result_posRef[4,1]
    freedom_interaction = result_posRef[4,3]
    t_interaction = result_posRef[4,4]
    p_interaction = result_posRef[4,5]
    
    if(i == 1)
    {
      df_neg = data.frame(c(beta_neg, freedom_neg, t_neg, p_neg))
      df_pos = data.frame(c(beta_pos, freedom_pos, t_pos, p_pos))
      df_interaction = data.frame(c(beta_interaction, freedom_interaction, t_interaction, p_interaction))
    }
    else
    {
      df_neg = cbind(df_neg, c(beta_neg, freedom_neg, t_neg, p_neg))
      df_pos = cbind(df_pos, c(beta_pos, freedom_pos, t_pos, p_pos))
      df_interaction = cbind(df_interaction, c(beta_interaction, freedom_interaction, t_interaction, p_interaction))
    }
    i = i + 1
  }
  
  names(df_neg) = f_list_str
  names(df_pos) = f_list_str
  names(df_interaction) = f_list_str
  
  rownames(df_neg) = c("beta", "df","t","p")
  rownames(df_pos) = c("beta", "df","t","p")
  rownames(df_interaction) = c("beta", "df","t","p")
  
  lme_list = list("lme_neg" = df_neg, "lme_pos" = df_pos, "lme_interaction" = df_interaction)
  return(lme_list)
}


#--------------------------------------LME plot----------------------------------------------------------
model = model_negRef

model_coef = summary(model)
str(model_coef)
coef_fixed = model_coef$coefficients   # fixed effects
Intercept_neg = coef_fixed[1,1]
Intercept_pos = coef_fixed[3,1] + Intercept_neg
slope_neg = coef_fixed[2,1]
slope_pos = coef_fixed[4,1] + slope_neg

Intercept = c(Intercept_neg, Intercept_pos)
slope = c(slope_neg, slope_pos)
g = c(0,1)
df_coef_fixed <- data.frame(Intercept, slope, valence = factor(g))     # based on the same valence to facet


#-----------Facet plot is better----------------------
valence_labs = c("Neg", "Pos")
names(valence_labs) = c("0", "1")

model_plot <- ggplot(data = df, mapping = aes(Mag1_f, Mag2_f, color = sub))+
  geom_point()+
  facet_wrap(~ valence, scale = "fixed", labeller = labeller(valence = valence_labs))+
  geom_abline(data = df_coef_fixed, aes(intercept = Intercept, slope = slope) , size = 1, color = "darkcyan")+   
  labs(x = "NE", y = "5-HT", title = "3.16Hz")+
  #theme_light()+
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 11,face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 7)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 15)))  # margin: t, r, b, l

model_plot

root_path = "/Users/liangyingliu/Desktop/Results/Fourier"
path = file.path(root_path, "LME")
ggsave(plot = model_plot, filename = str_c("NE_5HT_3.16hz", ".PNG"), path = path)



#------------------------------Heatmap t value plot------------------------------------------------
#-----------------make t value whose p value > 0.05 as NA------------------------------------
p_FDR = function(p)
{
  tmp = as.vector(t(p))    # row1, row2, row3...
  p_adjust  <- p.adjust(tmp, method = "fdr", n = length(tmp))
  p_adjust  = matrix(p_adjust, n_pairs,n_frequency, byrow = TRUE)  
  p_adjust = as.data.frame(p_adjust)
  return(p_adjust)
}

p_heatmap_neg_adjust = p_FDR(p_heatmap_neg)
p_heatmap_pos_adjust = p_FDR(p_heatmap_pos)
p_heatmap_interaction_adjust = p_FDR(p_heatmap_interaction)

indx  = which(p_heatmap_neg_adjust < 0.05, arr.ind = TRUE)
head(indx)

indx  = which(p_heatmap_neg < 0.05, arr.ind = TRUE)
head(indx)

t_heatmap_neg_threshold = t_heatmap_neg
t_heatmap_neg_threshold[p_heatmap_neg_adjust > 0.05] = NA   # no need to use which
df_heatmap = t_heatmap_neg_threshold

t_heatmap_pos_threshold = t_heatmap_pos
t_heatmap_pos_threshold[p_heatmap_pos_adjust > 0.05] = NA   # no need to use which
df_heatmap = t_heatmap_pos_threshold

t_heatmap_interaction_threshold = t_heatmap_interaction
t_heatmap_interaction_threshold[p_heatmap_interaction > 0.05] = NA   # no need to use which
df_heatmap = t_heatmap_interaction_threshold


which(p_heatmap_interaction_adjust == min(p_heatmap_interaction_adjust), arr.ind = TRUE)

p_vec = as.vector(t(p_heatmap_interaction))
top5 = p_vec[order(p_vec, decreasing = FALSE)[1:5]]
loc_top1 = which(p_heatmap_interaction == top5[1], arr.ind = TRUE)
loc_top2 = which(p_heatmap_interaction == top5[2], arr.ind = TRUE)
loc_top3 = which(p_heatmap_interaction == top5[3], arr.ind = TRUE)
loc_top4 = which(p_heatmap_interaction == top5[4], arr.ind = TRUE)
loc_top5 = which(p_heatmap_interaction == top5[5], arr.ind = TRUE)


#--------------------------------------------------------------------------------------------
cols = colnames(df_heatmap)
df_heatmap_arr = melt(df_heatmap, measure.vars = cols, variable_name = "freq", value.name = "t")
colnames(df_heatmap_arr)[which(names(df_heatmap_arr) == "value")] = "t"

n_pairs = dim(df_heatmap)[1]
n_frequency = dim(df_heatmap)[2]

x_fre = rep(1:n_frequency, each = n_pairs)   # each
df_heatmap_arr$x_fre = x_fre

y_pair = rep(1:n_pairs, time = n_frequency)   # time
df_heatmap_arr$y_pair = y_pair

options(repr.plot.width = 1, repr.plot.height = 1)

#pair_name2 = str_sub(pair_name, start = 1, end = -5)
pair_name2 = c("DA_5HT", "5HT_DA", "DA_NE", "NE_DA", "DA_pH", "pH_DA", "5HT_NE", "NE_5HT", "5HT_pH", "pH_5HT", "NE_pH", "pH_NE")

min(df_heatmap_arr$t)


t_heatmap <- ggplot(data = df_heatmap_arr, mapping = aes(factor(x_fre), factor(y_pair), fill = t))+
  #geom_tile()+
  #scale_fill_gradient(low="red", high="white") 
  geom_raster()+
  geom_tile(data = df_heatmap_arr[212,], aes(factor(x_fre), factor(y_pair)), color = "black", fill = "transparent", size = 1)+     # add black border
  #scale_fill_distiller(palette = "Spectral", direction = 1, breaks = c(0.05, 0.50, 0.99), labels = c("0.05", "0.50", "1.00"))+
  scale_fill_distiller(palette = "Spectral", na.value = 'white')+      # make t value whose p value > 0.05 as white 
  scale_x_discrete(expand = c(0,0), breaks = 1:n_frequency, labels = f_list_str)+    # discrete must be factor, not numeric
  scale_y_discrete(expand = c(0,0), breaks = 1:n_pairs, labels = pair_name2)+
  labs(x = "Frequency", y = "", fill = "t value", title = "Neg - Pos")+
  theme_light()+
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = -30, vjust = 1, hjust = 0),
        axis.text.y = element_text(size = 10, face = "bold"),
        # axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 7)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 20)),
        aspect.ratio = 1 / 2.5)
t_heatmap


root_path = "/Users/liangyingliu/Desktop/Results/Fourier"
path = file.path(root_path, "Colinearity")
ggsave(plot = t_heatmap, filename = str_c("Neg-Pos_uncor_blackborder", ".PNG"), path = path)



#----------------From statistics view, this is wrong, but just for comparison with LME------------------------
cor_DA_5HT = Colinearity(Mag_DA, Mag_5HT)
cor_DA_NE = Colinearity(Mag_DA, Mag_NE)
cor_DA_pH = Colinearity(Mag_DA, Mag_pH)

cor_5HT_NE = Colinearity(Mag_5HT, Mag_NE)
cor_5HT_pH = Colinearity(Mag_5HT, Mag_pH)

cor_NE_pH = Colinearity(Mag_NE, Mag_pH)

#----------correlation across subs-----------
f_list_str = as.character(f_list)

Colinearity = function(Mag1, Mag2)
{
  Mag1_neg = Mag1 %>% filter(valence == 0)
  Mag1_pos = Mag1 %>% filter(valence == 1)
  
  Mag2_neg = Mag2 %>% filter(valence == 0)
  Mag2_pos = Mag2 %>% filter(valence == 1)
  
  i = 1
  for(fre in f_list_str)
  {
    cor_neg = cor.test(Mag1_neg[[fre]], Mag2_neg[[fre]])     # not $
    cor_pos = cor.test(Mag1_pos[[fre]], Mag2_pos[[fre]])
    
    r_neg = cor_neg$statistic
    freedom_neg = cor_neg$parameter
    p_neg = cor_neg$p.value
    
    r_pos = cor_pos$statistic
    freedom_pos = cor_pos$parameter
    p_pos = cor_pos$p.value
    
    if(i == 1)
    {
      df_neg = data.frame(c(r_neg, freedom_neg, p_neg))
      df_pos = data.frame(c(r_pos, freedom_pos, p_pos))
    }
    else
    {
      df_neg = cbind(df_neg, c(r_neg, freedom_neg, p_neg))
      df_pos = cbind(df_pos, c(r_pos, freedom_pos, p_pos))
    }
    i = i +1
  }
  
  names(df_neg) = f_list_str
  names(df_pos) = f_list_str
  
  cor_list = list("cor_neg" = df_neg, "cor_pos" = df_pos)
  return(cor_list)
}




#-----------------------------------------------------Other tries-----------------------------------------------------------
group = group_DA
group = group_5HT
group = group_NE
group = group_pH

Mag_all = group$df_Mag_half_valence
colnames(Mag_all)

new_names = c(f_half_str, "valence")
names(Mag_all) = new_names

Mag_DA_all_arr = melt(Mag_all, id.vars = "valence", measure.vars = f_half_str, variable_name = "freq", value.name = "Mag")
colnames(Mag_DA_all_arr)[which(names(Mag_DA_all_arr) == "value")] = "DA"

Mag_5HT_all_arr = melt(Mag_all, id.vars = "valence", measure.vars = f_half_str, variable_name = "freq", value.name = "Mag")
colnames(Mag_5HT_all_arr)[which(names(Mag_5HT_all_arr) == "value")] = "5HT"

Mag_NE_all_arr = melt(Mag_all, id.vars = "valence", measure.vars = f_half_str, variable_name = "freq", value.name = "Mag")
colnames(Mag_NE_all_arr)[which(names(Mag_NE_all_arr) == "value")] = "NE"

Mag_pH_all_arr = melt(Mag_all, id.vars = "valence", measure.vars = f_half_str, variable_name = "freq", value.name = "Mag")
colnames(Mag_pH_all_arr)[which(names(Mag_pH_all_arr) == "value")] = "pH"


#---------------------------Combine dimensions------------------------------
Mag_DA_5HT_NE_pH = cbind(Mag_DA_5HT_NE, pH = Mag_pH_all_arr[,3])
f_list = as.numeric(f_half)

for(f in f_list)
{
  d_tmp = Mag_DA_5HT_NE_pH %>% filter(freq == 3.16) %>% arrange(valence)  # arrange valence to match colors
  colors = c("#E69F00", "#56B4E9")
  colors = colors[as.factor(d_tmp$valence)]
  
  source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
  
  f_str = as.character(f)
  scatterplot3d(d_tmp[, 4:5], xlab = "DA", ylab = "5HT", zlab = "NE", main = f, color = colors, pch = 16, grid = FALSE, box = FALSE)
  addgrids3d(d_tmp[, 3:5], grid = c("xy", "xz", "yz"))
}

ggplot(d_tmp, mapping = aes(x = HT, y = NE, color = factor(valence)))+
  geom_point(size = 1)+ 
  scale_color_manual(values=c("#7FAEFF","#808A87"),name = "Valence",labels = c("Neg", "Pos"))+
  labs(title = "0.7Hz")+
  xlab('5-HT')+
  ylab("NE")+
  #scale_x_discrete(expand = c(0,0))+
  #scale_y_continuous(expand = c(0,0))+
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



for(f in f_list)
{
  d_tmp = Mag_DA_5HT_NE_pH %>% filter(freq == f) %>% arrange(valence)  # arrange valence to match colors
  
  pca = prcomp(d_tmp[, 3:6], center = TRUE, scale. = TRUE)
  #summary(pca)
  
  ggbiplot(pca, ellipse = TRUE,groups = factor(d_tmp$valence), var.axes = FALSE)+
    labs(color = "Cluster")+
    #scale_color_manual(values = c("#868686FF", "#CD534CFF"))+
    theme_light()+
    theme(axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 7)),
          axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10, face = "bold"))
}









