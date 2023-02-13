# written by Liangying, 1/26/2023

library(dtw)
library(dtwclust)
library(tidyverse)
library(reshape)
library(BBmisc)
library(plyr)
library(TSdist)
library(readxl)

#--------------------------------------Data Read + Data Arrangement-------------------------------------------------
path = "/Volumes/proj/montague_trainee_1/L_Liu"
valence_path = file.path(path, "valence.csv")
DA_path = file.path(path, "DA.csv")
NE_path = file.path(path, "NE.csv")
HT_path = file.path(path, "5HT.csv")
pH_path = file.path(path, "pH.csv")
session_id_path = file.path(path, "session_id.csv")
id_mapping_path = file.path(path, "Emostroop_map.xlsx")

valence = read.csv(valence_path, header = FALSE)
DA = read.csv(DA_path, header = FALSE)
NE = read.csv(NE_path, header = FALSE)
HT = read.csv(HT_path, header = FALSE)
pH = read.csv(pH_path, header = FALSE)
session_id = read.csv(session_id_path, header = FALSE)

dDA = t(apply(DA, MARGIN = 1, FUN = Derivative)); dDA = as.data.frame(dDA)
dNE = t(apply(NE, MARGIN = 1, FUN = Derivative)); dNE = as.data.frame(dNE)
dHT = t(apply(HT, MARGIN = 1, FUN = Derivative)); dHT = as.data.frame(dHT)
dpH = t(apply(pH, MARGIN = 1, FUN = Derivative)); dpH = as.data.frame(dpH)

names = c("session_id", "unknown", "unknown", "unknown", "sub_id", "direction", "region", "unknown", "sex", "disease", "unknown")
id_mapping = read_excel(id_mapping_path, col_names = names)
sub_session_id = data.frame(sub_id = id_mapping$sub_id, session_id = id_mapping$session_id)
sub_session_id = na.omit(sub_session_id)
sub_id_num = str_sub(sub_session_id$sub_id, start = 8)
sub_session_id = data.frame(sub_session_id, sub_id_num)   # mapping between sub_id and session_id

n_timepoints = 57
colnames(DA) = c(1:n_timepoints); colnames(dDA) = c(1:n_timepoints)
colnames(NE) = c(1:n_timepoints); colnames(dNE) = c(1:n_timepoints)
colnames(HT) = c(1:n_timepoints); colnames(dHT) = c(1:n_timepoints)
colnames(pH) = c(1:n_timepoints); colnames(dpH) = c(1:n_timepoints)
colnames(valence) = "valence"
colnames(session_id) = "session_id"

#------------------long format------------------
DA_arr = pivot_longer(DA, cols = everything(), names_to = "Timepoint", values_to = "DA")
NE_arr = pivot_longer(NE, cols = everything(), names_to = "Timepoint", values_to = "NE")
HT_arr = pivot_longer(HT, cols = everything(), names_to = "Timepoint", values_to = "HT")
pH_arr = pivot_longer(pH, cols = everything(), names_to = "Timepoint", values_to = "pH")

dDA_arr = pivot_longer(dDA, cols = everything(), names_to = "Timepoint", values_to = "dDA")
dNE_arr = pivot_longer(dNE, cols = everything(), names_to = "Timepoint", values_to = "dNE")
dHT_arr = pivot_longer(dHT, cols = everything(), names_to = "Timepoint", values_to = "dHT")
dpH_arr = pivot_longer(dpH, cols = everything(), names_to = "Timepoint", values_to = "dpH")

valence_arr = rep(valence[['valence']], each = n_timepoints)      # Learn this trick, 1111110000011111100000
session_id_arr = rep(session_id[['session_id']], each = n_timepoints)

df = data.frame(Timepoint = as.numeric(DA_arr$Timepoint),   # 注意这一步非常的关键，画图的x和y一定要是数值型，如果是字符型，1后面跟的不是2，而是10，就会错乱
                DA = DA_arr$DA, 
                NE = NE_arr$NE, 
                HT = HT_arr$HT, 
                pH = pH_arr$pH, 
                valence = valence_arr, 
                session_id = session_id_arr)
df = left_join(df, sub_session_id, by = "session_id")    # combine sub_id info, SQL style

df_d = df
df_d = data.frame(df_d, dDA = dDA_arr$dDA, dNE = dNE_arr$dNE, dHT = dHT_arr$dHT, dpH = dpH_arr$dpH)

#------------------Wide format--------------------
df_trials = data.frame(valence, session_id)     # each row is a trial
df_trials = left_join(df_trials, sub_session_id, by = "session_id") 

df_trials = df_trials %>% group_by(session_id,sub_id_num) %>%
  do(BlockAdd(.))       # combine blocks info

df_trials = df_trials %>% group_by(session_id,sub_id_num,blocks) %>% 
  do(TrialAdd(.))       # combine trials info

# Attention! group_by will automatically rank the data based on the first parameter sub_id_num, so now df_trials have no longer the same rank
# therefore, put session_id in the first parameter instead of sub_id_num

#------------------long format----------------------
blocks_arr = rep(df_trials$blocks, each = n_timepoints)
trials_arr = rep(df_trials$trials, each = n_timepoints)
df = data.frame(df, blocks = blocks_arr, trials = trials_arr)  # combine blocks and trials info to long format df
df$blocks = blocks_arr
df$trials = trials_arr
df_d$blocks = blocks_arr
df_d$trials = trials_arr

# df = within(df, rm("blocks", "trials","blocks.1","trials.1"))   delete certain columns in dataframe

#------------------save data------------------------
path_long_csv = file.path(path, "Long_format.csv")
path_long_d_csv = file.path(path, "long_format_derivative.csv")
path_wide_csv = file.path(path, "Wide_format.csv")
write.csv(df, path_long_csv, row.names = TRUE)
write.csv(df_trials, path_wide_csv, row.names = TRUE)
write.csv(df_d, path_long_d_csv, row.names = TRUE)

path_long_RData = file.path(path, "Long_format.RData")
path_long_d_RData = file.path(path, "long_format_derivative.RData")
path_wide_RData = file.path(path, "Wide_format.RData")
save(df, file = path_long_RData)
save(df_d, file = path_long_d_RData)
save(df_trials, file = path_wide_RData)


#p = file.path(path, "sub_session_id.RData")    # will return error if the absolute path is too long
save(sub_session_id, file = "/Users/liangyingliu/Desktop/sub_session_id.RData")
save(df_d, file = "/Users/liangyingliu/Desktop/sub_session_id.RData")

#-----------Functions to Rearrange the data----------------------
BlockAdd = function(df)
{
  valence = df$valence
  start = valence[1]
  rows = length(valence)
  blocks = c()
  block_count = 1
  
  for(i in 1:rows)
  {
    if (valence[i] == start)
    {
      blocks[i] = block_count
    }
    else
    {
      block_count = block_count + 1
      blocks[i] = block_count
      start = valence[i]
    }
  }
  df = data.frame(df, as.data.frame(blocks))      # must dataframe
}

TrialAdd = function(df)
{
  n_trials = length(df$blocks)
  trials = seq(1, n_trials, 1)
  df = data.frame(df, as.data.frame(trials))      # must dataframe
}


#---------------------plot concentration time series to inspect shape similarity---------------------------
#--------pick subject and session and neurotransmitter-----------
sub = "010"
sub_name = str_c("sub_", sub)
session = filter(sub_session_id, sub_id_num == sub)$session_id
y = "DA"

#----- ---change facet labels---------------
valence_labs = c("Neg", "Pos")
names(valence_labs) = c("0", "1")

session_labs = c(str_c("session", session))
names(session_labs) = c(session)

blocks_labs = c("block1", "block2", "block3", "block4")
names(blocks_labs) = c("1", "2", "3", "4")

#----------plot with single trials------------------
plot_valence_all <- ggplot(df_d %>% filter(sub_id_num == sub), mapping = aes(x = Timepoint, y = DA, color = factor(trials)))+
  geom_line()+
  labs(title = sub_name, color = "Trials")+
  facet_wrap(~ valence + session_id + blocks, scales = "fixed", 
             labeller = labeller(valence = valence_labs, session_id = session_labs, blocks = blocks_labs))+
  geom_vline(aes(xintercept = 26),
             color = "black",
             alpha = 0.7)+
  theme_light()+
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 7)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))

plot_valence_all <- ggplot(df_d %>% filter(sub_id_num == sub), mapping = aes(x = Timepoint, y = DA, color = factor(valence)))+
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

#---------save plots to corresponding subject and session files-----------------
root_path = "/Users/liangyingliu/Desktop/Results"
sub_path = file.path(root_path, sub_name)
sub_session_path = file.path(sub_path, str_c("session_", session))

if(! dir.exists(sub_path))
{
  dir.create(sub_path)
}

if(! dir.exists(sub_session_path))
{
  dir.create(sub_session_path)
}

ggsave(plot = plot_valence_all, filename = str_c("DA_single_trial", ".PNG"), path = sub_session_path)
ggsave(plot = plot_valence_all, filename = "DA_Session_119_120_121", path = sub_session_path)


