# written by Liangying, 1/27/2023

library(dtw)
library(dtwclust)
library(tidyverse)
library(reshape)
library(BBmisc)
library(plyr)
library(TSdist)
library(readxl)
library(stats)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

#-----------------------------------------------DTW----------------------------------------------------------
# attention that dtw package only recognize wide format data, i.e., each row is a time series
df_DA = cbind(DA, df_trials)
df_DA = data.frame(DA, df_trials)

#--------pick subject and session-----------
sub = "014"
sub_name = str_c("sub_","014")
session = filter(sub_session_id, sub_id_num == sub)$session_id

df_DA_dtw = filter(df_DA, session == session, sub_id_num == sub)[, 1:n_timepoints]
df_DA_dtw_all = filter(df_DA, session == session, sub_id_num == sub)

df_DA_dtw = filter(df_DA, session == session, sub_id_num == sub)[, 26:36]
df_DA_dtw = DA[, 26:36]

# https://rdrr.io/cran/dtwclust/man/tsclust.html
# change the step pattern -- registering a custom distance with the "proxy" package and using it
dtwstep <- function(x, y, ...) {
  dtw::dtw(x, y, step.pattern = symmetric1,
           distance.only = TRUE, ...)$distance
}       # for the PAM, the cross-distance matrix is the only input to the DTW algorithm

dtwstep_asymmetric <- function(x, y, ...) {
  dtw::dtw(x, y, step.pattern = rabinerJuangStepPattern(6,"c"),
           distance.only = TRUE, ...)$distance
}       # for the PAM, the cross-distance matrix is the only input to the DTW algorithm

# Registering the function with 'proxy'
if (!pr_DB$entry_exists("dtwSTEP2")) # only register once, you cannot directly modify it with the same name
  proxy::pr_DB$set_entry(FUN = dtwstep, names=c("dtwSTEP2"),
                         loop = TRUE, type = "metric", distance = TRUE,
                         description = "DTW with different step patterns")

# Registering the function with 'proxy'
if (!pr_DB$entry_exists("dtwSTEP1"))
  proxy::pr_DB$set_entry(FUN = dtwstep, names=c("dtwSTEP1"),
                         loop = TRUE, type = "metric", distance = TRUE,
                         description = "DTW with different step patterns")

# Registering the function with 'proxy'
if (!pr_DB$entry_exists("dtwAsymmetric"))
  proxy::pr_DB$set_entry(FUN = dtwstep_asymmetric, names=c("dtwAsymmetric"),
                         loop = TRUE, type = "metric", distance = TRUE,
                         description = "DTW with different step patterns")

pr_DB$get_entry_names() 
pr_DB$delete_entry("dtwSTEP2")
summary(pr_DB)
summary(pr_DB,"long")

load("/Users/liangyingliu/Desktop/Results/sub_013/session_133/dtw_DA.RData")
load("/Users/liangyingliu/Desktop/Results/sub_014/session_158/dtw_DA.RData")

load("/Users/liangyingliu/Desktop/Results/dtw_DA_allsubs.RData")

clust_pam_DA = tsclust(df_DA_dtw, type = "partitional", k= 2L:9L, distance = "dtw", centroid = "pam")

clust_pam_DA = tsclust(DA, type = "partitional", k= 2L:9L, distance = "dtw", centroid = "pam")

clust_pam_DA = tsclust(df_DA_dtw, type = "partitional", k= 2L:9L, distance = "Euclidean", centroid = "pam")

clust_pam_DA = tsclust(df_DA_dtw, type = "partitional", k= 2L:9L, distance = "dtwAsymmetric", centroid = "pam")
tmp3 = clust_pam_DA@distmat

SSE = sapply(clust_pam_DA, function(cl) { sum(cl@cldist ^ 2) })

root_path = "/Users/liangyingliu/Desktop/Results"
sub_path = file.path(root_path, sub_name)
sub_session_path = file.path(sub_path, str_c("session_", session))
p = file.path(sub_session_path, "dtw_DA_Euclidean.RData")
save(clust_pam_DA, file = p)

save(clust_pam_DA, file = "/Users/liangyingliu/Desktop/Results/dtw_DA_allsubs.RData")

#plot(clust_pam_DA, type = "sc")
plot(clust_pam_DA[[1]], type = "sc")    # clust_pam[[1]], clust_pam[[2]]...代表不同的cluster
plot(clust_pam_DA, type = "sc", clus = 2L)
plot(clust_pam_DA, type = "series", clus = 2L)
plot(clust_pam_DA[[1]], type = "centroids")

p2 = file.path(sub_session_path, "Euclidean/cluster9.pdf")
pdf(p2)
plot(clust_pam_DA[[2]], type = "sc")    # clust_pam[[1]], clust_pam[[2]]...代表不同的cluster
dev.off()

#------------------------------------------Derivative DTW----------------------------------------------------------
Derivative = function(ts)
{
  dts = c()
  n = length(ts) 
  for (i in 2:(n - 1))
  {
    dts[i] = ((ts[i] - ts[i-1]) + ((ts[i+1] - ts[i-1]) / 2)) / 2     # First-order derivative
    if (i == 2)
    {
      dts[1] = dts[i]
    }
    if(i == (n - 1))
    {
      dts[n] = dts[i]
    }
  }
  return(t(dts))
}

dDA = apply(df_DA_dtw, MARGIN = 1, FUN = Derivative)

# apply()函数是最常用的代替for循环的函数，可以对矩阵、数据框、多维数组，按行或列进行循环计算，但是只能是行或列
# apply(x, MARGIN, FUN, ...)       MARGIN: 1代表行，2代表列
# 如果是按照因子分组统计的话，可以用tapply(data, sex, mean)
# 如果是list这种不规则数据结构呢？用lapply()，可以对list逐成分进行计算

#---------------------------------------Clustering evaluation------------------------------------------------------
n = 8    # try times
cvi_sum = c()
clust_pam = clust_pam_DA

for (i in 1:n)
{
  tmp = unname(cvi(clust_pam[[i]]))   # named num -> unname
  cvi_sum = rbind(cvi_sum, t(tmp))    # add rows to tibble by using rbind()
}

cvi_sum = as_tibble(cvi_sum)
colnames(cvi_sum) = c("Sil", "SF", "CH","DB","DBstar","D","COP")
cvi_sum = cvi_sum %>% mutate(k = 2:(n+1))
cvi_sum = as.data.frame(cvi_sum) # if you want to use melt, you have to use DataFrame instead of tibble, the output from dplyr is not just a data frame, it inherits from dataframe. It's a known issue with dplyr.
cvi_sum_arr = melt(cvi_sum, id.vars = "k", measure.vars = c("Sil", "SF", "CH","DB","DBstar","D","COP"), variable_name = "CVI", value.name = "value")

plot <- ggplot(data = cvi_sum_arr, mapping = aes(x = k, y = value))+
  geom_line(linewidth = 0.6)+
  facet_wrap(~ CVI, scale = "free")+
  #labs(x = "k")+
  theme_light()+
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 7)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
        legend.position = "none")

elbow = data.frame(SSE = SSE, k = c(2:(n+1)))
plot <- ggplot(data = elbow, mapping = aes(x = k, y = SSE))+
  geom_line(linewidth = 0.6)+
  #labs(x = "k")+
  theme_light()+
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 7)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
        legend.position = "none")

plot

ggsave(plot = plot, filename = "cvi_Euclidean.PNG", path = sub_session_path)

#-------------------------------Piecharts plots for each cluster---------------------------------------------------
index = t(cbind(df_DA_dtw[,0], cluster = clust_pam_DA[[6]]@cluster))
index     # 一共有1000行，假设前500个trials是negative, 后500个trials是positive

neg_trials = which(df_DA_dtw_all$valence == 0)
pos_trials = which(df_DA_dtw_all$valence == 1)

c1_neg = length(neg_trials[which(index[neg_trials] == 1)])
c2_neg = length(neg_trials[which(index[neg_trials] == 2)])
c1_pos = length(pos_trials[which(index[pos_trials] == 1)])
c2_pos = length(pos_trials[which(index[pos_trials] == 2)])

c3_neg = length(neg_trials[which(index[neg_trials] == 3)])
c3_pos = length(pos_trials[which(index[pos_trials] == 3)])

c1 = c(c1_neg, c1_pos)
c2 = c(c2_neg, c2_pos)
labs = c("neg", "pos")

c3 = c(c3_neg, c3_pos)

pie_data = tibble(group = c("Neg", "Pos"),
                  value = c(c1_neg, c1_pos))

pie_data = tibble(group = c("Neg", "Pos"),
                  value = c(c2_neg, c2_pos))

pie_data = tibble(group = c("Neg", "Pos"),
                  value = c(c3_neg, c3_pos))

# Compute the position of labels
pie_data <- pie_data %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(pie_data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop)


# Piecharts plot
plot <- ggplot(pie_data, aes(x = "", y = prop, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  #facet_wrap(~ Cluster, scale = "fixed")+
  theme(legend.position = "none") +
  geom_text(aes(y = ypos, label = value), color = "white", size = 6) +
  scale_fill_manual(values = c("#868686FF", "#CD534CFF"))+
  labs(fill = "Valence", title = "Cluster1")+
  #guides(fill = guide_legend(title = "Valence"))+
  theme_void()+
  theme(legend.text = element_text(size = 17),
        legend.title = element_text(size = 17, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
#scale_fill_brewer(palette="Set1")
plot

ggsave(plot = plot, filename = "Cluster1.PNG", path = sub_session_path)
ggsave(plot = plot, filename = "elbow_plot.PNG", path = "/Users/liangyingliu/Desktop/Results/sub_013/session_133")
ggsave(plot = plot, filename = "elbow_plot.PNG", path = "/Users/liangyingliu/Desktop/Results/sub_014/session_158")

plot(clust_pam, type = "sc", clus = 2L)
plot(clust_pam, type = "series", clus = 1L)
plot(clust_pam, type = "centroids", clus = 1L)

clust_hier = tsclust(data_cluster, type = "h", k = 2L, distance = "dtw", trace = TRUE)   # preproc = zscore
plot(clust_hier, type = "sc")
plot(clust_hier)


#--------------------------------------------PCA plot------------------------------------------------------
load("/Users/liangyingliu/Desktop/Results/sub_013/session_133/dtw_DA.RData")
load("/Users/liangyingliu/Desktop/Results/sub_014/session_158/dtw_DA.RData")

index = t(cbind(df_DA_dtw[,0], cluster = clust_pam_DA[[8]]@cluster))

index = t(cbind(DA[,0], cluster = clust_pam_DA[[1]]@cluster))
cluster = t(index)

colnames(df_DA_dtw) = str_c("t", 1:57)
rownames = df_DA_dtw_all$valence
rownames[which(rownames == 0)] = "neg"
rownames[which(rownames == 1)] = "pos"

pca = prcomp(df_DA_dtw, center = TRUE, scale. = TRUE)
#summary(pca)

plot = ggbiplot(pca, ellipse = TRUE,groups = factor(cluster), var.axes = FALSE, labels = rownames)+
       labs(color = "Cluster", title = "Sub_014")+
       #scale_color_manual(values = c("#868686FF", "#CD534CFF"))+
       theme_light()+
       theme(axis.text.x = element_text(size = 10, face = "bold"),
            axis.text.y = element_text(size = 10, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 7)),
            axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
            plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10, face = "bold"))
plot

ggsave(plot = plot, filename = "pca_cluster9.PNG", path = "/Users/liangyingliu/Desktop/Results/sub_014/session_158")

ggsave(plot = plot, filename = "pca_cluster2.PNG", path = "/Users/liangyingliu/Desktop/Results/sub_013/session_133")




























