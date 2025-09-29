
library(ggplot2)
library(plotly)
library(plyr)
library(ggpubr)
library(grid)
library(ggthemes)
library(R.matlab)
library(latex2exp)



My_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 10),
  axis.title.y = element_text(size = 12),
  axis.text.y = element_text(size = 10),
  plot.title = element_text(hjust = 0.5, size = 14),
  legend.title=element_text(size=14), 
  legend.text=element_text(size=12)
  # panel.grid.major.x = element_blank(),
  #panel.grid.minor.x = element_blank()
)



col=c("black", "red2", "purple2")
shapes = c(16, 17, 15)


load("~/Conrank_noshift.RData")


nReps <- ncol(Err_SCP)
nns <- rep(nS_Seq , nReps)    
nns <- rep(nns, 3)  # 3 different Methods

Method <- c(rep("Cox-kmeans", nnS*nReps),
            rep("SCP", nnS*nReps), 
            rep("TL-SCP", nnS*nReps))
Method <- factor(Method, levels=c('Cox-kmeans', 'SCP', "TL-SCP"))
#######################
se <- c(c(Err_cox), c(Err_SCP),  c(Err_tSCP))
#se <- log(se)
nmi <- c(c(NMI_M_cox), c(NMI_M_SCP), c(NMI_M_tSCP))

nclusters <- c(c(k_cox), c(k_SCP), c(k_tSCP))

dt <- data.frame(Method, nns, se, nmi, nclusters)
dt$nns <- as.factor(dt$nns)

###############
###  Relative err 
##############
p1 <- ggplot(dt, aes(x=nns, y=se, color=Method))+ #, fill = Method)) +
  geom_boxplot( outlier.shape = NA)+
  labs(#title = "MSE of coefficients for method 3",
    y="relErr", x=expression(source~sample~size~n[S]))+ theme_bw() +
  scale_shape_manual(values=shapes)+
  scale_color_manual(values=col)+
  #scale_fill_manual(values=col)+
  #scale_fill_brewer(palette="Dark2")+
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16,
                               angle = 45, vjust = 0.5, hjust=0.5),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    legend.title=element_text(size=18), 
    legend.text=element_text(size=16)
  )


###############
###  NMI 
##############

p2 <- ggplot(dt, aes(x=nns, y=nmi, color=Method)) +
  geom_boxplot( outlier.shape = NA)+
  labs(y='NMI', x=expression(source~sample~size~n[S]))+ theme_bw() +
  scale_shape_manual(values=shapes)+
  scale_color_manual(values=col)+ #ylim(0.5, 3)+
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16,
                               angle = 45, vjust = 0.5, hjust=0.5),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    legend.title=element_text(size=18), 
    legend.text=element_text(size=16)
  )




ggarrange(p1, p2, nrow=1, ncol = 2, common.legend = TRUE, legend = "top")



