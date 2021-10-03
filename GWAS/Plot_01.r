#pos_dis
ggplot(data, aes(x=Pair, y=Distance.M., group=Chr)) + 
  geom_line() + ggtitle("Distance(M) of each pair of locus") + 
  facet_grid(Chr ~ .,scales="free") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),axis.title.x = element_text(size = 17),axis.title.y = element_text(size = 17))

#box_plot
ggplot(data=data, aes(x=Chr, y=Distance.M.,group=Chr))+
  geom_boxplot() + 
  theme_classic() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),axis.title.x = element_text(size = 17),axis.title.y = element_text(size = 17))
geom_dotplot(binaxis='y', stackdir='center', dotsize=1)


ggplot(data=data, aes(x=Region, y=weight,group=Region))+
  geom_boxplot() + 
  theme_classic() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),axis.title.x = element_text(size = 17),axis.title.y = element_text(size = 17))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
