
TFC_bar_data$Genotype <- as.factor(TFC_bar_data$Genotype)

aggregate(TFC_bar_data[, 3:6], list(TFC_bar_data$Genotype), mean)

#ggplot(data=TFC_bar_data, aes(x=Subject, y=TFC2_baseline, fill=Genotype)) +
#  geom_bar(stat="identity", fill="steelblue")+
#  theme_minimal()

ggplot(data = TFC_bar_data,
       aes(x = Genotype, y = TFC2_trace, fill = factor(Genotype))) +
  geom_rain()
  # stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
  #              geom = "point", shape = 18, size = 3,
  #              show.legend = FALSE)

