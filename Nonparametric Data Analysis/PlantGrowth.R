soil = read.csv('data/soil.csv')

ggplot(soil, aes(respiration, fill = condition)) + geom_histogram(bins = 20)
ggplot(soil[soil$condition == 'growth',], aes(y=respiration)) + geom_boxplot() + 
  labs(x = 'Soil Condition', y = 'Respiration (mols)', 
       title = 'Boxplot of CO2 in Soil for Growth Group')
ggsave('growthbox.png', width = 5, height = 10)

ggplot(soil[soil$condition == 'gap',], aes(y=respiration)) + geom_boxplot() + 
  labs(x = 'Soil Condition', y = 'Respiration (mols)', 
       title = 'Boxplot of CO2 in Soil for Gap Group') 
ggsave('gapbox.png', width = 5, height = 10)


med_growth = median(soil[soil$condition == 'growth',]$respiration)
med_gap = median(soil[soil$condition == 'gap',]$respiration)

d_obs = med_gap - med_growth

all.perms = sapply(1:10000,function(i){
  the.numbers = soil$respiration
  the.groups = soil$condition
  change.groups = sample(the.groups,length(the.groups),replace = FALSE) # shuffles groups
  group.1.med =  median(the.numbers[change.groups == levels(the.groups)[1]]) # finds median for group 1
  group.2.med = median(the.numbers[change.groups == levels(the.groups)[2]]) # finds median for group 2
  difference.in.meds= group.1.med-group.2.med #finds difference in means
  return(difference.in.meds)
})

sample.meds = aggregate(respiration ~ condition, soil,median)[,2] # finds mean per group
difference = sample.meds[1] - sample.meds[2] #finds difference in means
p.value.two = mean(abs(all.perms) >= abs(difference)) #calculates two-sided p-value
p.value.two 
p_val_ci = round(c(p.value.two - 1.96 * sqrt(p.value.two * (1 - p.value.two) / 2000), 
                   p.value.two + 1.96 * sqrt(p.value.two * (1 - p.value.two) / 2000)), 5)
