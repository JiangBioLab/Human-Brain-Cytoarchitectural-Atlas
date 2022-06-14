library(ggplot2)

meta <- readRDS('seu.harmony.anno.meta.v2.rds')

aa <- table(meta$hicat_cluster_subclasses, meta$region)

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
region <- c('FPPFC','DLPFC','VLPFC','M1','S1','S1E','PoCG','SPL','SMG','AG','VISP','ITG','STG','ACC')
regionID <- c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
names(region_color) <- region
names(region) <- regionID

library(ggmosaic)
library(wesanderson)

region_color <- c('#3F4587','#8562AA','#EC8561','#B97CB5','#D43046','#F0592B','#ED4A96','#593C97','#A54486','#FBDE13','#299FAE','#75CCE3','#0C6939','#0D9547')
# subclasses_region.pdf 5x15
ggplot(data = meta) +
  geom_mosaic(aes(x = product(hicat_cluster_subclasses), fill = region), size = 0.3, color = 'black') +
  scale_fill_manual(values = region_color) +
  #scale_fill_manual(values = c(wes_palette("Darjeeling2"),wes_palette("Darjeeling2"),wes_palette("Darjeeling2"))) +
  theme_classic() + RotatedAxis()


Exc <- meta[meta$hicat_cluster_merge_level2 == 'EXC',]
Exc$hicat_cluster_merge_newID <- factor(Exc$hicat_cluster_merge_newID,
                                        levels = c("104","60", "116","140","119","111","80", "107","81",
                                                   "53", "63", "114","72", "142","124","115","137","61", "106","74", "47", "43", "22", "13", "41", "108","18", "30", "34", "65",
                                                   "62", "12", "14", "48", "39", "26", "7",  "17", "9",  "3",  "38", "10", "5",  "11", "16", "71", "133","98", "131","136","128",
                                                   "123","85", "139","93", "91", "138","134", "76"))
# Exc_cluster_region.pdf 5x18
ggplot(data = Exc) +
  geom_mosaic(aes(x = product(hicat_cluster_merge_newID), fill = region), size = 0.3, color = 'black') +
  scale_fill_manual(values = region_color) +
  #scale_fill_manual(values = c(wes_palette("Darjeeling2"),wes_palette("Darjeeling2"),wes_palette("Darjeeling2"))) +
  theme_classic() + RotatedAxis()


Inh <- meta[meta$hicat_cluster_merge_level2 == 'INH',]
Inh$hicat_cluster_merge_newID <- factor(Inh$hicat_cluster_merge_newID,
                                        levels = c("83", "19", "109","32", "33", "86", "54", "75", "73", "122","120","87", "101",
                                                   "45", "118","126","130","40", "90", "99", "94", "20", "44", "51", "49", "132",
                                                   "92", "78", "29", "52", "68", "46", "15", "77",
                                                   "121","100","97", "67", "25", "24", "102","55", "112","110","28", "37", "21", 
                                                   "58", "113","8",  "88", "36", "50", "64", "105",
                                                   "141","23", "135"))
# Inh_cluster_region.pdf 5x18
ggplot(data = Inh) +
  geom_mosaic(aes(x = product(hicat_cluster_merge_newID), fill = region), size = 0.3, color = 'black') +
  scale_fill_manual(values = region_color) +
  #scale_fill_manual(values = c(wes_palette("Darjeeling2"),wes_palette("Darjeeling2"),wes_palette("Darjeeling2"))) +
  theme_classic() + RotatedAxis()



Non <- meta[meta$hicat_cluster_merge_level1 == 'Non-neuron',]
Non$hicat_cluster_merge_newID <- factor(Non$hicat_cluster_merge_newID,
                                        levels = c("66", "4",  "89", "57", "95", "1",  "56", "127","42", "27", "2",  "69",
                                                   "70", "82", "6",  "84", "35", "125","79", "117","129","31", "103","96", "59"))
# NonNeuron_cluster_region.pdf 5x15
ggplot(data = Non) +
  geom_mosaic(aes(x = product(hicat_cluster_merge_newID), fill = region), size = 0.3, color = 'black') +
  scale_fill_manual(values = region_color) +
  #scale_fill_manual(values = c(wes_palette("Darjeeling2"),wes_palette("Darjeeling2"),wes_palette("Darjeeling2"))) +
  theme_classic() + RotatedAxis()





library(ggmosaic)
library(wesanderson)


# NOT RUN {
data(titanic)

ggplot(data = titanic) +
  geom_mosaic(aes(x = product(Class), fill = Survived))
# good practice: use the 'dependent' variable (or most important variable)
# as fill variable

ggplot(data = titanic) +
  geom_mosaic(aes(x = product(Class, Age), fill = Survived))

ggplot(data = titanic) +
  geom_mosaic(aes(x = product(Class), conds = product(Age), fill = Survived))

ggplot(data = titanic) +
  geom_mosaic(aes(x = product(Survived, Class), fill = Age))

# Just excluded for timing. Examples are included in testing to make sure they work
# }
# NOT RUN {
data(happy)

ggplot(data = happy) + geom_mosaic(aes(x = product(happy)), divider="hbar")

ggplot(data = happy) + geom_mosaic(aes(x = product(happy))) +
  coord_flip()

# weighting is important
ggplot(data = happy) +
  geom_mosaic(aes(weight=wtssall, x=product(happy)))

ggplot(data = happy) + geom_mosaic(aes(weight=wtssall, x=product(health), fill=happy)) +
  theme(axis.text.x=element_text(angle=35))

ggplot(data = happy) +
  geom_mosaic(aes(weight=wtssall, x=product(health), fill=happy), na.rm=TRUE)

ggplot(data = happy) +
  geom_mosaic(aes(weight=wtssall, x=product(health, sex, degree), fill=happy),
              na.rm=TRUE)

# here is where a bit more control over the spacing of the bars is helpful:
# set labels manually:
ggplot(data = happy) +
  geom_mosaic(aes(weight=wtssall, x=product(age), fill=happy), na.rm=TRUE, offset=0) +
  scale_x_productlist("Age", labels=c(17+1:72))

# thin out labels manually:
labels <- c(17+1:72)
labels[labels %% 5 != 0] <- ""
ggplot(data = happy) +
  geom_mosaic(aes(weight=wtssall, x=product(age), fill=happy), na.rm=TRUE, offset=0) +
  scale_x_productlist("Age", labels=labels)

ggplot(data = happy) +
  geom_mosaic(aes(weight=wtssall, x=product(age), fill=happy, conds = product(sex)),
              divider=mosaic("v"), na.rm=TRUE, offset=0.001) +
  scale_x_productlist("Age", labels=labels)

ggplot(data = happy) +
  geom_mosaic(aes(weight=wtssall, x=product(age), fill=happy), na.rm=TRUE, offset = 0) +
  facet_grid(sex~.) +
  scale_x_productlist("Age", labels=labels)

ggplot(data = happy) +
  geom_mosaic(aes(weight = wtssall, x = product(happy, finrela, health)),
              divider=mosaic("h"))

ggplot(data = happy) +
  geom_mosaic(aes(weight = wtssall, x = product(happy, finrela, health)), offset=.005)

# Spine example
ggplot(data = happy) +
  geom_mosaic(aes(weight = wtssall, x = product(health), fill = health)) +
  facet_grid(happy~.)
# }
# NOT RUN {
# end of don't run
# }