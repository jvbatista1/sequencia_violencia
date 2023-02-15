#write.csv(painel, "painel_v2.csv", row.names = F)

library(dplyr)
library(tidyverse)
library(TraMineR)

load("C:/Users/Victor/Documents/sequencia_violencia/painel_v7.RData")

colnames(painel)

w2016 <- painel %>% 
  filter(year == 2016)

w2017 <- painel %>% 
  filter(year == 2017)

painel <- full_join(w2016, w2017, by = "ID3", suffix = c("_2016", "_2017"))

painel_emo <- painel %>% 
  drop_na(vio_emo_2016, vio_emo_2017, vio_emo_12_2017)

painel_fis <- painel %>% 
  drop_na(vio_fis_2016, vio_fis_2017, vio_fis_12_2017)

painel_sex <- painel %>% 
  drop_na(vio_sex_2016, vio_sex_2017, vio_sex_12_2017)

########################################
### Painel com tr?s per?odos
########################################
#painel <- rename(painel,
#       "2016" = violence12_2016,
#       "2017" = violence12_2017,
#       "2019" = violence12_2019)

### creating sequence object
painel.labels <- c('sem ocorrencia','ocorrencia')
painel.scode <- c('0', '1')

painelseq_emo_2 <- seqdef(painel_emo, c('vio_emo_2016', 'vio_emo_12_2017'), states = painel.scode, labels = painel.labels, xtstep = 1)
painelseq_emo_3 <- seqdef(painel_emo, c('vio_emo_2016', 'vio_emo_12_2017', 'vio_emo_2017'), states = painel.scode, labels = painel.labels, xtstep = 1)

painelseq_fis_2 <- seqdef(painel_fis, c('vio_fis_2016', 'vio_fis_12_2017'), states = painel.scode, labels = painel.labels, xtstep = 1)
painelseq_fis_3 <- seqdef(painel_fis, c('vio_fis_2016', 'vio_fis_12_2017', 'vio_fis_2017'), states = painel.scode, labels = painel.labels, xtstep = 1)

painelseq_sex_2 <- seqdef(painel_sex, c('vio_sex_2016', 'vio_sex_12_2017'), states = painel.scode, labels = painel.labels, xtstep = 1)
painelseq_sex_3 <- seqdef(painel_sex, c('vio_sex_2016', 'vio_sex_12_2017', 'vio_sex_2017'), states = painel.scode, labels = painel.labels, xtstep = 1)

### first 10 sequences
seqiplot(painelseq_emo_2, withlegend = T, title = 'Index plot (10 first sequences)', border = NA)

### 10 most frequent sequences
seqfplot(painelseq_emo_2, withlegend = T, border = NA, title = 'Sequence frequency plot')

seqtab(painelseq_emo_2, format = 'STS', tlim = 1:27)
seqtab(painelseq_emo_3, format = 'STS', tlim = 1:27)

seqtab(painelseq_fis_2, format = 'STS', tlim = 1:27)
seqtab(painelseq_fis_3, format = 'STS', tlim = 1:27)

seqtab(painelseq_sex_2, format = 'STS', tlim = 1:27)
seqtab(painelseq_sex_3, format = 'STS', tlim = 1:27)

### state distribution by time points
seqdplot(painel.seq, withlegend = F, border = NA, title = 'State distribution plot')

### entropy in each time point
seqHtplot(painel.seq, title = 'Entropy index')

### turbulence
Turbulence <- seqST(painel.seq)
summary(Turbulence)
hist(Turbulence, col = "cyan", mais = 'Sequence turbulence')

### optimal matching distance
submat <- seqsubm(painel.seq, method = "TRATE")
dist.om1 <- seqdist(painel.seq, method = "OM", indel = 1, sm = submat, with.missing = T)

### typology of the trajectories
library(cluster)
clusterward1 <- agnes(dist.om1, diss = T, method = 'ward')
plot(clusterward1)

### tr?s clusters
cli.3 <- cutree(clusterward1, k = 3)
cli.3 <- factor(cli.3, labels = c("Tipo 1", "Tipo 2", "Tipo 3"))
#table(cli.3fac)
seqfplot(painel.seq, group = cli.3, pbarw = T)
seqmtplot(painel.seq, group = cli.3)

### quatro clusters
cli.4 <- cutree(clusterward1, k = 4)
cli.4 <- factor(cli.4, labels = c("Tipo 1", "Tipo 2", "Tipo 3", "Tipo 4"))
#table(cli.4fac)
seqfplot(painel.seq, group = cli.4, pbarw = T)
seqmtplot(painel.seq, group = cli.4)
seqHtplot(painel.seq, group = cli.4, title = 'Entropy index')

### tr?s clusters definido pelo autor
seqfplot(painel.seq, group = painel$tipo, pbarw = T)
seqmtplot(painel.seq, group = painel$tipo)
seqHtplot(painel.seq, group = painel$tipo, title = 'Entropy index')

########################################
### Logistic regression
########################################
library(stats)
library(mfx)
### criando dummy
tipo_1 <- painel$tipo == 'tipo_1'
tipo_2 <- painel$tipo == 'tipo_2'
tipo_3 <- painel$tipo == 'tipo_3'

tipo_1_reglog <- glm(tipo_1 ~ idade + raca + religiao + escolaridade, family = binomial(link = logit), data = painel)
summary(tipo_1_reglog)
logitmfx(tipo_1_reglog, data = painel)

tipo_2_reglog <- glm(tipo_2 ~ idade + raca + religiao + escolaridade, family = binomial(link = logit), data = painel)
summary(tipo_2_reglog)

tipo_3_reglog <- glm(tipo_3 ~ idade + raca + religiao + escolaridade, family = binomial(link = logit), data = painel)
summary(tipo_3_reglog)

########################################
### Multinomial Logit Model
########################################

mlogit(tipo ~ idade + raca + religiao + escolaridade, data = painel)


########################################
### Multinomial Logit Model Tipo 2
########################################
table(painel$tipo)
table(painel$idade)
table(painel$raca)
table(painel$escolaridade)
table(painel$religiao)

painel$raca <- factor(painel$raca, levels = c(1, 2, 3, 5))
painel$religiao <- factor(painel$religiao, levels = c(1, 2, 3, 4, 5, 6))
painel$escolaridade <- factor(painel$escolaridade, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9))
painel$tipo <- factor(painel$tipo)

painel$tipo <- relevel(painel$tipo, ref = 'tipo_1')

library(nnet)
library(stargazer)

model <- multinom(tipo ~ idade + raca + religiao + escolaridade, data=painel)
stargazer(model, type="html", out="model.htm")

####### ocultar #######
#zvalues <- summary(model)$coefficients / summary(model)$standard.errors
#zvalues

#pvalues <- (1 - pnorm(abs(zvalues), 0, 1)) * 2
#pnorm(abs(zvalues), lower.tail=FALSE)*2
#pvalues

library(margins)
margins(model, category = 'tipo')

####### relative risk ratios #######
model_rrr <-exp(coef(model))
stargazer(model, type="text", coef=list(model_rrr), p.auto=F, out="modelrrr.htm")

logitmfx(model, data = painel)

########################################
### Painel com dois per?odos
########################################
### creating sequence object
########################################
painel.labels <- c('sem ocorr?ncia','ocorr?ncia')
painel.scode <- c('0', '1')
painel.seq <- seqdef(painel, c('violence12_2016', 'violence12_2017'), states = painel.scode, labels = painel.labels, xtstep = 1)

### first 10 sequences
seqiplot(painel.seq, withlegend = T, title = 'Index plot (10 first sequences)', border = NA)

### 10 most frequent sequences
seqfplot(painel.seq, withlegend = T, border = NA, title = 'Sequence frequency plot')
seqtab(painel.seq, format = 'STS', tlim = 1:27)

### state distribution by time points
seqdplot(painel.seq, withlegend = F, border = NA, title = 'State distribution plot')

####### clusteriza?oes alternativas #######
##### 2 per?odos #####
#painel_2p <- read.csv("~/pessoal/traminer/painel_v10_2p.csv", sep=";")
#painel.seq <- seqdef(painel_2p, c('violence12_2016', 'violence12_2017'), states = painel.scode, labels = painel.labels, xtstep = 1)
#painel.seq <- seqdef(filter(painel_2p, qChange_partner_2019 == 1), c('violence12_2016', 'violence12_2017'), states = painel.scode, labels = painel.labels, xtstep = 1)

### first 10 sequences
seqiplot(painel.seq, withlegend = T, title = 'Index plot (10 first sequences)', border = NA)

### 10 most frequent sequences
seqfplot(painel.seq, withlegend = T, border = NA, title = 'Sequence frequency plot')
seqtab(painel.seq, format = 'STS', tlim = 1:27)

### state distribution by time points
seqdplot(painel.seq, withlegend = F, border = NA, title = 'State distribution plot')

### optimal matching distance
submat <- seqsubm(painel.seq, method = "TRATE")
dist.om1 <- seqdist(painel.seq, method = "OM", indel = 1, sm = submat, with.missing = T)

### typology of the trajectories
#library(cluster)
clusterward1 <- agnes(dist.om1, diss = T, method = 'ward')
plot(clusterward1)

### tr?s clusters
cli.3 <- cutree(clusterward1, k = 3)
cli.3 <- factor(cli.3, labels = c("Tipo 1", "Tipo 2", "Tipo 3"))

#table(cli.4fac)
seqfplot(painel.seq, group = cli.3, pbarw = T)
#seqmtplot(painel.seq, group = cli.3)

##### 3 per?odos #####
#painel_3p <- read.csv("~/pessoal/traminer/painel_v10_3p.csv", sep=";")
#painel.seq <- seqdef(painel_3p, c('violence12_2016', 'violence12_2017', 'violence12_2019'), states = painel.scode, labels = painel.labels, xtstep = 1)
painel.seq <- seqdef(filter(painel_3p, qChange_partner_2019 == 1), c('violence12_2016', 'violence12_2017', 'violence12_2019'), states = painel.scode, labels = painel.labels, xtstep = 1)

### first 10 sequences
seqiplot(painel.seq, withlegend = T, title = 'Index plot (10 first sequences)', border = NA)

### 10 most frequent sequences
seqfplot(painel.seq, withlegend = T, border = NA, title = 'Sequence frequency plot')
seqtab(painel.seq, format = 'STS', tlim = 1:27)

### state distribution by time points
seqdplot(painel.seq, withlegend = F, border = NA, title = 'State distribution plot')

### optimal matching distance
submat <- seqsubm(painel.seq, method = "TRATE")
dist.om1 <- seqdist(painel.seq, method = "OM", indel = 1, sm = submat, with.missing = T)

### typology of the trajectories
#library(cluster)
clusterward1 <- agnes(dist.om1, diss = T, method = 'ward')
plot(clusterward1)

### tr?s clusters
cli.3 <- cutree(clusterward1, k = 3)
cli.3 <- factor(cli.3, labels = c("Tipo 1", "Tipo 2", "Tipo 3"))

#table(cli.4fac)
seqfplot(painel.seq, group = cli.3, pbarw = T)
#seqmtplot(painel.seq, group = cli.3)



########################################
### TESTE DE MC NEMAR
########################################

painel2 <- painel %>% 
  select(violence12_2016, violence12_2017, dur_relac_2016, q_intro_2_2016, q201_a, q101f, q101e, q722_d_2016) %>% 
  na.omit()

table(painel2$violence12_2016, painel2$violence12_2017)

test <- mcnemar.test(table(painel2$violence12_2016, painel2$violence12_2017))
test
                                                                                                                                                                                                                      