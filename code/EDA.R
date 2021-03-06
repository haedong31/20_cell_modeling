library(readxl)
library(tidyverse)


## current-densities (normalized) ----
path <- './MGAT1_Data/Final MGAT1KO Data for Hui/JMCC paper/Nav Currents/IV Recordings/EB docs 3-22-19/Mgat1KO-Control INa Vr - 3-22-19.xlsx'

# wild type (control)
normI_wt <- read_excel(path, sheet = 'Control I-V', range = cell_limits(c(1, 1), c(24, 22)))  # cellranger::as.cell_limits('V24')
colnames(normI_wt)[2] <- 'cell_ID'
normI_wt$`Normalized (pA/pF)`[3:5] <- 16.90
normI_wt$`Normalized (pA/pF)`[7] <- 16.90
normI_wt$`Normalized (pA/pF)`[9:11] <- 17.30
normI_wt$`Normalized (pA/pF)`[13:15] <- 17.40
normI_wt$`Normalized (pA/pF)`[17:19] <- 18.90
normI_wt$`Normalized (pA/pF)`[22:23] <- 18.70

normI_wt_2 <- normI_wt %>% 
  pivot_longer(colnames(normI_wt)[3:22], names_to = 'voltage', values_to = 'INa') %>% 
  arrange(cell_ID)
normI_wt_2$voltage <- normI_wt_2$voltage %>% as.numeric()

normI_wt_mean <- normI_wt_2 %>% 
  group_by(voltage) %>% 
  summarise(mean_INa = mean(INa))

# WT plottings
normI_wt_2 %>% 
  ggplot(aes(x = voltage, y = INa, group = cell_ID)) +
  geom_point(aes(color = cell_ID)) +
  geom_line(aes(color = cell_ID)) +
  xlab('Voltage (mV)') +
  ylab('Normalized INa (pA/pF) WT') +
  theme(legend.position = 'bottom')

normI_wt_mean %>% 
  ggplot(aes(x = voltage, y = mean_INa, group = 1)) +
  geom_line() +
  geom_point() +
  xlab('Voltage (mV)') +
  ylab('Mean Normalized INa (pA/pF) WT')

# MGAT1KO 
normI_ko <- read_excel(path, sheet = 'KO IV', range = cell_limits(c(1, 1), c(17, 21)))  # cellranger::as.cell_limits('U17')
normI_ko_2 <- normI_ko %>% 
  pivot_longer(colnames(normI_ko)[2:21], names_to = 'voltage', values_to = 'INa')
normI_ko_2$voltage <- normI_ko_2$voltage %>% as.numeric()

normI_ko_mean <- normI_ko_2 %>% 
  group_by(voltage) %>% 
  summarise(mean_INa = mean(INa))

# MGAT1KO plottings 
normI_ko_2 %>% 
  ggplot(aes(x = voltage, y = INa, group = Normalized)) +
  geom_point(aes(color = Normalized)) +
  geom_line(aes(color = Normalized)) +
  xlab('Voltage (mV)') +
  ylab('Normalized INa (pA/pF) MGAT1KO') +
  theme(legend.position = 'bottom')

normI_ko_mean %>% 
  ggplot(aes(x = voltage, y = mean_INa, group = 1)) +
  geom_point() + 
  geom_line() +
  xlab('Voltage (mV)') +
  ylab('Mean Normalized INa (pA/pF) MGAT1KO') 

normI_mean <- inner_join(normI_wt_mean, normI_ko_mean, "voltage")
colnames(normI_mean) <- c('voltage', 'WT', 'MGAT1KO')
normI_mean <- normI_mean %>% 
  pivot_longer(c('WT', 'MGAT1KO'), names_to = 'Group', values_to = 'density')
normI_mean %>% 
  ggplot(aes(x = voltage, y = density, group = Group)) +
  geom_point(aes(color = Group)) +
  geom_line(aes(color = Group)) +
  xlab('Voltage (mV)') +
  ylab('Density (pA/pF)')


## Nav Mgat1KO IV Data ----
path <- './MGAT1_Data/Final MGAT1KO Data for Hui/JMCC paper/Nav Currents/IV Recordings/Mgat1KO IV/01-31-2018 Na.xlsx'

# steady-state activation
mgat1ko_ivs <- read_excel(path, sheet = 'IVS Cell ', range = cell_limits(c(1, 1), c(21, 6)))
mgat1ko_ivs %>% 
  filter(mV <= -20) %>% 
  ggplot(aes(x = mV, y = `G/Gmax`, group = 1)) +
  geom_point() +
  geom_line()

# steady-state inactivation
mgat1ko_ssi <- read_excel(path, sheet = 'SSI Cell', range = cell_limits(c(1, 1), c(17, 3)))
mgat1ko_ssi %>% 
  ggplot(aes(x = mV, y = I2, group = 1)) +
  geom_line() +
  geom_point()
mgat1ko_ssi %>% 
  ggplot(aes(x = mV, y = `I2/I2max`, group = 1)) +
  geom_line() +
  geom_point()

# recovery
mgat1ko_rev <- read_excel(path, sheet = 'Recovery -100', range = cell_limits(c(1, 1), c(31, 4)))
mgat1ko_rev %>% 
  ggplot(aes(x = `Time mS`, y = I1, group = 1)) +
  geom_point() +
  geom_line()
mgat1ko_rev %>% 
  ggplot(aes(x = `Time mS`, y = I2, group = 1)) +
  geom_point() +
  geom_line()
mgat1ko_rev %>% 
  ggplot(aes(x = `Time mS`, y = `I2/I1`, group = 1)) +
  geom_point() +
  geom_line()


## Nav WT Data ----
path <- './MGAT1_Data/Final MGAT1KO Data for Hui/JMCC paper/Nav Currents/IV Recordings/WT IV/Nav_WT.xlsx'

ivs_wt <- read_excel(path, sheet = 'IVS')
ssi_wt <- read_excel(path, sheet = 'SSI')
recovery_wt <- read_excel(path, sheet = 'recovery-100')

# aggregate the data
ggmax_mean_wt <- ivs_wt %>% 
  group_by(voltage) %>% 
  summarise(mean_ggmax = mean(`G/Gmax`))

iimax_mean_wt <- ssi_wt %>% 
  group_by(voltage) %>% 
  summarise(mean_iimax = mean(`I2/I2max`))

# plotting
# activation & inactivation

ggmax_mean_wt %>%
  filter(voltage <= -20) %>% 
  ggplot(aes(x = voltage, y = mean_ggmax, group = 1)) +
  geom_point() +
  geom_line() +
  ylab('SS Activation')

# inactivation
iimax_mean_wt %>% 
  ggplot(aes(x = voltage, y = mean_iimax, group = 1)) +
  geom_point() +
  geom_line() +
  ylab('SS Inactivation')

