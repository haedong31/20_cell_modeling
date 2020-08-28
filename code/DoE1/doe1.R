## design of experiment -----
# maximum resolution minimum abberation design

# Ito
set.seed(7981)

fnames <- seq(1, 7) %>% as.character()
dgn <- FrF2(nruns = 128, nfactors = 7, factor.names = fnames)

dgn_info_tbl <- design.info(dgn)
run_order <- run.order(dgn)
write_csv(dgn, "./dgn_Ito_matrix.csv")

# import experiment result
res <- read_csv('doe_Ito_tau.csv', col_names = FALSE)
res <- res$X1

# log transformation
# log_res <- log(res)
# res <- res/10000
# resp <- sort(res, decreasing = TRUE)

# estimate effects 
dgn_res <- add.response(dgn, res)
my_anv <- lm(dgn_res)
summary(my_anv)
half_mes <- my_anv$coefficients[2:8]
mes <- 2*half_mes
me_names <- names(mes)
names(mes) <- str_sub(me_names, 1, nchar(me_names)-1)

# visualization
hnplot(mes, method = 'Lenth', half = TRUE, ID = TRUE, main = "Tau Ito")
MEPlot(dgn_res, select = c(1, 2, 3, 4, 6, 7), main = "Effects plot for tau Ito")
# DanielPlot(dgn_res, half = TRUE, alpha = 0.15, 
#            main = "Half Normal Plot")

# IKslow
set.seed(8760)
fnames <- seq(1, 8) %>% as.character()
dgn <- FrF2(nruns = 256, nfactors = 8, factor.names = fnames)

dgn_info_tbl <- design.info(dgn)
run_order <- run.order(dgn)
write_csv(dgn, "./dgn_IKslow_matrix.csv")

res <- read_csv('doe_IKslow_tau.csv', col_names = FALSE)
res <- res$X1

# estimate effects 
dgn_res <- add.response(dgn, res)
my_anv <- lm(dgn_res)
summary(my_anv)
half_mes <- my_anv$coefficients[2:9]
mes <- 2*half_mes
me_names <- names(mes)
names(mes) <- str_sub(me_names, 1, nchar(me_names)-1)

# visualization
hnplot(mes, method = 'Lenth', half = TRUE, ID = TRUE, main = "Tau IKslow")
MEPlot(dgn_res, select = seq(1,8), main = "Effects plot for tau IKslow")
