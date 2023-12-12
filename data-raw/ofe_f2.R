## code to prepare `ofe_f2` dataset goes here

ofe_f2 <- sf::read_sf('data-raw/F2_BD_P_2020.gpkg')
ofe_f2 <- ofe_f2[, c('Trat', 'Rinde_kg')]
ofe_f2[['Rinde_kg']] <- ofe_f2[['Rinde_kg']] / 1000
ofe_f2 <- ofe_f2[ofe_f2[['Rinde_kg']] < 9.8, ]
colnames(ofe_f2)[1:2] <- c('Treatment', 'Yield_tn')
ofe_f2$Treatment <- gsub("Aplicado", "Fertilized", ofe_f2$Treatment)
ofe_f2$Treatment <- gsub("Testigo", "Control", ofe_f2$Treatment)

usethis::use_data(ofe_f2, overwrite = TRUE)
