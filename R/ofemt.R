#' Method to compare treatments in unreplicated agricultural field trials
#' with georeferenced data
#'
#' @description
#' Approach to statistically analyze unreplicated OFE to promote field-specific
#' inference of treatment effects. Statistical tools for spatial data are
#' coupled with permutation tests to determine the statistical significance
#' between treatment means.
#'
#' @param data an \code{sf} object with projected coordinates. If the object
#' geometries are in geographical coordinates, a valid \code{crs} must be specified.
#' @param y character of length 1; specifying \code{data} column name containing
#'   (numeric) response variable.
#' @param x character of length 1; specifying \code{data} column name containing
#'   (character or factor) classification variable (treatments to compare).
#' @param cellsize numeric of length 1 or 2 with target cell size. Argument
#'    passed to \code{\link[sf]{st_make_grid}}.
#' @param nmin_cell numeric of length 1; specifying minimum number of data points
#'    that accommodate each grid cell.
#' @param n_p numeric positive value specifying the number of permutations.
#' @param n_s numeric positive value specifying the number of random
#'     samples of grid cells to carry out perumational ANOVAs.
#' @param alpha numeric of length 1, alpha value (between 0 and 1) for meant test.
#' @param shift numeric of length 1 or 2; specifying shift from lower.
#'   left corner coordinates (x, y) of \code{data}.
#' @param alpha_bonferroni logical specifying if p-value must be adjusted using
#'   Bonferroni correction for multiple comparisons.
#' @param crs coordinate reference system: integer with the EPSG code,
#'   or character with proj4string to convert coordinates if \code{data} has
#'   longitude/latitude data.
#' @param ... additional argument for \code{\link[sf]{st_make_grid}}.
#'
#'
#' @details
#' The methodology involves:
#' \enumerate{
#'   \item calculation of effective sample size (ESS) given the underlying spatial
#' structure.
#'   \item ANOVA permutation test on a random sample of ESS.
#'   \item generation of the empirical distribution of p-values from repetition of
#' step two. The median of this empirical distribution is regarded as the
#' p-value associated with the non-treatment effect hypothesis.
#' }
#' The test can be easily extended to cover scenarios with more than two
#' treatments by employing pairwise comparisons of OFE across multiple
#' treatments, with p-values can be adjusted for multiplicity
#' using Bonferroni correction
#'
#' @return a list of length 4 with mean test comparisson results
#'
#' @references A proposed method to compare treatments in unreplicated
#' agricultural field trials with georeferenced data. Córdoba M.,
#' Paccioretti P.,Balzarini M. Under review.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' data("ofe_f2", package = "ofemeantest")
#' ofemt(data = ofe_f2,
#'       y = "Yield_tn",
#'       x = "Treatment")
#'
ofemt <- function(data,
                  y,
                  x,
                  cellsize = 10,
                  nmin_cell = 4,
                  n_p = 1000,
                  n_s = 200,
                  alpha = 0.05,
                  shift = 0,
                  alpha_bonferroni = FALSE,
                  crs = NULL,
                  ...) {
  # Test
  if (!inherits(data, "sf")) {
    stop("data must be an sf object")
  }
  if (sf::st_is_longlat(data) & !is.na(sf::st_crs(data))) {
    stopifnot(
      'coordinates provided in crs are longlat degrees' =
        !sf::st_is_longlat(crs),
      'crs must be provided' =
        !is.null(crs)
    )
    data <- sf::st_transform(data, crs = crs)

  }
  # y
  if (missing(y)) {
    stop('y must be a valid column name')
  }
  stopifnot('y must be a vald column name' =
              y %in% colnames(data))
  if (!is.numeric(data[[y]])) {
    stop('column', y, 'must be a numeric column.')
  }
  # x
  if (missing(x)) {
    stop('x must be a valid column name')
  }
  if (length(x) != 1) {
    stop(paste0('x must be a vector of length 1.'))
  }
  stopifnot('x must be a vald column name' =
              x %in% colnames(data))
  if (!(is.character(data[[x]]) | is.factor(data[[x]]))) {
    stop('column ', x, ' must be a character or factor column.')
  }

  # shift
  if (!is.numeric(shift)) {
    stop(paste0('shift (', utils::head(shift), ') must be a number.'))
  }
  if (length(shift) < 1 | length(shift) > 2) {
    stop(paste0('shift (', utils::head(shift), ') must be a vector of length 1 or 2.'))
  }
  # dim celda
  if (!is.numeric(cellsize)) {
    stop(paste0('cellsize (', utils::head(cellsize), ') must be a number.'))
  }
  #length(cellsize) mus be 1 or 2
  if (length(cellsize) < 1 | length(cellsize) > 2) {
    stop(paste0(
      'cellsize (',
      utils::head(cellsize),
      ') must be a vector of length 1 or 2.'
    ))
  }
  # bonferroni
  if (!is.logical(alpha_bonferroni) |
      length(alpha_bonferroni) != 1) {
    stop(paste0('alpha_bonferroni must be a logical vector of length 1.'))
  }
  # alpha

  if (length(alpha) != 1) {
    stop(paste0('alpha must be a vector of length 1.'))
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop('alpha (',
         alpha,
         ') must be a unique numeric value between 0 and 1.')
  }

  if (!is.numeric(n_p)) {
    stop(paste0('n_p (', utils::head(n_p), ') must be a number.'))
  }
  if (n_p <= 0) {
    stop(paste0('n_p (', n_p, ') must be a positive number.'))
  }

  raw_data <- data

  # Removes spaces in tratment variable
  data[['gvar']] <- gsub(" ", ".", data[[x]])

  treatments <- length(unique(data$gvar))

  if (treatments < 2) {
    stop('Two or more treatments must be compared')
  }

  grid <- sf::st_as_sf(
    sf::st_make_grid(
      data,
      cellsize = cellsize,
      flat_topped = TRUE,
      offset = sf::st_bbox(data)[c("xmin", "ymin")] -
        shift,
      # ...
    )
  )

  grid$IDcelda <- seq_len(nrow(grid))


  datos <-
    sf::st_join(data,
                grid,
                join = sf::st_intersects,
                left = FALSE)

  grid3 <-
    sf::st_join(grid,
                data,
                join = sf::st_intersects,
                left = FALSE)


  my_unique_n <- stats::aggregate(datos$gvar,
                           by = list("IDcelda" = datos$IDcelda),
                           function(x) {
                             c(length(unique(x)),
                               ifelse(length(unique(x)) == 1,
                                      unique(x),
                                      NA),
                               length(x))
                           }, simplify = TRUE)
  my_unique_n <- do.call(data.frame, my_unique_n)
  my_unique_n <-
    stats::setNames(my_unique_n, c("IDcelda", "unique_trat", "trat", "n"))


  my_row_to_keep <-
    my_unique_n$unique_trat == 1 & my_unique_n$n >= nmin_cell

  my_IDcelda_keeped <- my_unique_n[my_row_to_keep,]


  # if all are FALSE (!TRUE)
  if (all(!my_row_to_keep)) {
    stop(paste0(
      'No cells with more than nmin_cell (',
      nmin_cell,
      ') observations.'
    ))
  }


  datos_celda_sel <- subset(datos,
                            datos[["IDcelda"]] %in% my_IDcelda_keeped[["IDcelda"]])#datos#merge(datos, my_unique_n, by = "IDcelda")


  # if (plot_grilla == T) {
  #   # set.seed(semilla)
  #   plot1 <-
  #     #mapview::mapview(grid3, alpha.regions = 0, layer.name = "Grilla") +
  #     mapview::mapview(
  #       datos,
  #       #col.regions = RColorBrewer::brewer.pal(n = treatments, name = "Set3"),
  #       zcol = "gvar",
  #       cex = 4,
  #       alpha = 0,
  #       layer.name = "Data"
  #     )
  #
  #   # mapview::mapviewGetOption("vector.palette")
  #
  #   gridsel <-
  #     subset(grid3, grid3$IDcelda %in% datos_celda_sel$IDcelda)
  #
  #   plot2 <-  mapview::mapview(
  #     gridsel,
  #     color = "blue",
  #     cex = 4,
  #     alpha.regions = 0,
  #     layer.name = "Selected cells"
  #   )
  #
  #   print(plot1 + plot2)
  #
  # }


  #browser()



  compar <- utils::combn(unique(data$gvar), 2, simplify = TRUE)

  variables <- c(y, "X", "Y")

  # Que diferencia hay aquii!?!?!? ----
  datos_celda_sel_mediana <-
    data.frame(sf::st_coordinates(datos_celda_sel),
               sf::st_drop_geometry(datos_celda_sel))  %>%
    dplyr::group_by(gvar, IDcelda) %>%
    dplyr::summarise_at(variables, stats::median) %>%
    dplyr::ungroup()

  datos_celda_sel_mediana <-
    dplyr::left_join(datos_celda_sel_mediana,
                     unique(sf::st_drop_geometry(datos_celda_sel[, c("IDcelda", "gvar", x)])),
                     by = c("IDcelda", "gvar"))

  my_model <- stats::lm(stats::as.formula(paste(y,
                                  paste(x, collapse = " + "),
                                  sep = " ~ ")),
                 data = datos_celda_sel_mediana)
  datos_celda_sel_mediana$residuos <-
    stats::aov(my_model)[['residuals']]

  datos_celda_sel_mediana <-
    sf::st_as_sf(
      datos_celda_sel_mediana,
      coords = c("X", "Y"),
      crs = sf::st_crs(datos_celda_sel)
    )

  # Matriz de vecinos
  k1 <- spdep::knn2nb(spdep::knearneigh(datos_celda_sel_mediana))
  # Minima distancia para asegurarse de tener un dato vecino.
  dmax <- max(unlist(spdep::nbdists(k1, datos_celda_sel_mediana)))
  # Armado de grilla
  gri <- spdep::dnearneigh(datos_celda_sel_mediana, 0, dmax)
  dist <- spdep::nbdists(gri, datos_celda_sel_mediana)
  # peso de vecino ponderado por la distancia
  fdist <- lapply(dist, function(x) {
    1 / (x / dmax)
  })
  # Matriz de pesos espaciales
  lw <-
    tryCatch({
      spdep::nb2listw(gri, glist = fdist, style = "W")
    },
    error = function(e) {
      spdep::set.ZeroPolicyOption(TRUE)
      spdep::nb2listw(gri, glist = fdist, style = "W")
    })

  # Magnitud de la autocorrelacion espacial de los residuos
  rho <- spatialreg::aple(datos_celda_sel_mediana$residuos, lw)
  IM <-
    spdep::moran.test(datos_celda_sel_mediana$residuos, lw)[[3]][[1]]
  n <- nrow(datos_celda_sel_mediana)
  ess <-  round(geostan::n_eff(rho = rho, n = n), 0)

  tablamuestra <- data.frame(
    "n" = n,
    "ESS" = ess,
    "Rho" = rho,
    "IM" = IM
  )


  # Medias Tratamientos
  medias_tratamientos <-
    sf::st_drop_geometry(datos_celda_sel_mediana)  %>%
    dplyr::group_by_at(x) %>%
    dplyr::summarize_at(y, median)


  # Permutacion multiple
  multipermutacion <- function(p) {
    # p=1
    base_permutacion <- datos_celda_sel_mediana %>%
      dplyr::group_by(gvar) %>%
      dplyr::slice_sample(n = ceiling(ess / treatments))

    permt_trat <- function(x_compar) {
      #x=1
      tabla_permt_parcial <-
        base_permutacion[base_permutacion$gvar %in% compar[,x_compar],]
      suppressWarnings(
        tabla_permt_parcial_2 <-
          permuco::aovperm(stats::as.formula(paste(y, x, sep = "~")),
                           data = tabla_permt_parcial,
                           np = n_p)
      )

      tabla_permt_parcial_2 <-
        data.frame(
          'Trt_1' = compar[, x_compar][1],
          'Trt_2' = compar[, x_compar][2],
          "Comparison" = paste(unique(tabla_permt_parcial$gvar), collapse = " vs. "),
          "p-valor" = tabla_permt_parcial_2$table$`resampled P(>F)`[1]
        )

    }

    tabla_permutacion_multiple <-
      do.call(rbind, lapply(seq_len(ncol(compar)), permt_trat))
    tabla_permutacion_multiple$Corrida <- p
    tabla_permutacion_multiple

  }
  resultadosmultipermutacion <-
    withr::with_seed(7,
                     do.call(rbind, lapply(seq_len(n_s), multipermutacion)))


  #write.table(resultadosmultipermutacion, paste0(getwd(),"/",field,"Valores_p_",y,"_",x,".txt"))

  # # Histograma p valor
  # plotpvalor <-
  #   ggplot2::ggplot(
  #     resultadosmultipermutacion,
  #     ggplot2::aes(x = p.valor,
  #                  fill = Comparison,
  #                  colour = Comparison)
  #   ) +
  #   ggplot2::geom_histogram(alpha = 0.5, position = "identity") +
  #   ggplot2::labs(y = "Absolute Frecuency", x = "p-value") +
  #   ggplot2::facet_wrap(~ Comparison) +
  #   ggplot2::xlim(0, 1) +
  #   ggplot2::ylim(0, 50)



  # Medias p valor multi permutaciones
  medias_resultadosmultipermutacion <-
    resultadosmultipermutacion  %>%
    dplyr::group_by(Comparison) %>%
    dplyr::summarise(
      "p-value" = median(p.valor, na.rm = TRUE),
      "Trt_1" = unique(Trt_1),
      "Trt_2" = unique(Trt_2)
    )


  medias_resultadosmultipermutacion <-
    merge(
      medias_resultadosmultipermutacion,
      medias_tratamientos,
      by.x = 'Trt_1',
      by.y = x
    )
  medias_resultadosmultipermutacion <-
    medias_resultadosmultipermutacion[order(medias_resultadosmultipermutacion[[y]]),]
  diffpermutacion <- medias_resultadosmultipermutacion$`p-value`
  names(diffpermutacion) <-
    gsub(" vs. ",
         "-",
         medias_resultadosmultipermutacion$Comparison)
  #browser()
  if (alpha_bonferroni) {
    alpha <- ifelse(treatments > 2, alpha / treatments, alpha)
  }

  letras_comp <-
    multcompView::multcompLetters(diffpermutacion, threshold = alpha)

  my_letters <- letras_comp$Letters
  if ("monospacedLetters" %in% names(letras_comp)) {
    my_letters <- letras_comp$monospacedLetters
  }
  # monospacedLetters
  letras <-
    data.frame(
      'Treatment' = names(letras_comp$Letters),
      'letters' =  my_letters,
      check.names = FALSE
    )

  colnames(letras)[colnames(letras) == 'Treatment'] <- x
  diffpermutacion_letras <- merge(medias_tratamientos,
                                  letras)
  diffpermutacion_letras <-
    diffpermutacion_letras[order(diffpermutacion_letras[[y]], decreasing = TRUE), ]

  # plotcomparacion <-
  #   ggplot2::ggplot(diffpermutacion_letras, ggplot2::aes(x = reorder(.data[[x]], .data[[y]]), y =
  #                                                          .data[[y]])) +
  #   ggplot2::geom_bar(
  #     width = 0.5,
  #     stat = "identity",
  #     color = "black",
  #     position = ggplot2::position_dodge(width = 0.8),
  #   ) +
  #   ggplot2::geom_text(
  #     ggplot2::aes(
  #       x = reorder(.data[[x]], .data[[y]]),
  #       y = .data[[y]] ,
  #       label = as.matrix(letters)
  #     ),
  #     position = ggplot2::position_dodge(width = 0.8),
  #     vjust = -(0.5),
  #     size = 5
  #   ) +
  #   ggplot2::labs(x = xlab,
  #                 y = ylab) +
  #   ggplot2::theme(text = ggplot2::element_text(size = 15))

  # print(plotcomparacion)

  colnames(diffpermutacion_letras)[colnames(diffpermutacion_letras) == y] <-
    paste(y, 'mean', sep = '_')
  #browser()

  infogral <-
    data.frame(
      "Cellsize" = paste(cellsize, collapse = "_"),
      "Min.Obs/Cell" = nmin_cell,
      "Max.Obs/Cell" = max(my_IDcelda_keeped$n),
      "Median.Obs/Cell" = median(my_IDcelda_keeped$n),
      tablamuestra
    )
  resultadotabla <-
    list(
      "General information" = infogral,
      "Cells per treatment" = table(my_IDcelda_keeped$trat),
      "ANOVA permutation test" = medias_resultadosmultipermutacion[, c("Comparison", "p-value")],
      "Means comparison" = diffpermutacion_letras
    )
  resultadotabla
}
