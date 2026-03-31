#' OFE permutation analysis
#'
#' Runs a cell-based OFE analysis with permutation ANOVA and spatial diagnostics.
#'
#' @param data sf points; must contain columns `y` and `x`.
#' @param y response column (numeric).
#' @param x treatment column (factor/character).
#' @param cellsize,shift,angle_deg,buffer,min_per_cell grid settings (used when `grid` is NULL).
#' @param n_p number of permutations per ANOVA run.
#' @param n_s number of sampling runs.
#' @param alpha significance threshold for letters.
#' @param p_adjust_method p-value adjustment method across pairwise comparisons.
#' @param keep_components `"none"`, `"light"`, or `"full"` to embed geometries in output.
#' @param grid optional: an `ofe_grid`, a list with `grid_sel`, or an `sf` object of polygons.
#' @param crs optional target projected CRS if `data` is in lon/lat.
#' @param nmin_cell,alpha_bonferroni **Deprecated.**
#'   Use `min_per_cell` instead of `nmin_cell`, and `p_adjust_method = "bonferroni"`
#'   instead of `alpha_bonferroni`.
#'
#' @return An object of class `ofemt_result`.
#' @export
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
#' @references A new method to compare
#' treatments in unreplicated on-farm experimentation. Córdoba M.,
#' Paccioretti P., Balzarini M. Under review.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#'   my_data <- ofe_f2
#'   res <- ofemt(my_data, y = "Yield", x = "Treatment",
#'                cellsize = 10, min_per_cell = 4, alpha = 0.05)
#' }
ofemt <- function(
  data,
  y,
  x,
  cellsize = 10,
  min_per_cell = 4,
  nmin_cell = NULL,
  n_p = 1000,
  n_s = 200,
  alpha = 0.05,
  shift = c(0, 0),
  p_adjust_method = c("none", "bonferroni", "holm", "BH"),
  crs = NULL,
  keep_components = c("none", "light", "full"),
  grid = NULL,
  angle_deg = 0,
  buffer = 0,
  alpha_bonferroni = NULL
) {
  # --- Deprecation handling ---
  if (!is.null(nmin_cell)) {
    warning(
      "`nmin_cell` is deprecated; use `min_per_cell` instead.",
      call. = FALSE
    )
    min_per_cell <- nmin_cell
  }

  if (!is.null(alpha_bonferroni)) {
    warning(
      "`alpha_bonferroni` is deprecated; use `p_adjust_method = 'bonferroni'` instead.",
      call. = FALSE
    )
    p_adjust_method <- "bonferroni"
  }

  keep_components <- match.arg(
    if (is.logical(keep_components)) {
      if (keep_components) "light" else "none"
    } else {
      keep_components
    },
    c("none", "light", "full")
  )
  p_adjust_method <- match.arg(
    if (is.logical(p_adjust_method)) {
      if (p_adjust_method) "bonferroni" else "none"
    } else {
      p_adjust_method
    },
    c("none", "bonferroni", "holm", "BH")
  )

  stopifnot(inherits(data, "sf"))
  stopifnot(
    length(y) == 1,
    length(x) == 1,
    y %in% names(data),
    x %in% names(data)
  )
  if (!is.numeric(data[[y]])) {
    stop("`y` must be numeric.")
  }
  if (!(is.character(data[[x]]) || is.factor(data[[x]]))) {
    stop("`x` must be factor/character.")
  }

  if (sf::st_is_longlat(data)) {
    if (is.null(crs)) {
      stop("`data` is lon/lat. Provide a projected `crs`.")
    }
    crs_obj <- sf::st_crs(crs)
    if (sf::st_is_longlat(crs_obj)) {
      stop("Target `crs` must be projected.")
    }
    data <- sf::st_transform(data, crs_obj)
  }

  data[[".trt"]] <- gsub(" ", ".", as.character(data[[x]]))
  treatments <- length(unique(data$.trt))
  if (treatments < 2) {
    stop("Need at least two treatments.")
  }

  grid_source <- "generated"
  used_params <- list(
    cellsize = if (length(cellsize) == 1) {
      c(cellsize, cellsize)
    } else {
      cellsize[1:2]
    },
    shift = shift,
    angle_deg = angle_deg,
    buffer = buffer,
    min_per_cell = min_per_cell
  )
  # min_applied <- TRUE
  min_used <- min_per_cell
  info_note <- NA_character_

  # --- grid handling

  if (!is.null(grid) & !inherits(grid, "ofe_grid")) {
    stop("`grid` must be ofe_grid.")
  }

  if (!is.null(grid) & inherits(grid, "ofe_grid")) {
    stopifnot(inherits(grid$grid_sel, "sf"))
    grid_source <- "provided_ofe_grid_sel"
    used_params <- grid$params
    grid <- grid$grid_sel
    # class(grid) <- c('ofe_grid', class(grid))
    # min_applied <- FALSE
    min_used <- used_params$min_per_cell
  }

  # --- generated path: call make_ofe_grid and DO NOT re-apply min threshold here
  if (grid_source == "generated") {
    gint <- make_ofe_grid(
      data = data,
      x = x,
      cellsize = used_params$cellsize,
      angle_deg = used_params$angle_deg,
      buffer = used_params$buffer,
      shift = used_params$shift #,
      # min_per_cell = used_params$min_per_cell,
      # return_both = TRUE
    )
    grid <- gint$grid_sel
    ap <- attr(grid, "ofe_params", exact = TRUE)
    if (!is.null(ap)) {
      used_params <- ap
    }
    # min_applied <- FALSE
    min_used <- used_params$min_per_cell
  } else {
    if (is.na(sf::st_crs(grid))) {
      stop("CRS mismatch between `grid` and `data`.")
    }
    if (!any(sf::st_geometry_type(grid) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("`grid` must be polygons.")
    }
    if (sf::st_crs(grid) != sf::st_crs(data)) {
      warning("`grid` CRS was trasformed to `data` CRS.", call. = FALSE)
      grid <- sf::st_transform(grid, sf::st_crs(data))
    }
  }

  # if (!"CellID" %in% names(grid)) {
  #   grid$CellID <- seq_len(nrow(grid))
  # }

  # # join points to cells
  jdat <- sf::st_join(data, grid, join = sf::st_intersects, left = FALSE)

  my_unique_n <-
    stats::aggregate(
      datos$gvar,
      by = list("IDcelda" = datos$IDcelda),
      function(x) {
        c(
          length(unique(x)),
          ifelse(length(unique(x)) == 1, unique(x), NA),
          length(x)
        )
      },
      simplify = TRUE
    )
  my_unique_n <- do.call(data.frame, my_unique_n)
  my_unique_n <-
    stats::setNames(my_unique_n, c("IDcelda", "unique_trat", "trat", "n"))

  my_row_to_keep <-
    my_unique_n$unique_trat == 1 & my_unique_n$n >= nmin_cell

  my_IDcelda_keeped <- my_unique_n[my_row_to_keep, ]

  # if all are FALSE (!TRUE)
  if (all(!my_row_to_keep)) {
    stop(paste0(
      'No cells with more than nmin_cell (',
      nmin_cell,
      ') observations.'
    ))
  }

  datos_celda_sel <- subset(
    datos,
    datos[["IDcelda"]] %in% my_IDcelda_keeped[["IDcelda"]]
  )

  compar <- utils::combn(unique(data$gvar), 2, simplify = TRUE)

  variables <- c(y, "X", "Y")

  datos_celda_sel_mediana <-
    data.frame(
      sf::st_coordinates(datos_celda_sel),
      sf::st_drop_geometry(datos_celda_sel)
    ) %>%
    dplyr::group_by(gvar, IDcelda) %>%
    dplyr::summarise_at(variables, stats::median) %>%
    dplyr::ungroup()

  datos_celda_sel_mediana <-
    dplyr::left_join(
      datos_celda_sel_mediana,
      unique(sf::st_drop_geometry(datos_celda_sel[, c("IDcelda", "gvar", x)])),
      by = c("IDcelda", "gvar")
    )

  my_model <- stats::lm(
    stats::as.formula(paste(y, paste(x, collapse = " + "), sep = " ~ ")),
    data = datos_celda_sel_mediana
  )
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
    tryCatch(
      {
        spdep::nb2listw(gri, glist = fdist, style = "W")
      },
      error = function(e) {
        spdep::set.ZeroPolicyOption(TRUE)
        spdep::nb2listw(gri, glist = fdist, style = "W")
      }
    )

  # Magnitud de la autocorrelacion espacial de los residuos
  rho <- spatialreg::aple(datos_celda_sel_mediana$residuos, lw)
  # rho <-  ifelse(rho < 0, 0, rho)
  # Constrain to [0, 1]
  rho <- min(1, max(0, rho))
  MI <-
    spdep::moran.test(datos_celda_sel_mediana$residuos, lw)[[3]][[1]]
  n <- nrow(datos_celda_sel_mediana)
  my_n_eff <- n_eff(rho = rho, n = n)
  ess <- round(my_n_eff, 0)

  tablamuestra <- data.frame(
    "n" = n,
    "ESS" = ess,
    "Rho" = rho,
    "MI" = MI
  )

  # Medias Tratamientos
  medias_tratamientos <-
    sf::st_drop_geometry(datos_celda_sel_mediana) %>%
    dplyr::group_by_at(x) %>%
    dplyr::summarize_at(y, stats::median)

  # Permutacion multiple
  multipermutacion <- function(p) {
    # p=1
    base_permutacion <- datos_celda_sel_mediana %>%
      dplyr::group_by(gvar) %>%
      dplyr::slice_sample(n = ceiling(ess / treatments))

    permt_trat <- function(x_compar) {
      #x=1
      tabla_permt_parcial <-
        base_permutacion[base_permutacion$gvar %in% compar[, x_compar], ]
      suppressWarnings(
        tabla_permt_parcial_2 <-
          permuco::aovperm(
            stats::as.formula(paste(y, x, sep = "~")),
            data = tabla_permt_parcial,
            np = n_p
          )
      )

      tabla_permt_parcial_2 <-
        data.frame(
          'Trt_1' = compar[, x_compar][1],
          'Trt_2' = compar[, x_compar][2],
          "Comparison" = paste(
            unique(tabla_permt_parcial$gvar),
            collapse = " vs. "
          ),
          "p-valor" = tabla_permt_parcial_2$table$`resampled P(>F)`[1]
        )
    }

    tabla_permutacion_multiple <-
      do.call(rbind, lapply(seq_len(ncol(compar)), permt_trat))
    tabla_permutacion_multiple$Corrida <- p
    tabla_permutacion_multiple
  }
  resultadosmultipermutacion <-
    withr::with_seed(7, do.call(rbind, lapply(seq_len(n_s), multipermutacion)))

  # Medias p valor multi permutaciones
  medias_resultadosmultipermutacion <-
    resultadosmultipermutacion %>%
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
    medias_resultadosmultipermutacion[
      order(medias_resultadosmultipermutacion[[y]]),
    ]
  diffpermutacion <- medias_resultadosmultipermutacion$`p-value`
  names(diffpermutacion) <-
    gsub(" vs. ", "-", medias_resultadosmultipermutacion$Comparison)
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
      'letters' = my_letters,
      check.names = FALSE
    )

  colnames(letras)[colnames(letras) == 'Treatment'] <- x
  diffpermutacion_letras <- merge(medias_tratamientos, letras)
  diffpermutacion_letras <-
    diffpermutacion_letras[
      order(diffpermutacion_letras[[y]], decreasing = TRUE),
    ]

  colnames(diffpermutacion_letras)[colnames(diffpermutacion_letras) == y] <-
    paste(y, 'mean', sep = '_')

  infogral <-
    data.frame(
      "Cellsize" = paste(cellsize, collapse = "_"),
      "Min.Obs/Cell" = nmin_cell,
      "Max.Obs/Cell" = max(as.numeric(my_IDcelda_keeped$n)),
      "Median.Obs/Cell" = median(as.numeric(my_IDcelda_keeped$n)),
      tablamuestra
    )
  resultadotabla <-
    list(
      "General information" = infogral,
      "Cells per treatment" = table(my_IDcelda_keeped$trat),
      "ANOVA permutation test" = medias_resultadosmultipermutacion[, c(
        "Comparison",
        "p-value"
      )],
      "Means comparison" = diffpermutacion_letras
    )
  resultadotabla
}


#' Effective sample size
#'
#' @description An approximate calculation for the effective sample size for
#' spatially autocorrelated data. Only valid for approximately normally
#' distributed data.
#'
#' @param n Number of observations.
#' @param rho Spatial autocorrelation parameter from a simultaneous
#' autoregressive model.
#'
#' @return Returns effective sample size n*, a numeric value.
#'
#' @details
#'
#' Implements Equation 3 from Griffith (2005).
#'
#' @source
#'
#' Griffith, Daniel A. (2005). Effective geographic sample size in the presence of spatial autocorrelation. Annals of the Association of American Geographers. Vol. 95(4): 740-760.
#'
n_eff <- function(n, rho) {
  a <- 1 / (1 - exp(-1.92369))
  b <- (n - 1) / n
  c <- (1 - exp(-2.12373 * rho + 0.20024 * sqrt(rho)))
  n_eff <- n * (1 - a * b * c)
  return(n_eff)
}
