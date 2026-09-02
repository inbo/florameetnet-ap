# R/functions.R
# Functions used by _targets.R. One function per pipeline step, so each
# step is independently testable and targets can cache it correctly.

# ---- Data ingestion --------------------------------------------------------

read_clustering <- function(path) {
  read.csv2(path) |>
    janitor::clean_names() |>
    tibble::as_tibble()
}

build_cyclus_lookup <- function(
    cycle_length = 5, n_cycles = 2,
    start_year = 2016
) {
  # cyclus 1: jaar 1-5 = 2016-2020
  # cyclus 2: jaar 1-5 = 2021-2025
  years <- start_year + (seq_len(cycle_length * n_cycles) - 1)
  tibble::tibble(
    cyclus = factor(rep(seq_len(n_cycles), each = cycle_length)),
    cyclus_jaar = factor(rep(seq_len(cycle_length), times = n_cycles)),
    jaar = years
  )
}

read_meetnetdesign_raw <- function(path, cyclus_lookup) {
  read.csv2(path) |>
    janitor::clean_names() |>
    tibble::as_tibble() |>
    dplyr::rename(cyclus_jaar = jaar) |>
    dplyr::mutate(
      cyclus_jaar = factor(as.numeric(cyclus_jaar)),
      group = factor(group)
    ) |>
    dplyr::inner_join(
      cyclus_lookup, by = "cyclus_jaar", relationship = "many-to-many"
    )
}

compute_popsizes <- function(clustering) {
  clustering |>
    dplyr::select(ifbl = kmhok, group) |>
    dplyr::count(group, name = "pop_size") |>
    dplyr::mutate(group = factor(group))
}

compute_samplesizes <- function(meetnetdesign_raw) {
  meetnetdesign_raw |>
    dplyr::count(group, cyclus, cyclus_jaar, name = "sample_size_panel") |>
    dplyr::mutate(
      sample_size_cyclus = sum(sample_size_panel),
      .by = c(group, cyclus)
    )
}

build_meetnetdesign <- function(meetnetdesign_raw, popsizes, samplesizes) {
  meetnetdesign_raw |>
    dplyr::inner_join(popsizes, by = dplyr::join_by(group)) |>
    dplyr::inner_join(samplesizes, by = dplyr::join_by(cyclus_jaar, group, cyclus)) |>
    dplyr::mutate(
      weight_cyclus = pop_size / sample_size_cyclus,
      weight_panel = pop_size / sample_size_panel
    )
}

# ---- Florabank occurrence data ---------------------------------------------

fetch_florabank_data <- function(
    starting_year = 2015,
    ifbl_resolution = "1km-by-1km",
    taxongroup = "Vaatplanten") {
  db_connectie <- inbodb::connect_inbo_dbase("D0152_00_Flora")
  on.exit(inbodb::dbDisconnect(db_connectie), add = TRUE)

  inbodb::get_florabank_taxon_ifbl_year(
    connection = db_connectie,
    starting_year = starting_year,
    ifbl_resolution = ifbl_resolution,
    taxongroup = taxongroup,
    collect = TRUE
  ) |>
    janitor::clean_names()
}

fetch_florabank_hokken <- function() {
  db_connectie <- inbodb::connect_inbo_dbase("D0152_00_Flora")
  on.exit(inbodb::dbDisconnect(db_connectie), add = TRUE)
  hokken <- tbl(db_connectie, "Hok") |>
    collect()
  return(hokken)
}

join_survey_data <- function(meetnetdesign, meetnet_data) {
  meetnetdesign |>
    dplyr::left_join(
      meetnet_data,
      by = dplyr::join_by(ifbl == hok, jaar == jaar)
    )
}

filter_survey_data <- function(
  svydata,
  ifbl_remove,
  min_n_taxa
) {
  svydata |>
    dplyr::filter(
      !(ifbl %in% ifbl_remove),
      n() >= min_n_taxa,
      .by = c(ifbl, cyclus_jaar, group, cyclus, jaar)
    )
}

# ---- Taxon clustering (subspecies/varieties -> main species) --------------

build_soortenlijst <- function(svydata_ap) {
  svydata_ap |>
    dplyr::select(
      parent_taxon_id, parent_taxoncode,
      parent_naam_wetenschappelijk, parent_naam_nederlands
    ) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(parent_naam_wetenschappelijk))
}

cluster_taxon_names <- function(soortenlijst) {
  rgbif::name_parse(soortenlijst$parent_naam_wetenschappelijk) |>
    tibble::as_tibble() |>
    dplyr::select(
      parent_naam_wetenschappelijk = scientificname,
      genusorabove, specificepithet
    ) |>
    dplyr::mutate(clustertaxon = paste(genusorabove, specificepithet))
}

build_svydata_clustered <- function(svydata_ap, names_clustered) {
  svydata_ap |>
    dplyr::select(
      group, ifbl, cyclus_jaar, cyclus, jaar, maand, pop_size,
      sample_size_panel, sample_size_cyclus, weight_cyclus, weight_panel,
      parent_naam_wetenschappelijk
    ) |>
    dplyr::left_join(
      names_clustered |>
        dplyr::select(parent_naam_wetenschappelijk, clustertaxon),
      by = "parent_naam_wetenschappelijk"
    ) |>
    dplyr::distinct(
      group, ifbl, cyclus_jaar, cyclus, jaar, maand, pop_size,
      sample_size_panel, sample_size_cyclus, weight_cyclus, weight_panel,
      clustertaxon
    )
}

build_svydata_wide <- function(svydata_ap_clustered) {
  svydata_ap_clustered |>
    dplyr::mutate(value = ifelse(!is.na(clustertaxon), 1, NA)) |>
    tidyr::pivot_wider(
      names_from = clustertaxon,
      values_from = value,
      values_fill = 0
    )
}

get_species_list <- function(svydata_ap_clustered) {
  sort(unique(svydata_ap_clustered$clustertaxon))
}

# ---- Survey designs ---------------------------------------------------------

build_survey_design <- function(svydata_wide, weights_var) {
  # fpc is intentionally omitted: only well-surveyed squares were used to
  # define the strata, so strata membership of squares that are not
  # well-surveyed is unknown.
  ReGenesees::e.svydesign(
    data = svydata_wide,
    ids = ~ifbl, # PSU is the site (handles correlation t and t+5)
    strata = ~group,
    weights = stats::as.formula(paste0("~", weights_var))
    # fpc = ~pop_size intentionally omitted
  )
}

# ---- Status & change: overall ----------------------------------------------

get_status <- function(design, x) {
  form_y <- stats::as.formula(paste0("~ `", x, "`"))
  df <- ReGenesees::svystatTM(
    design = design, estimator = "Mean", y = form_y,
    by = ~cyclus, conf.int = TRUE
  )
  names(df) <- c("cyclus", "mean", "se", "lcl", "ucl")
  df$species <- x
  df
}

get_all_status <- function(design, species) {
  purrr::map(species, \(x) get_status(design, x)) |> purrr::list_rbind()
}

get_change <- function(design, x) {
  form_model <- stats::as.formula(paste0("`", x, "` ~ cyclus"))
  mod_out <- ReGenesees::svystatB(
    design = design, model = form_model, conf.int = TRUE
  )
  # row 1 = Intercept, row 2 = cyclus2 dummy (the change estimate)
  change_row <- mod_out[2, ]
  data.frame(
    species = x,
    change = change_row[[1]],
    se = change_row[[2]],
    lcl = change_row[[3]],
    ucl = change_row[[4]]
  )
}

classify_abs_change <- function(df) {
  dplyr::mutate(
    df,
    abs_significance = dplyr::case_when(
      lcl <= 0 & ucl >= 0 ~ "niet significant",
      ucl < 0            ~ "significante daling",
      lcl > 0            ~ "significante stijging"
    ),
    abs_change_type = dplyr::case_when(
      round(change, 8) > 0 ~ "Stijging",
      round(change, 8) < 0 ~ "Daling",
      round(change, 8) == 0 ~ "Geen wijziging"
    )
  )
}

get_all_change <- function(design, species) {
  purrr::map(species, \(x) get_change(design, x)) |>
    purrr::list_rbind() |>
    classify_abs_change()
}

plot_change <- function(change_results) {

  p <- change_results |>
    dplyr::filter(abs_significance != "niet significant") |>
    dplyr::mutate(species = stats::reorder(species, change)) |>
    ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_pointrange(
      ggplot2::aes(
        x = change, xmin = lcl, xmax = ucl, y = species,
        colour = abs_significance
      )
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format()) +
    ggplot2::labs(x = "Absolute wijziging")

  return(p)
}

# ---- Status & change: by stratum -------------------------------------------

get_status_stratum <- function(design, x) {
  form_y <- stats::as.formula(paste0("~ `", x, "`"))
  df <- ReGenesees::svystatTM(
    design = design, estimator = "Mean", y = form_y,
    by = ~ group + cyclus, conf.int = TRUE
  )
  names(df)[1:6] <- c("group", "cyclus", "mean", "se", "lcl", "ucl")
  df$species <- x
  df
}

get_all_status_stratum <- function(design, species) {
  purrr::map(species, \(x) get_status_stratum(design, x)) |> purrr::list_rbind()
}

get_change_stratum <- function(design, x) {
  form_model <- stats::as.formula(
    paste0("`", x, "` ~ 0 + group + group:cyclus")
  )
  mod_out <- ReGenesees::svystatB(
    design = design, model = form_model, conf.int = TRUE
  )

  res <- as.data.frame(mod_out)
  res$term <- rownames(res)
  change_rows <- res |> dplyr::filter(grepl("cyclus2", term))

  data.frame(
    species = x,
    group = sub("^group(.*):cyclus2$", "\\1", change_rows$term),
    change = change_rows[[1]],
    se = change_rows[[2]],
    lcl = change_rows[[3]],
    ucl = change_rows[[4]]
  )
}

get_all_change_stratum <- function(design, species) {
  purrr::map(species, \(x) get_change_stratum(design, x)) |>
    purrr::list_rbind() |>
    classify_abs_change()
}

plot_change_stratum <- function(change_stratum_results) {
  p <- change_stratum_results |>
    dplyr::filter(
      abs_significance != "niet significant",
      !is.na(abs_significance),
      abs(change) > 1e-3
    ) |>
    tidyr::nest(.by = group) |>
    dplyr::mutate(
      data = purrr::map(
        data, \(x) x |> dplyr::mutate(species = reorder(species, change))
      )
    )

  p <- p |>
    dplyr::mutate(
      plots = purrr::map2(
        .x = data,
        .y = group,
        .f = \(x, y) {
          ggplot2::ggplot(data = x) +
            ggplot2::geom_vline(
              xintercept = 0, linetype = "dashed", color = "gray50"
            ) +
            ggplot2::geom_pointrange(
              ggplot2::aes(
                x = change, xmin = lcl, xmax = ucl, y = species,
                colour = abs_significance
              )
            ) +
            ggplot2::scale_x_continuous(labels = scales::percent_format()) +
            ggplot2::labs(
              title = y,
              x = "Absolute wijziging", y = "Soort"
            )
        }
      )
    )
  return(p)
}

# ---- Relative change (ratio, delta method on log scale) -------------------

delta_method_ci <- function(ratio, se_ratio) {
  lcl <- numeric(length(ratio))
  ucl <- numeric(length(ratio))
  for (i in seq_along(ratio)) {
    if (is.na(ratio[i]) || is.infinite(ratio[i])) {
      lcl[i] <- NA
      ucl[i] <- NA
    } else if (ratio[i] == 0) {
      # complete disappearance: ratio 0, se 0; log breaks here so we
      # assign the exact bounds manually
      lcl[i] <- 0
      ucl[i] <- 0
    } else {
      log_ratio <- log(ratio[i])
      se_log_ratio <- se_ratio[i] / ratio[i]
      lcl[i] <- exp(log_ratio - 1.96 * se_log_ratio)
      ucl[i] <- exp(log_ratio + 1.96 * se_log_ratio)
    }
  }
  list(lcl = lcl, ucl = ucl)
}

get_rel_change_stratum <- function(design, x) {
  des_temp <- design
  species_vec <- as.numeric(as.character(des_temp$variables[[x]]))
  cyc_vec <- as.character(des_temp$variables$cyclus)
  des_temp$variables$num_var <- as.numeric(species_vec * (cyc_vec == "2"))
  des_temp$variables$den_var <- as.numeric(species_vec * (cyc_vec == "1"))

  df <- ReGenesees::svystatR(
    design = des_temp, num = ~num_var, den = ~den_var,
    by = ~group, conf.int = TRUE
  )
  ratio <- df[[2]]
  se_ratio <- df[[3]]
  ci <- delta_method_ci(ratio, se_ratio)

  data.frame(
    species = x, group = df[[1]],
    rel_change = ratio - 1, se = se_ratio,
    lcl = ci$lcl - 1, ucl = ci$ucl - 1
  )
}

get_rel_change <- function(design, x) {
  des_temp <- design
  species_vec <- as.numeric(as.character(des_temp$variables[[x]]))
  cyc_vec <- as.character(des_temp$variables$cyclus)
  des_temp$variables$num_var <- as.numeric(species_vec * (cyc_vec == "2"))
  des_temp$variables$den_var <- as.numeric(species_vec * (cyc_vec == "1"))

  df <- ReGenesees::svystatR(
    design = des_temp, num = ~num_var, den = ~den_var, conf.int = TRUE
  )
  ratio <- df[["Ratio"]]
  se_ratio <- df[["SE"]]
  ci <- delta_method_ci(ratio, se_ratio)

  data.frame(
    species = x, rel_change = ratio - 1, se = se_ratio,
    lcl = ci$lcl - 1, ucl = ci$ucl - 1
  )
}

classify_rel_changes <- function(df) {
  dplyr::mutate(
    df,
    rel_significance = dplyr::case_when(
      # New colonizations
      is.infinite(rel_change)         ~ "significante stijging",
      # Exact same sites      lcl < 0 & ucl > 0 ~ "niet significant"
      rel_change == 0 & se == 0   ~ "niet significant",
      lcl <= 0 & ucl >= 0     ~ "niet significant",
      ucl < 0            ~ "significante daling",
      lcl > 0            ~ "significante stijging"
    ),
    rel_change_type = dplyr::case_when(
      rel_change == -1 & se == 0  ~ "Enkel aanwezig in eerste survey",
      is.na(rel_change)           ~ "Afwezig in beide surveys",
      is.infinite(rel_change)     ~ "Enkel aanwezig in tweede survey",
      rel_change == 0 & se > 0    ~ "Geen netto verandering, wel andere hokken",
      rel_change < 0              ~ "Daling",
      rel_change > 0              ~ "Stijging",
      rel_change == 0 & se == 0   ~
        "Aanwezig in beide surveys in exact dezelfde hokken",
      TRUE                        ~ "Andere gevallen"
    )
  )
}

# n_workers controls the mirai daemon pool used to parallelise across
# species. Daemons are spun up and torn down within the call so the target
# itself stays a plain, cacheable function from targets' point of view.
get_all_rel_change_stratum <- function(design, species, n_workers = 6) {
  mirai::daemons(n_workers)
  on.exit(mirai::daemons(0), add = TRUE)

  purrr::map(
    species,
    purrr::in_parallel(
      \(x) get_rel_change_stratum(design, x),
      design = design,
      get_rel_change_stratum = get_rel_change_stratum,
      delta_method_ci = delta_method_ci
    )
  ) |>
    purrr::list_rbind() |>
    classify_rel_changes()
}

get_all_rel_change <- function(design, species, n_workers = 6) {
  mirai::daemons(n_workers)
  on.exit(mirai::daemons(0), add = TRUE)

  purrr::map(
    species,
    purrr::in_parallel(
      \(x) get_rel_change(design, x),
      design = design,
      get_rel_change = get_rel_change,
      delta_method_ci = delta_method_ci
    )
  ) |>
    purrr::list_rbind() |>
    classify_rel_changes()
}

plot_rel_change_stratum <- function(rel_change_results_stratum) {

  p <- rel_change_results_stratum |>
    dplyr::filter(
      rel_significance != "niet significant",
      rel_change_type %in% c("Daling", "Stijging")
    ) |>
    dplyr::mutate(
      ratio = rel_change + 1, lcl_ratio = lcl + 1, ucl_ratio = ucl + 1
    ) |>
    tidyr::nest(.by = group) |>
    dplyr::mutate(
      data = purrr::map(data, \(x) x |> dplyr::mutate(species = reorder(species, rel_change)))
    )

  p <- p |>
    dplyr::mutate(
      plots = purrr::map2(
        .x = data,
        .y = group,
        .f = \(x, y) {
          ggplot2::ggplot(data = x) +
            ggplot2::geom_vline(
              xintercept = 1, linetype = "dashed", color = "gray50"
            ) +
            ggplot2::geom_pointrange(
              ggplot2::aes(
                x = ratio, xmin = lcl_ratio, xmax = ucl_ratio, y = species,
                colour = rel_significance)
            ) +
            ggplot2::scale_x_log10(
              name = "Relatieve wijziging (Index)",
              breaks = c(1/50, 1/20, 1/10, 1/2, 1/5, 1, 2, 5, 10, 20, 50),
              labels = c("x1/50", "x1/20", "x1/10", "x1/2", "x1/5", "x1",
                         "x2", "x5", "x10", "x20", "x50")
            ) +
            ggplot2::labs(title = y) +
            ggplot2::labs(x = "Relatieve wijziging", y = "Soort")
        }
      )
    )

  return(p)
}

plot_rel_change <- function(rel_change_results) {
  p <- rel_change_results |>
    dplyr::filter(
      rel_significance != "niet significant",
      rel_change_type %in% c("Daling", "Stijging")
    ) |>
    dplyr::mutate(
      ratio = rel_change + 1, lcl_ratio = lcl + 1, ucl_ratio = ucl + 1,
      species = stats::reorder(species, rel_change)
    ) |>
    ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    ggplot2::geom_pointrange(
      ggplot2::aes(x = ratio, xmin = lcl_ratio, xmax = ucl_ratio, y = species,
                   colour = rel_significance)
    ) +
    ggplot2::scale_x_log10(
      name = "Relatieve wijziging (Index)",
      breaks = c(1/50, 1/20, 1/10, 1/2, 1/5, 1, 2, 5, 10, 20, 50),
      labels = c("x1/50", "x1/20", "x1/10", "x1/2", "x1/5", "x1",
                 "x2", "x5", "x10", "x20", "x50")
    ) +
    ggplot2::labs(x = "Relatieve wijziging", y = "Soort")

  return(p)
}
