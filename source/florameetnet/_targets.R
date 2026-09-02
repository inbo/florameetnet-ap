library(targets)

tar_option_set(
  packages = c(
    "ReGenesees", "janitor", "inbodb", "dplyr", "tidyr", "rgbif",
    "purrr", "mirai", "tibble", "ggplot2", "scales", "tidytext",
    "rprojroot", "sf"
  ),
  format = "qs"
)

tar_source("_targets_functions.R")

inaccessible_ifbl <- c("b4-45-31", "c7-47-42")
richness_threshold <- 50

list(
  tar_target(
    git_root,
    find_root(is_git_root)
  ),
  tar_target(
    clustering_file,
    file.path(git_root, "data/fbl_hok_clustering.csv"),
    format = "file"
  ),
  tar_target(clustering, read_clustering(clustering_file)),
  tar_target(cyclus_lookup, build_cyclus_lookup()),
  tar_target(
    meetnetdesign_file,
    file.path(git_root, "data/meetnet_ap_hokken.csv"),
    format = "file"
  ),
  tar_target(
    meetnetdesign_raw,
    read_meetnetdesign_raw(meetnetdesign_file, cyclus_lookup)
  ),
  tar_target(popsizes, compute_popsizes(clustering)),
  tar_target(samplesizes, compute_samplesizes(meetnetdesign_raw)),
  tar_target(
    meetnetdesign,
    build_meetnetdesign(meetnetdesign_raw, popsizes, samplesizes)
  ),
  tar_target(
    inaccessible,
    inaccessible_ifbl
  ),
  # If the underlying data change, invalidate this target
  tar_target(florabank_data, fetch_florabank_data()),
  tar_target(hokken, fetch_florabank_hokken()),
  tar_target(
    kmhokken,
    hokken |>
      dplyr::filter(HoktypeID == 1) |>
      sf::st_as_sf(wkt = "GeomText", crs = 4326)
  ),
  tar_target(
    uurhokken,
    hokken |>
      dplyr::filter(HoktypeID != 1) |>
      sf::st_as_sf(wkt = "GeomText", crs = 4326)
  ),
  tar_target(
    vlaanderen_sf,
    read_sf(
      paste0(
        "OAPIF:", "https://geo.api.vlaanderen.be/VRBG2025/ogc/features/v1/"
      ),
      layer = "Refgew"
    )
  ),
  tar_target(
    fig_stratificatie,
    kmhokken |>
      inner_join(clustering, by = join_by(Code == kmhok)) |>
      ggplot() +
      geom_sf(data = uurhokken, fill = NA, colour = "lightgrey") +
      geom_sf(
        aes(fill = group, colour = group)
      ) +
      geom_sf(data = vlaanderen_sf, fill = NA) +
      theme_minimal()
  ),
  tar_target(
    fig_steekproef,
    kmhokken |>
      inner_join(meetnetdesign, by = join_by(Code == ifbl)) |>
      filter(cyclus == 1) |>
      ggplot() +
      geom_sf(data = uurhokken, fill = NA, colour = "lightgrey") +
      geom_sf(
        aes(fill = cyclus_jaar, colour = cyclus_jaar)
      )  +
      geom_sf(data = vlaanderen_sf, fill = NA) +
      facet_wrap(~ group) +
      theme_minimal()
  ),
  tar_target(
    svydata_ap,
    join_survey_data(meetnetdesign, florabank_data)
  ),
  tar_target(
    svydata_filtered,
    filter_survey_data(
      svydata = svydata_ap,
      ifbl_remove = inaccessible,
      min_n_taxa = richness_threshold
    )
  ),
  # ---- taxon clustering (subspecies/varieties -> main species) --------
  tar_target(soortenlijst, build_soortenlijst(svydata_filtered)),
  tar_target(names_clustered, cluster_taxon_names(soortenlijst)),
  tar_target(
    svydata_clustered,
    build_svydata_clustered(svydata_filtered, names_clustered)
  ),
  tar_target(svydata_clus_wide, build_svydata_wide(svydata_clustered)),
  tar_target(species, get_species_list(svydata_clustered)),
  # ---- survey designs ---------------------------------------------------
  tar_target(
    design_weight_panel,
    build_survey_design(svydata_clus_wide, "weight_panel")
  ),
  tar_target(
    design_weight_cyclus,
    build_survey_design(svydata_clus_wide, "weight_cyclus")
  ),
  # ---- status & change: overall -----------------------------------------
  tar_target(cyclus_results, get_all_status(design_weight_cyclus, species)),
  tar_target(change_results, get_all_change(design_weight_cyclus, species)),
  tar_target(
    absolute_wijziging_plot,
    plot_change(
      change_results
    )
  ),
  # ---- status & change: by stratum ---------------------------------------
  tar_target(
    cyclus_stratum_results,
    get_all_status_stratum(design_weight_cyclus, species)
  ),
  tar_target(
    change_stratum_results,
    get_all_change_stratum(design_weight_cyclus, species)
  ),
  tar_target(
    absolute_wijziging_stratum_plot,
    plot_change_stratum(
      change_stratum_results
    )
  ),
  # ---- relative change (ratio, species run in parallel via mirai) -------
  tar_target(
    rel_change_results_stratum,
    get_all_rel_change_stratum(design_weight_cyclus, species)
  ),
  tar_target(
    rel_change_results,
    get_all_rel_change(design_weight_cyclus, species)
  ),
  tar_target(
    relatieve_wijziging_stratum_plot,
    plot_rel_change_stratum(
      rel_change_results_stratum
    )
  ),
  tar_target(
    relatieve_wijziging_plot,
    plot_rel_change(
      rel_change_results
    )
  )
)
