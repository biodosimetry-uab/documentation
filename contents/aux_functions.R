theory_plot_estimated_dose_curve <- function(est_doses, fit_coeffs, fit_var_cov_mat,
                                             protracted_g_value, conf_int_curve,
                                             aberr_name, show_projections = FALSE) {
  # Validate est_doses names
  assessments <- names(est_doses)
  biodosetools:::match_names(assessments, c("whole", "partial", "hetero"))

  # Parse dose estimation list
  est_full_doses <- data.frame(
    dose = c(
      if ("whole" %in% assessments) est_doses$whole$est_doses[["dose"]],
      if ("partial" %in% assessments) est_doses$partial$est_doses[["dose"]],
      if ("hetero" %in% assessments) est_doses$hetero$est_doses[["dose1"]],
      if ("hetero" %in% assessments) est_doses$hetero$est_doses[["dose2"]]
    ),
    yield = c(
      if ("whole" %in% assessments) est_doses$whole$est_doses[["yield"]],
      if ("partial" %in% assessments) est_doses$partial$est_doses[["yield"]],
      if ("hetero" %in% assessments) est_doses$hetero$est_yields[["yield1"]],
      if ("hetero" %in% assessments) est_doses$hetero$est_yields[["yield2"]]
    ),
    type = c(
      if ("whole" %in% assessments) rep("Whole-body", 3),
      if ("partial" %in% assessments) rep("Partial-body", 3),
      if ("hetero" %in% assessments) rep("Heterogeneous 1", 3),
      if ("hetero" %in% assessments) rep("Heterogeneous 2", 3)
    ),
    level = c("Lower", "Estimate", "Upper")
  )

  # Parse confidence intervals
  conf_int_text_whole <- biodosetools:::parse_conf_int_text(est_doses$whole$conf_int)
  conf_int_text_partial <- biodosetools:::parse_conf_int_text(est_doses$partial$conf_int)
  conf_int_text_hetero <- biodosetools:::parse_conf_int_text(est_doses$hetero$conf_int)

  # Rightmost limit of the plot
  max_dose <- 1.05 * est_full_doses[["dose"]] %>%
    ifelse(is.na(.), 0, .) %>%
    max()

  # Generalised fit coefficients and variance-covariance matrix
  general_fit_coeffs <- biodosetools:::generalise_fit_coeffs(fit_coeffs[, "estimate"])
  general_fit_var_cov_mat <- biodosetools:::generalise_fit_var_cov_mat(fit_var_cov_mat)

  # Correct CIs
  # TODO: parse from est_doses
  conf_int_curve <- conf_int_curve %>%
    biodosetools:::correct_conf_int(general_fit_var_cov_mat, protracted_g_value, type = "curve")

  # Plot data from curves
  curves_data <- data.frame(
    dose = seq(0, max_dose, length.out = 100)
  ) %>%
    dplyr::mutate(
      yield = biodosetools:::calculate_yield(.data$dose, type = "estimate", general_fit_coeffs, NULL, protracted_g_value, 0),
      yield_low = biodosetools:::calculate_yield(.data$dose, type = "lower", general_fit_coeffs, general_fit_var_cov_mat, protracted_g_value, conf_int_curve),
      yield_upp = biodosetools:::calculate_yield(.data$dose, type = "upper", general_fit_coeffs, general_fit_var_cov_mat, protracted_g_value, conf_int_curve)
    )

  # Parse assessment legend
  color_breaks <- c("Whole-body", "Partial-body", "Heterogeneous 1", "Heterogeneous 2")
  color_labels <- c(
    paste("Whole-body", conf_int_text_whole),
    paste("Partial-body", conf_int_text_partial),
    paste("Heterogeneous 1", conf_int_text_hetero),
    paste("Heterogeneous 2", conf_int_text_hetero)
  )
  color_indices <- est_full_doses$type %>%
    unique() %>%
    paste(collapse = "|") %>%
    grep(color_breaks)

  # Make base plot
  gg_curve <- ggplot2::ggplot(curves_data) +
    # Fitted curve
    ggplot2::stat_function(
      data = data.frame(x = c(0, max_dose)),
      mapping = ggplot2::aes(x = .data$x),
      fun = function(x) biodosetools:::yield_fun(x, general_fit_coeffs, protracted_g_value),
      linetype = "dashed"
    ) +
    ggplot2::labs(x = "Dose (Gy)", y = paste0(aberr_name, "/cells")) +
    ggplot2::theme_bw()

  if (conf_int_curve > 0) {
    # Confidence bands (Merkle, 1983)
    gg_curve <- gg_curve +
      ggplot2::geom_ribbon(
        data = curves_data,
        mapping = ggplot2::aes(x = .data$dose, ymin = .data$yield_low, ymax = .data$yield_upp),
        color = "black",
        linetype = "dashed",
        alpha = 0
      )
  }

  # Add doses to plot
  gg_curve <- gg_curve +
    # Estimated whole-body doses
    ggplot2::geom_point(
      data = est_full_doses,
      mapping = ggplot2::aes(x = .data$dose, y = .data$yield, shape = .data$level),
      size = 2, na.rm = TRUE
    ) +
    # Assessment
    ggplot2::scale_color_manual(
      values = grDevices::hcl(
        h = seq(15, 375, length = 4 + 1),
        l = 65,
        c = 100
      ) %>%
        .[1:4] %>%
        `names<-`(c("Partial-body", "Heterogeneous 1", "Heterogeneous 2", "Whole-body")),
      breaks = color_breaks[color_indices],
      labels = color_labels[color_indices]
    ) +
    # Estimation level
    ggplot2::scale_shape_manual(
      values = c("Lower" = 15, "Estimate" = 16, "Upper" = 17),
      breaks = c("Lower", "Estimate", "Upper")
    ) +
    ggplot2::guides(
      color = "none",
      shape = ggplot2::guide_legend(order = 2),
      fill = ggplot2::guide_legend(order = 3)
    ) +
    ggplot2::labs(color = "Assessment", shape = "Estimation") +
    # Tweak legend
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8),
      legend.spacing.y = ggplot2::unit(5, "points"),
      legend.key.height = ggplot2::unit(12, "points")
    )

  # Show dose projections
  if (show_projections) {
    gg_curve <- gg_curve +
      ggplot2::geom_segment(
        data = est_full_doses,
        mapping = ggplot2::aes(x = 0, xend = dose, y = yield, yend = yield),
        linetype = "dotted"
      ) +
      ggplot2::geom_segment(
        data = est_full_doses,
        mapping = ggplot2::aes(x = dose, xend = dose, y = 0, yend = yield),
        linetype = "dotted"
      )
  }

  gg_curve <- gg_curve +
    ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::number_format(accuracy = 0.1)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))

  return(gg_curve)
}


fix_estimated_dose_curve_color_labels <- function(gg, breaks, labels) {
  gg <- gg +
    ggplot2::scale_color_manual(
      values = grDevices::hcl(
        h = seq(15, 375, length = 4 + 1),
        l = 65, c = 100
      ) %>%
        .[1:4] %>%
        `names<-`(c("Partial-body", "Heterogeneous 1", "Heterogeneous 2", "Whole-body")),
      breaks = breaks,
      labels = labels
    )

  return(gg)
}


plot_fit_dose_curve_alt <- function(fit_results_list, aberr_name, max_dose, show_count_data = TRUE, conf_int = 0.95) {
  # Read objects from fit results list
  count_data <- as.data.frame(fit_results_list[["fit_raw_data"]])
  fit_coeffs <- fit_results_list[["fit_coeffs"]]
  fit_var_cov_mat <- fit_results_list[["fit_var_cov_mat"]]

  # Generalised fit coefficients
  general_fit_coeffs <- biodosetools:::generalise_fit_coeffs(fit_coeffs[, "estimate"])

  # Generalised variance-covariance matrix
  general_fit_var_cov_mat <- biodosetools:::generalise_fit_var_cov_mat(fit_var_cov_mat)

  # Generalised curves
  yield_fun <- function(d) {
    general_fit_coeffs[["coeff_C"]] +
      general_fit_coeffs[["coeff_alpha"]] * d +
      general_fit_coeffs[["coeff_beta"]] * d^2
  }

  chisq_df <- nrow(fit_coeffs)
  R_factor <- sqrt(stats::qchisq(conf_int, df = chisq_df))

  yield_error_fun <- function(d) {
    sqrt(
      general_fit_var_cov_mat[["coeff_C", "coeff_C"]] +
        general_fit_var_cov_mat[["coeff_alpha", "coeff_alpha"]] * d^2 +
        general_fit_var_cov_mat[["coeff_beta", "coeff_beta"]] * d^4 +
        2 * general_fit_var_cov_mat[["coeff_C", "coeff_alpha"]] * d +
        2 * general_fit_var_cov_mat[["coeff_C", "coeff_beta"]] * d^2 +
        2 * general_fit_var_cov_mat[["coeff_alpha", "coeff_beta"]] * d^3
    )
  }

  # Plot data
  plot_data <- count_data %>%
    dplyr::mutate(
      yield = X / N,
      dose = D
    ) %>%
    dplyr::select(dose, yield)

  curves_data <- data.frame(dose = seq(0, max_dose, length.out = 100)) %>%
    dplyr::mutate(
      yield = yield_fun(dose),
      yield_low = yield_fun(dose) - R_factor * yield_error_fun(dose),
      yield_upp = yield_fun(dose) + R_factor * yield_error_fun(dose)
    )

  # Make plot
  gg_curve <- ggplot2::ggplot(plot_data)

  if (show_count_data) {
    # Observed data
    gg_curve <- gg_curve +
      ggplot2::geom_point(
        mapping = ggplot2::aes(x = dose, y = yield)
      )
  }
  # Fitted curve
  gg_curve <- gg_curve +
    ggplot2::stat_function(
      data = data.frame(x = c(0, max(plot_data[["dose"]]))),
      mapping = ggplot2::aes(x),
      fun = function(x) yield_fun(x),
      linetype = "dashed"
    ) +
    # Confidence bands (Merkle, 1983)
    ggplot2::geom_ribbon(
      data = curves_data,
      ggplot2::aes(x = dose, ymin = yield_low, ymax = yield_upp),
      alpha = 0.25
    ) +
    ggplot2::labs(
      x = "Dose (Gy)",
      y = paste0(aberr_name, "/cells")
    ) +
    ggplot2::theme_bw()

  # Return object
  return(gg_curve)
}
