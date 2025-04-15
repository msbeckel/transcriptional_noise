#' @title GLMM Transcriptional Noise Analysis
#' @description Performs generalized linear mixed modeling of transcriptional noise
#' @param seurat Seurat object with metadata containing:
#' - transcriptional_noise: Response variable (positive continuous)
#' - Age: Fixed effect (numeric)
#' - hash.ID: Random effect (subject ID)
#' - Celltype.LowRes: Cell type annotations
#' @param family GLM family (default: Gamma(link = "log"))
#' @return List containing:
#' - model_results: Dataframe with model coefficients and statistics
#' - diagnostics: Residual diagnostic plots
#' - comparisons: Model fit comparisons (AIC/BIC)
#' @details
#' - Uses GLMM to handle non-normal distributions without transformations
#' - Implements comprehensive model diagnostics
#' - Compares multiple family/link function combinations
#' - Creates output directory with results
#' - Generates interactive plots for better visualization
#' 
#.libPaths("/ngs/mllorens/mbeckel/R_libs")

run_glmm_noise_analysis <- function(seurat,  
                                    response_var,
                                    output_dir = "",
                                    fixed_effect = "Age",
                                    random_effect = "hash.ID",
                                    cell_type_col = "Celltype.LowRes",
                                    celltype_order = c("Microglia", "Oligodendro", 
                                                       "Astrocyte_qNSC", "aNSC_NPC", 
                                                       "Neuroblast"),
                                    family = Gamma(link = "log")) {
  
  # 1. Input Validation
  required_cols <- c(response_var, fixed_effect, random_effect, cell_type_col)
  if(!all(required_cols %in% colnames(seurat@meta.data))) {
    stop("Missing required metadata columns")
  }
  
  # 2. Create Output Directory
  output_dir <- paste0(output_dir, "GLMM_Results_", response_var,"_", family$family, "_", family$link)
  dir.create(output_dir, showWarnings = FALSE)
  
  # 3. Data Preparation
  metadata <- seurat@meta.data %>%
    filter(.data[[response_var]] > 0) %>%
    mutate(Celltype.LowRes = factor(.data[[cell_type_col]], levels = celltype_order))
  
  # Scale response and fixed effect variables
  metadata[paste0(c(response_var, fixed_effect), "_scaled")] <- scale(metadata[c(response_var, fixed_effect)], scale = TRUE, center = FALSE)
  
  # Update variable names
  response_var <- paste0(response_var, "_scaled")
  fixed_effect <- paste0(fixed_effect, "_scaled")
  
  # 4. Define GLMM Function
  fit_glmm <- function(data, formula, family) {
    glmer(
      formula = formula,
      data = data,
      family = family,
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
  }
  
  formula <- paste0(response_var," ~ ", fixed_effect, "+ (1 |", random_effect, ")")
  # 5. Fit Models by Cell Type
  model_results <- metadata %>%
    group_by(Celltype.LowRes) %>%
    group_map(~{
      if(n_distinct(.x$hash.ID) < 3) return(NULL)
      
      tryCatch({
        model <- fit_glmm(.x, formula, family)
        
        # Calculate pseudo-R2
        r2 <- MuMIn::r.squaredGLMM(model)
        
        # DHARMa diagnostics
        simulation <- DHARMa::simulateResiduals(model)
        
        list(
          model = model,
          summary = broom.mixed::tidy(model, conf.int = TRUE),
          r2 = r2,
          diagnostics = simulation
        )
      }, error = function(e) print(e))
    }) %>%
    setNames(levels(metadata$Celltype.LowRes))
  
  # 6. Process Results
  results_df <- map_dfr(model_results, ~{
    if(is.null(.x)) return(NULL)
    .x$summary %>%
      filter(effect == "fixed") %>%
      mutate(
        Significance = case_when(
          p.value < 0.001 ~ "***",
          p.value < 0.01 ~ "**",
          p.value < 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        R2_conditional = .x$r2[1],
        R2_marginal = .x$r2[2]
      )
  }, .id = "Celltype")
  
  # 7. Generate Diagnostics
  #diagnostic_plots <- map(model_results, ~{
  #  if(is.null(.x)) return(NULL)
  #  DHARMa::plotResiduals(.x$diagnostics)
  #})
  
  # 8. Forest Plot
  forest_plot <- results_df %>%
    filter(term == "Age_scaled") %>%
    mutate(cell_type = factor(Celltype, levels = celltype_order)) %>%
    ggplot(aes(x = cell_type, y = estimate)) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
    scale_x_discrete(limits = rev) +
    labs(title = 'Age Effect Sizes with 95% CI') +
    theme_bw()
  
  # 9. Scatter Plot
  scatter_data <- metadata %>%
    mutate(cell_type = factor(.data[[cell_type_col]], levels = celltype_order))
  
  scatter_plot <- ggplot(scatter_data, aes(x = .data[[fixed_effect]], y = .data[[response_var]])) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x) +
    facet_wrap(~cell_type, scales = "free_y") +
    labs(x = fixed_effect, y = response_var) +
    theme_bw()
  
  # Interactive Versions
  forest_interactive <- ggplotly(forest_plot)
  scatter_interactive <- ggplotly(scatter_plot)
  
  plots <- list(
    forest = forest_plot,
    scatter = scatter_plot,
    forest_interactive = forest_interactive,
    scatter_interactive = scatter_interactive
  )
  
  # 10. Save Outputs
  save_glmm_outputs(output_dir, results_df, model_results, diagnostic_plots, plots)
  
  # 11. Return Object
  list(
    results = results_df,
    models = map(model_results, "model")
  )
}

# Helper Functions --------------------------------------------------------
save_glmm_outputs <- function(dir, results, models, diag_plots, plots) {
  # Save results
  write.csv(results, file.path(dir, "glmm_results.csv"), row.names = FALSE)
  
  # Save model objects
  saveRDS(models, file.path(dir, "glmm_models.rds"))
  
  # Save diagnostic plots
  #pdf(file.path(dir, "diagnostic_plots.pdf"))
  #walk(diag_plots, print)
  #dev.off()
  
  # Save plots
  ggsave(file.path(dir, "forest_plot.png"), plots$forest, width = 10, height = 8)
  ggsave(file.path(dir, "scatter_plot.png"), plots$scatter, width = 12, height = 10)
  
  # Save interactive plots
  htmlwidgets::saveWidget(plots$forest_interactive, 
                          file.path(dir, "interactive_forest.html"))
  htmlwidgets::saveWidget(plots$scatter_interactive, 
                          file.path(dir, "interactive_scatter.html"))
  
  # Create summary report
  render_glmm_report(dir, results)
}

render_glmm_report <- function(output_dir, results_df, diagnostic_plots, model_objects) {
  # Create temporary Rmd file
  rmd_content <- paste0(
    "---\n",
    "title: 'GLMM Transcriptional Noise Analysis Report'\n",
    "output: html_document\n",
    "---\n\n",
    
    "```{r setup, include=FALSE}\n",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)\n",
    "library(ggplot2)\n",
    "library(dplyr)\n",
    "```\n\n",
    
    "# GLMM Analysis Results\n\n",
    
    "## Model Summary\n",
    "```{r summary}\n",
    "results <- read.csv('", file.path(output_dir, "glmm_results.csv"), "') %>% \n",
    "  filter(term == 'Age_scaled') %>% \n",
    "  select(-c(term, effect, group))\n",
    "knitr::kable(results, caption = 'Model Coefficients and Significance')\n",
    "```\n\n",
    
    "## Forest plot\n",
    "```{r forest, echo=FALSE, out.width = '100%'}\n",
    "knitr::include_graphics('", file.path(output_dir, "forest_plot.png"), "')\n",
    "```\n\n",
    
    "## Scatter plot\n",
    "```{r scatter, echo=FALSE, out.width = '100%'}\n",
    "knitr::include_graphics('", file.path(output_dir, "scatter_plot.png"), "')\n",
    "```\n\n",
    
    "## Model Diagnostics\n",
    "```{r diagnostics}\n",
    "models <- readRDS('", file.path(output_dir, "glmm_models.rds"), "')\n",
    "celltypes <- names(models)[!sapply(models, is.null)]\n",
    "```\n\n",
    
    "### Residual Diagnostics\n",
    "```{r residual-plots}\n",
    "for(ct in celltypes) {\n",
    "  print(ct)\n",
    "  dharma_res <- DHARMa::simulateResiduals(models[[ct]]$model)\n",
    "  print(DHARMa::plotResiduals(dharma_res))\n",
    "}\n",
    "```\n\n",
    
    "### Model Fit Statistics\n",
    "```{r fit-stats}\n",
    "fit_stats <- data.frame(\n",
    "  Celltype = character(),\n",
    "  AIC = numeric(),\n",
    "  BIC = numeric(),\n",
    "  stringsAsFactors = FALSE\n",
    ")\n",
    "\n",
    "for(ct in celltypes) {\n",
    "  fit_stats <- rbind(fit_stats, data.frame(\n",
    "    Celltype = ct,\n",
    "    AIC = AIC(models[[ct]]$model),\n",
    "    BIC = BIC(models[[ct]]$model)\n",
    "  ))\n",
    "}\n",
    "\n",
    "knitr::kable(fit_stats, caption = 'Model Fit Statistics')\n",
    "```\n"
  )
  
  rmd_file <- file.path(output_dir, "report.Rmd")
  writeLines(rmd_content, rmd_file)
  
  # Render the report
  rmarkdown::render(
    input = rmd_file,
    output_format = "html_document",
    output_dir = output_dir,
    quiet = TRUE
  )
  
  # Cleanup temporary Rmd
  file.remove(rmd_file)
}

