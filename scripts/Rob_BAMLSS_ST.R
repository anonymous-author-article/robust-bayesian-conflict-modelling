# ==============================================================================
# ROBUST BAYESIAN SPATIO-TEMPORAL MODELLING OF CONFLICT DATA
# ==============================================================================

rm(list = ls()); gc()

# ==============================================================================
# 1. ENVIRONMENT SETUP AND PACKAGE INITIALIZATION
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  bamlss, gamlss.dist, dplyr, tidyr, ggplot2, ggrepel, ggspatial,
  lubridate, sf, spdep, scales, gridExtra, patchwork, ncvreg, mgcv,
  parallel, RColorBrewer, kableExtra, viridis, scoringRules, statmod,
  boot, coda,stringr,cowplot,geodata,here
)


theme_set(theme_bw(base_size = 12) + 
            theme(panel.grid.minor = element_blank(),
                  text = element_text(family = "serif")))

set.seed(1234) 

# ==============================================================================
# 2. DATA PRE-PROCESSING AND SPATIAL WEIGHT MATRIX CONSTRUCTION
# ==============================================================================

# The 'here' function automatically finds the root of the project.
df_raw <- readRDS(here("data", "FinalDATASET.rds"))
shapefile <- readRDS(here("data", "cabo_delgado_boundary.rds"))

if(inherits(df_raw, "sf")) df_raw <- st_drop_geometry(df_raw)

# Ensure consistency between administrative boundaries and data records
common_districts <- intersect(unique(df_raw$admin2), shapefile$NAME_2)
if(length(common_districts) == 0) stop("ERRO: District names do not match.")

shape_subset <- shapefile[shapefile$NAME_2 %in% common_districts, ]
shape_subset <- shape_subset[order(shape_subset$NAME_2), ] 

df_clean <- df_raw %>% 
  filter(admin2 %in% common_districts) %>%
  distinct(event_id, date, admin2, fatalities, .keep_all = TRUE)

# Constructing the Spatial Adjacency Matrix (Queen's contiguity) and Penalty Matrix
nb_map <- poly2nb(shape_subset, row.names = shape_subset$NAME_2)
adj_mat <- nb2mat(nb_map, style = "B", zero.policy = TRUE)
penalty_mat <- diag(rowSums(adj_mat)) - adj_mat
rownames(penalty_mat) <- colnames(penalty_mat) <- rownames(adj_mat)

# Extracting centroids for spatial interaction terms
centroids_sf <- st_centroid(shape_subset)
coords_mat <- st_coordinates(centroids_sf)
df_coords <- data.frame(NAME_2 = shape_subset$NAME_2, lon = coords_mat[,1], lat = coords_mat[,2])

vars_to_scale <- c("dist_montepuez_ruby", "dist_balama_graphite", 
                   "dist_mocimboa_port", "dist_afungi_lng",
                   "dist_mueda_Airport", "dist_pemba_Airport",
                   "SEXRATIO", "DEPND", "precip_anomaly")

#Time indexing, factor levelling, and variable scaling (Z-scores)
df_model <- df_clean %>%
  left_join(df_coords, by = c("admin2" = "NAME_2")) %>%
  mutate(
    date_obj = as.Date(paste0(year_month, "-01")),
    financing_phase = case_when(
      date_obj < as.Date("2017-11-01") ~ "1-Local",
      date_obj >= as.Date("2017-11-01") & date_obj < as.Date("2019-04-01") ~ "2-Transition",
      date_obj >= as.Date("2019-04-01") ~ "3-Regional"
    ) %>% factor(levels = c("1-Local", "2-Transition", "3-Regional")),
    t_index = as.numeric(factor(year_month, levels = unique(year_month[order(as.Date(paste0(year_month, "-01")))]))),
    admin2 = factor(admin2, levels = rownames(penalty_mat)),
    log_pop = log(POP + 1),
    classif_actor1 = tidyr::replace_na(as.character(classif_actor1), "Unidentified"),
    across(all_of(vars_to_scale), ~ as.numeric(scale(.)), .names = "{.col}_s")
  ) %>%
  filter(classif_actor1 != "Civilian groups") %>%
  droplevels()

message("Data prep complete. N: ", nrow(df_model))



# ==============================================================================
# 3. EXPLORATORY DATA ANALYSIS: RESPONSE VARIABLE DISTRIBUTION
#    Assessing the marginal distribution of fatality counts to inform family selection.
# ==============================================================================
par(mfrow=c(2,3))

hist(df_model$fatalities, main = "Raw fatality counts", xlab = "Fatalities")
hist(df_model$fatalities[df_model$fatalities>0], main = "Only fatalities > 0", xlab = "Fatalities")
hist(df_model$fatalities[df_model$fatalities>1], main = "Only fatalities > 1", xlab = "Fatalities")
hist(df_model$fatalities[df_model$fatalities>1 & df_model$fatalities<=60], main = "1 <Fatalities <50", xlab = "Fatalities")
hist(log(df_model$fatalities +1), main = "log transformation", xlab = "log(Fatality + 1)")
hist(log(df_model$fatalities/exp(df_model$log_pop)), main = "log rate transformation", xlab = "log(Fatality/Population)")


# ==============================================================================
# 4. DESCRIPTIVE TEMPORAL ANALYSIS AND PEAK EVENT IDENTIFICATION
#    Visualizing temporal trends, seasonality, and identifying high-fatality outlier events
#    to produce Figure 1 (Descriptive Analysis).
# ==============================================================================

theme_set(theme_bw(base_size = 12) + 
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(color = "gray40", size = 11)))

df_year <- df_model %>%
  group_by(year(date)) %>%
  summarise(total_fatalities = sum(fatalities, na.rm = TRUE)) %>%
  ungroup()

df_month_season <- df_model %>%
  mutate(month_lbl = month(date_obj, label = TRUE, abbr = TRUE)) %>% 
  group_by(month_lbl) %>%
  summarise(total_fatalities = sum(fatalities, na.rm = TRUE)) %>%
  ungroup()

df_timeline <- df_model %>%
  group_by(date_obj) %>%
  summarise(total_fatalities = sum(fatalities, na.rm = TRUE)) %>%
  ungroup()

# Thresholding for annotating peak violence months
threshold <- mean(df_timeline$total_fatalities) + 1.5 * sd(df_timeline$total_fatalities)

peak_months <- df_timeline %>%
  filter(total_fatalities > threshold) %>%
  arrange(desc(total_fatalities)) %>%
  head(7) 

get_headline_desc <- function(target_date, raw_data) {
  start_m <- floor_date(target_date, "month")
  end_m   <- ceiling_date(target_date, "month") - days(1)
  
  top_event <- raw_data %>%
    filter(date_obj >= start_m & date_obj <= end_m) %>%
    arrange(desc(fatalities)) %>%
    slice(1)
  
  if(nrow(top_event) == 0) return("Violence Peak")
  
  type <- if("sub_event_type" %in% names(top_event)) top_event$sub_event_type else top_event$event_type
  loc <- as.character(top_event$admin2)
  deaths <- top_event$fatalities
  
  return(paste0(type, " - ", loc, " (", deaths, "†)"))
}

annotations <- peak_months %>%
  rowwise() %>%
  mutate(legend_label = get_headline_desc(date_obj, df_model)) %>%
  ungroup() %>%
  mutate(plot_label = gsub(" - ", "\n", legend_label))


g_year <- ggplot(df_year, aes(x = factor(`year(date)`), y = total_fatalities)) +
  geom_col(fill = "#2c3e50", width = 0.7, alpha = 0.9) +
  geom_text(aes(label = comma(total_fatalities)), vjust = -0.5, size = 3.5, fontface="bold") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = NULL, y = "Fatalities") +
  theme(axis.text.x = element_text(face="bold"))+
  theme_minimal()

g_season <- ggplot(df_month_season, aes(x = month_lbl, y = total_fatalities)) +
  geom_col(fill = "#e67e22", width = 0.7, alpha = 0.8) +
  geom_text(aes(label = comma(total_fatalities)), vjust = -0.5, size = 3) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs( x = NULL, y = NULL)

unique_events <- unique(annotations$legend_label)
n_events <- length(unique_events)
palette_pts <- brewer.pal(n = max(3, n_events), name = "Dark2")[1:n_events]
names(palette_pts) <- unique_events

phase_dates <- as.Date(c("2017-11-01", "2019-04-01"))
min_date <- min(df_timeline$date_obj)
max_date <- max(df_timeline$date_obj)

df_phases_rect <- data.frame(
  xmin = c(min_date, phase_dates[1], phase_dates[2]),
  xmax = c(phase_dates[1], phase_dates[2], max_date),
  Phase = factor(c("1. Local Funding", "2. Transition", "3. Regional"),
                 levels = c("1. Local Funding", "2. Transition", "3. Regional")),
  y_text = max(df_timeline$total_fatalities) * 0.95
)%>%
mutate(Phase= case_when(
  Phase=="1. Local Funding" ~ "1",
  Phase=="2. Transition" ~ "2",
  Phase=="3. Regional" ~ "3"
  
))

cols_combined <- c(palette_pts, "Trend (Loess)" = "#c0392b") 

fills_combined <- c(
  "1" = "#e8f6f3",
  "2" = "#fcf3cf",
  "3" = "#fadbd8",
  "Observed Fatalities" = "#34495e"
)

# Annotation function for plot labels
get_headline_desc <- function(target_date, raw_data) {
  start_m <- floor_date(target_date, "month")
  end_m   <- ceiling_date(target_date, "month") - days(1)
  month_data <- raw_data %>%
    filter(date_obj >= start_m & date_obj <= end_m)
  total_month <- sum(month_data$fatalities, na.rm = TRUE)
  top_event <- month_data %>%
    arrange(desc(fatalities)) %>%
    slice(1)
  if(nrow(top_event) == 0) return("Violence Peak")
  type <- if("sub_event_type" %in% names(top_event)) top_event$sub_event_type else top_event$event_type
  loc <- as.character(top_event$admin2)
  deaths <- top_event$fatalities
  return(paste0(type, " - ", loc, "\nEvent: ", deaths, " | Month Total: ", total_month))
}

annotations <- peak_months %>%
  rowwise() %>%
  mutate(legend_label = get_headline_desc(date_obj, df_model)) %>%
  ungroup() %>%
  mutate(plot_label = legend_label)

g_timeline <- ggplot() + 
  geom_vline(xintercept = phase_dates, linetype = "dashed", color = "gray60", size = 0.4) +
  geom_col(data = df_timeline, aes(x = date_obj, y = total_fatalities), 
           width = 20, alpha = 0.9) +
  geom_smooth(data = df_timeline, aes(x = date_obj, y = total_fatalities, color = "Trend (Loess)"), 
              method = "loess", span = 0.2, 
              fill = "gray70", alpha = 0.3, size = 0.8, se = TRUE) +
  geom_point(data = annotations, 
             aes(x = date_obj, y = total_fatalities, color = legend_label), 
             size = 3, stroke = 1, shape = 21, fill = "white",
             show.legend = FALSE) +
  geom_label_repel(data = annotations, aes(x = date_obj, y = total_fatalities, label = plot_label),
                   size = 3,  box.padding = 0.6, 
                   nudge_y = max(df_timeline$total_fatalities)*0.15, 
                   segment.color = "gray40", segment.size = 0.3, min.segment.length = 0) +
  scale_x_date(limits = c(min(df_model$date_obj), max(df_model$date_obj)), 
               date_breaks = "1 year", date_labels = "%Y", expand = c(0.01,0)) +
  scale_y_continuous(labels = comma)  +
  scale_color_manual(name = NULL, 
                     values = c(cols_combined, "Trend (Loess)" = "#c0392b"),
                     breaks = c("Trend (Loess)")) + 
  labs(subtitle = "Observed fatalities and peak events",
    y = "Monthly Fatalities (Total)",
    x = NULL) + 
  theme(
    legend.position = c(0.85, 0.9),
    legend.box = "horizontal",
    legend.background = element_rect(fill = "white", color = "gray90", size = 0.2),
    legend.key = element_blank(),
    legend.text = element_text(size = 8),
    legend.spacing.x = unit(0.2, 'cm'),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 0, r = 10, l = 10) 
  )

g_phases <- ggplot(df_phases_rect) +
  geom_segment(aes(x = xmin, xend = xmax, y = 0.8, yend = 0.8), color = "black", size = 0.3) +
  geom_segment(aes(x = xmin, xend = xmin, y = 0.6, yend = 1.0), color = "black", size = 0.3) +
  geom_segment(aes(x = xmax, xend = xmax, y = 0.6, yend = 1.0), color = "black", size = 0.3) +
  geom_text(aes(x = xmin + (xmax - xmin)/2, y = 0.35, label = Phase), size = 3, color = "black") +
  
  scale_x_date(limits = c(min(df_model$date_obj), max(df_model$date_obj)), 
               date_breaks = "1 year", date_labels = "%Y", expand = c(0.01,0)) +
  scale_y_continuous(limits = c(0, 1.05)) + 
  
  labs(x = NULL, y = NULL, 
       caption = "Note: Financing phases (1=local, 2= transition and 3=regional) derived from Weiss et al. (2023).\nBars represent monthly aggregates; labels detail the deadliest event.") +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 9, color = "black", margin = margin(t=2)), # Margem do ano reduzida
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.length.x = unit(0.1, "cm"),
    plot.margin = margin(t = -5, r = 10, l = 10, b = 10),
    plot.caption = element_text(hjust = 0, color = "gray40", size = 9, margin = margin(t=10))
  )


final_plot <- (g_year + g_season) / 
  g_timeline / 
  g_phases + 
  plot_layout(heights = c(1, 1.5, 0.2)) + 
  plot_annotation(tag_levels = list(c('A', 'B', 'C', ''))) 

print(final_plot)

g_phases <- g_phases + scale_y_continuous(limits = c(0, 1.25))

a <- g_timeline / g_phases +
  plot_layout(heights = c(10, 1)) 
ggsave(here("figures","Figure_Descriptive_Final_Robust.pdf"), a, width = 12, height = 10)

# ==============================================================================
# 5. SPATIAL DISTRIBUTION OF OBSERVED AGGREGATE FATALITIES
# ==============================================================================
df_risk_total <- df_model %>%
  group_by(admin2) %>%
  summarise(total_risk = sum(fatalities, na.rm = TRUE)) %>%
  ungroup()

map_data_cd <- shape_subset %>%
  left_join(df_risk_total, by = c("NAME_2" = "admin2")) %>%
  mutate(risk_cat = cut(total_risk, 
                        breaks = c(0, 50, 100, 500, 1000, 1500, Inf),
                        labels = c("<50", "51-100", "101-500", "501-1000", "1001-1500", ">1500"),
                        include.lowest = TRUE))

moz_sf <- geodata::gadm(country = "MOZ", level = 1, path = tempdir()) %>% st_as_sf()
moz_sf$is_cd <- ifelse(moz_sf$NAME_1 == "Cabo Delgado", "Yes", "No")

p_main <- ggplot(map_data_cd) +
  geom_sf(data = shapefile, fill = "white", color = "black",size = 0.5) +
  geom_sf(aes(fill = risk_cat), color = "black", size = 0.5) +
  scale_fill_manual(
    values = c(
      "<50"       = "#ffffcc", 
      "51-100"    = "#c2e699", 
      "101-500"   = "#78c679", 
      "501-1000"  = "#31a354", 
      "1001-1500" = "#006837", 
      ">1500"     = "#003311"
    ),
    na.value = "white",
    name = ""
  ) +
  
  geom_sf_text(aes(label = NAME_2), size = 3, color = "black", fontface = "bold", check_overlap = TRUE) +
  annotation_scale(location = "br", width_hint = 0.3, pad_x = unit(0.2, "in")) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  
  theme_void() +
  theme(
    #legend.position = c(0.8, 0.2), 
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = NA, alpha = 0.8),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b=10)),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10)
  ) +
  labs(title = "", subtitle = "")

p_inset <- ggplot(moz_sf) +
  geom_sf(fill = "white", color = "gray60", size = 0.1) +
  geom_sf(data = subset(moz_sf, is_cd == "Yes"), fill = "black", color = NA) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

p_left_panel <- ggdraw(p_main) +
  draw_plot(p_inset, x = 0.01, y = 0.68, width = 0.25, height = 0.25)

ggsave(here("figures","Figure_map.pdf"), p_left_panel, width = 12, height = 10)


# ==============================================================================
# 6. Combining spatial and temporal descriptive plots
# ==============================================================================

p_right_panel <- g_timeline / g_phases + plot_layout(heights = c(4, 1))
ggsave(here("figures","Figure_tem.pdf"), p_right_panel, width = 12, height = 10)


p_left_panel_clean <- p_left_panel + theme(plot.margin = margin(0, 0, 0, 0)) 
p_right_panel_clean <- p_right_panel + theme(plot.margin = margin(0, 0, 0, 2)) 

final_figure <- plot_grid(
  p_left_panel_clean, 
  p_right_panel_clean, 
  ncol = 2, 
  rel_widths = c(1, 1.7), 
  labels = c("A", "B"),
  label_size = 14,
  align = "h",    
  axis = "tb"    
)

print(final_figure)

ggsave(here("figures","Figure_Descriptive.pdf"), final_figure, width = 12, height = 8)

# ==============================================================================
# 7. VARIABLE SELECTION AND OUTLIER DETECTION (PLE METHOD)
# ==============================================================================

X_mat <- model.matrix(~ type_of_Terrorism + classif_actor1 + log_pop + precip_anomaly+bord+
                        SEXRATIO_s + DEPND_s + dist_mueda_Airport_s + financing_phase+
                        dist_afungi_lng_s + cyclone_exposure_6m - 1, data = df_model)
Y_vec <- df_model$fatalities

fit_scad <- ncvreg(X_mat, Y_vec, family = "poisson", penalty = "SCAD")
best_lam <- fit_scad$lambda.min
coef_scad <- coef(fit_scad, lambda = best_lam)
selected_vars <- names(coef_scad)[coef_scad != 0]

preds_ple <- predict(fit_scad, X_mat, lambda = best_lam, type = "response")
res_ple <- (Y_vec - preds_ple) / sqrt(preds_ple + 1e-5)
outliers_ple_idx <- which(abs(res_ple) > 3.5)
outliers_table <- df_model[outliers_ple_idx, ] %>% select(date_obj, admin2, fatalities)

message(paste("PLE Detected", length(outliers_ple_idx), "potential outliers."))

# ==============================================================================
# 8. DEFINITION OF CUSTOM STATISTICAL FAMILIES
#    Implementing custom Zero-One Inflated Negative Binomial (ZOINB) and 
#    Rescaled Beta (RSB) mixture families, including derivatives for MCMC sampling.
# ==============================================================================
safe_log <- function(x, minval = 1e-300) log(pmax(x, minval))
clamp_vec <- function(x, lo = -10, hi = 10) pmax(pmin(x, hi), lo)

num_score_hess_factory <- function(d_fun) {
  score_fun <- function(y, par, id = NULL, ...) {
    delta <- 1e-5; par_p <- par; par_m <- par
    par_p[[id]] <- par[[id]] + delta; par_m[[id]] <- par[[id]] - delta
    lp <- safe_log(d_fun(y, par_p)); lm <- safe_log(d_fun(y, par_m))
    return((lp - lm)/(2*delta))
  }
  hess_fun <- function(y, par, id = NULL, ...) {
    s_p <- score_fun(y, lapply(par, function(x) x + 1e-5), id=id)
    s_m <- score_fun(y, lapply(par, function(x) x - 1e-5), id=id)
    return((s_p - s_m)/(2*1e-5) * -1)
  }
  list(score = score_fun, hess = hess_fun)
}

# ZOINB (Zero-One Inflated Negative Binomial) Implementation

f_ZOINB_fam <- function(...) {
  d_fun <- function(y, par, log = FALSE) {
    mu <- pmax(par$mu, 1e-8); sigma <- pmax(par$sigma, 1e-8)
    nu <- plogis(par$nu); tau <- plogis(par$tau)
    size <- 1/sigma
    
    p0 <- dnbinom(0, mu=mu, size=size)
    p1 <- dnbinom(1, mu=mu, size=size)
    denom <- pmax(1 - p0 - p1, 1e-10)
    
    d <- numeric(length(y))
    d[y==0] <- nu[y==0]
    d[y==1] <- (1 - nu[y==1]) * tau[y==1]
    
    idx_gt1 <- y > 1
    if(any(idx_gt1)) {
      d_nb <- dnbinom(y[idx_gt1], mu=mu[idx_gt1], size=size[idx_gt1])
      d[idx_gt1] <- (1 - nu[idx_gt1]) * (1 - tau[idx_gt1]) * (d_nb / denom[idx_gt1])
    }
    
    d <- pmax(d, 1e-300)
    if(log) return(log(d)) else return(d)
  }
  
  p_fun <- function(y, par, ...) {
    mu <- pmax(par$mu, 1e-8); sigma <- pmax(par$sigma, 1e-8)
    nu <- plogis(par$nu); tau <- plogis(par$tau)
    size <- 1/sigma
    
    p0 <- dnbinom(0, mu=mu, size=size)
    p1 <- dnbinom(1, mu=mu, size=size)
    denom <- pmax(1 - p0 - p1, 1e-10)
    
    out <- numeric(length(y))
    for(i in seq_along(y)) {
      yi <- y[i]
      if(yi < 0) { out[i] <- 0; next }
      if(yi < 0.9) { out[i] <- nu[i]; next }
      if(yi < 1.9) { out[i] <- nu[i] + (1-nu[i])*tau[i]; next }
      
      cdf_nb <- pnbinom(yi, mu=mu[i], size=size[i])
      base_accum <- p0[i] + p1[i]
      cdf_trunc <- (cdf_nb - base_accum) / denom[i]
      cdf_trunc <- max(0, min(1, cdf_trunc))
      
      out[i] <- nu[i] + (1-nu[i])*tau[i] + (1-nu[i])*(1-tau[i])*cdf_trunc
    }
    return(pmin(pmax(out, 1e-10), 1 - 1e-10))
  }
  
  sh_mu <- num_score_hess_factory(d_fun); sh_sigma <- num_score_hess_factory(d_fun)
  sh_nu <- num_score_hess_factory(d_fun); sh_tau <- num_score_hess_factory(d_fun)
  
  list(family="ZOINB", names=c("mu","sigma","nu","tau"), 
       links=c("log","log","logit","logit"), d=d_fun, p=p_fun,
       score=list(mu=function(y,p,...) sh_mu$score(y,p,id="mu"), sigma=function(y,p,...) sh_sigma$score(y,p,id="sigma"),
                  nu=function(y,p,...) sh_nu$score(y,p,id="nu"), tau=function(y,p,...) sh_tau$score(y,p,id="tau")),
       hess=list(mu=function(y,p,...) sh_mu$hess(y,p,id="mu"), sigma=function(y,p,...) sh_sigma$hess(y,p,id="sigma"),
                 nu=function(y,p,...) sh_nu$hess(y,p,id="nu"), tau=function(y,p,...) sh_tau$hess(y,p,id="tau")))
}

# RSB (Rescaled Beta Mixture) Implementation
f_RSB_mix <- function(a = 0.5, b = 0.5, n_quad = 60, grid_max = 50) {
  grid_vals <- exp(seq(log(0.001), log(grid_max), length.out = n_quad))
  
  dRSB_val <- function(eta) {
    eta <- pmax(eta, 1e-8)
    log_part <- log(1 + eta)
    num <- (log_part)^(a - 1)
    den <- (1 + eta) * (1 + log_part)^(a + b)
    val <- num / (den * beta(a, b))
    val[!is.finite(val)] <- 0
    return(val)
  }
  rsb_dens <- dRSB_val(grid_vals)
  
  grid_log <- log(grid_vals)
  w <- diff(grid_log)
  w <- c(w[1], w)
  w <- w * grid_vals 
  W_eff <- rsb_dens * w
  
  d_mix <- function(y, lambda, s) {
    lambda <- pmax(lambda, 1e-8)
    s <- plogis(s)
    d_pois_std <- dpois(y, lambda)
    
    int_val <- numeric(length(y))
    for(j in 1:length(grid_vals)) {
      lam_j <- lambda * grid_vals[j]
      dens_j <- dpois(y, lam_j)
      int_val <- int_val + dens_j * W_eff[j]
    }
    
    out <- (1 - s) * d_pois_std + s * int_val
    return(pmax(out, 1e-300))
  }
  
  d_fun <- function(y, par, log = FALSE) {
    dens <- d_mix(y, pmax(par$lambda, 1e-8), par$s)
    if (log) return(log(dens)) else return(dens)
  }
  
  p_fun <- function(y, par, ...) {
    n <- length(y)
    out <- numeric(n)
    
    for(i in 1:n) {
      yi <- y[i]
      if(yi < 0) { out[i] <- 0; next }
      
      lim_k <- min(yi, 300) 
      seq_k <- 0:lim_k
      
      lam_vec <- rep(par$lambda[i], length(seq_k))
      s_vec <- rep(par$s[i], length(seq_k))
      
      d_vals <- d_mix(seq_k, lam_vec, s_vec) 
      
      prob_accum <- sum(d_vals)
      if(yi > 300) prob_accum <- prob_accum + 1e-6 
      out[i] <- prob_accum
    }
    return(pmin(pmax(out, 1e-10), 1 - 1e-10))
  }
  
  sh_lambda <- num_score_hess_factory(d_fun); sh_s <- num_score_hess_factory(d_fun)
  
  list(family = "RSB", names = c("lambda", "s"), links = c("log", "logit"),
       d = d_fun, p = p_fun,
       score = list(lambda = function(y, par, ...) sh_lambda$score(y, par, id = "lambda"),
                    s      = function(y, par, ...) sh_s$score(y, par, id = "s")),
       hess = list(lambda = function(y, par, ...) sh_lambda$hess(y, par, id = "lambda"),
                   s      = function(y, par, ...) sh_s$hess(y, par, id = "s")))
}

# ==============================================================================
# 9. ALGORITHM FOR BAYESIAN MODEL FITTING
#    A wrapper function handling the optimization phase (posterior mode finding)
#    followed by MCMC sampling (GMCMC) with convergence checks.
# ==============================================================================

fit_robust_model <- function(label, formula_list, custom_fam_obj = NULL, df, 
                             n.iter = 6000, n.burn = 200, chains = 1, cores = 1,
                             thin = 10, step = 50, verbose = FALSE, 
                             clamp_limits = c(-10, 10), save_dir = "model_fits") {
  
  dir.create(save_dir, showWarnings = FALSE)
  message(sprintf("\n>>> Fitting Model: %s", label))
  
  if(!is.null(custom_fam_obj)) {
    fam_obj <- custom_fam_obj
  } else {
    if(label == "ZINB")        fam_obj <- gamlss.dist::ZINBI
    else if(label == "ZANBI")  fam_obj <- gamlss.dist::ZANBI
    else if(label == "PIG")    fam_obj <- gamlss.dist::PIG
    else if(label == "ZIPIG")  fam_obj <- gamlss.dist::ZIPIG
    else if(label == "DEL")    fam_obj <- gamlss.dist::DEL
    else if(label == "SICHEL") fam_obj <- gamlss.dist::SICHEL
    else stop(paste("Family not supported/mapped:", label))
  }
  
  message(" > Phase 1: Optimization (Stabilizing)...")
  opt <- tryCatch({
    bamlss(formula_list, family = fam_obj, data = df, sampler = FALSE, 
           optimizer = opt_bfit, control = list(maxit = 6000, trace = TRUE))
  }, error = function(e) { message("   Opt warn: ", e$message); return(NULL) })
  
  start_par <- NULL
  if(!is.null(opt)) {
    start_par <- coef(opt)
    start_par[is.na(start_par)] <- 0
    start_par <- clamp_vec(start_par, clamp_limits[1], clamp_limits[2])
    start_par <- start_par + rnorm(length(start_par), 0, 0.01) 
  }
  
  df_sub <- df %>% slice_sample(n = min(100, nrow(df)))
  try({
    if(!is.null(start_par)) {
    }
  }, silent = TRUE)
  
  message(" > Phase 2: MCMC Sampling...")
  
  run_mcmc <- function(starts) {
    tryCatch({
      bamlss(formula_list, family = fam_obj, data = df, 
             optimizer = FALSE, sampler = sam_GMCMC, start = starts,
             n.iter = n.iter, burnin = n.burn, thin = thin, 
             chains = chains, cores = cores, results = TRUE)
    }, error = function(e) return(NULL))
  }
  
  mod <- run_mcmc(start_par)
  
  if(is.null(mod)) {
    message("   Optimized start failed. Trying jittered start...")
    if(!is.null(start_par)) {
      mod <- run_mcmc(start_par + rnorm(length(start_par), 0, 0.1)) 
    }
  }
  
  if(is.null(mod)) {
    message("   Jitter failed. Trying random small start...")
    if(!is.null(start_par)) {
      start_rnd <- rnorm(length(start_par), 0, 0.05)
      names(start_rnd) <- names(start_par)
      mod <- run_mcmc(start_rnd)
    }
  }
  
  if(is.null(mod)) {
    message("   Random start failed. Trying cold start (NULL)...")
    mod <- run_mcmc(NULL) 
  }
  
  if(is.null(mod)) {
    message("   FATAL: Model failed to converge with all strategies.")
    return(NULL)
  }
  
  saveRDS(mod, file.path(save_dir, paste0(label, "_model.rds")))
  
  conv_report <- list(Rhat = NA, ESS = NA)
  try({
    smp <- bamlss::bamlss.get.samples(mod, combine = FALSE)

      if(!is.null(smp)) {
      mcmc_list <- as.mcmc.list(lapply(smp, as.mcmc))
      conv_report$ESS <- tryCatch(effectiveSize(mcmc_list), error=function(e) NA)
      if(chains > 1) {
        conv_report$Rhat <- tryCatch(gelman.diag(mcmc_list, autoburnin=FALSE)$psrf[,1], error=function(e) NA)
      }
    }
  }, silent = TRUE)
  
  mod$conv_report <- conv_report
  
  return(mod)
}

# ==============================================================================
# 10. MODEL SPECIFICATION AND ESTIMATION OF CANDIDATE DISTRIBUTIONAL FAMILIES
#     Defining the GAMLSS predictor formulas (Location, Scale, Shape) and iterating
#     through the candidate models (ZINB, ZIPIG, Sichel, etc.).
# ==============================================================================

f_mu_robust <- fatalities ~ 
  s(t_index, bs = "ps", k = 10) + 
  s(admin2, bs = "mrf", xt = list(penalty = penalty_mat)) + 
  ti(t_index, lon, lat, bs = c("ps", "tp"), d = c(1, 2), k = c(6, 6)) + 
  offset(df_model$log_pop) + classif_actor1 +
  s(dist_mueda_Airport_s, bs = "ps", k=10) +
  s(dist_afungi_lng_s, bs = "ps", k=10) +
  s(dist_mocimboa_port_s, bs = "ps", k=10) +
  financing_phase + precip_anomaly+bord+cyclone_exposure_6m

f_sigma <- ~ t_index + financing_phase
f_nu    <- ~ t_index + dist_afungi_lng_s
f_tau   <- ~ 1

f_2param <- list(mu = f_mu_robust, sigma = f_sigma)
f_3param <- list(mu = f_mu_robust, sigma = f_sigma, nu = f_nu)
f_4param <- list(mu = f_mu_robust, sigma = f_sigma, nu = f_nu, tau = f_tau)
f_rsb    <- list(lambda = f_mu_robust, s = f_nu)

ni <- 10000 
nb <- 2000
chains <- 1
cores <- 1

fitted_models <- list()

fitted_models$ZINB   <- fit_robust_model("ZINB", f_3param, NULL, df_model, ni, nb, chains, cores)
fitted_models$ZANBI  <- fit_robust_model("ZANBI", f_3param, NULL, df_model, ni, nb, chains, cores)
fitted_models$PIG    <- fit_robust_model("PIG", f_2param, NULL, df_model, ni, nb, chains, cores)
fitted_models$ZIPIG  <- fit_robust_model("ZIPIG", f_3param, NULL, df_model, ni, nb, chains, cores)
fitted_models$DEL    <- fit_robust_model("DEL", f_3param, NULL, df_model, ni, nb, chains, cores)
fitted_models$SICHEL <- fit_robust_model("SICHEL", f_3param, NULL, df_model, ni, nb, chains, cores)
fitted_models$ZOINB  <- fit_robust_model("ZOINB", f_4param, f_ZOINB_fam(), df_model, ni, nb, chains, cores)
fitted_models$RSB    <- fit_robust_model("RSB", f_rsb, f_RSB_mix(a=0.5, b=0.5, n_quad=60), df_model, ni, nb, chains, cores)

# ==============================================================================
# 11. COMPUTATION OF GOODNESS-OF-FIT CRITERIA AND RESIDUAL SPATIAL AUTOCORRELATION
#     Aggregating WAIC, DIC, Log-Likelihood, Convergence stats (Rhat/ESS), and 
#     performing Moran's I test on randomized quantile residuals.
# ==============================================================================

if(exists("nb_map")) {
  listw_obj <- nb2listw(nb_map, style = "W", zero.policy = TRUE)
} else {
  warning("nb_map not found. Moran's test will be skipped.")
}

get_all_metrics <- function(model_obj, model_name, df_data, listw) {
  
  message(paste("Calculating statistics for:", model_name))
  
  # A. Information Criteria and Fit
  waic_val <- tryCatch(bamlss::WAIC(model_obj)$WAIC1, error = function(e) NA)
  dic_val  <- tryCatch(bamlss::DIC(model_obj)$DIC, error = function(e) NA)
  log_lik  <- logLik(model_obj)
  df_model <- attr(log_lik, "df")
  
  # B. Convergence (Rhat and ESS) If you need change the number of "chains <- 1 and cores <- 1"
  rhat <- NA; ess <- NA
  if(!is.null(model_obj$conv_report)) {
    if(any(!is.na(model_obj$conv_report$Rhat))) rhat <- max(model_obj$conv_report$Rhat, na.rm=TRUE)
    if(any(!is.na(model_obj$conv_report$ESS))) ess <- min(model_obj$conv_report$ESS, na.rm=TRUE)
  }
  
  # C. Residuals and Outliers
  resids <- residuals(model_obj, type = "quantile")
  
  if(sum(is.na(resids)) > length(resids)/2) {
    return(NULL) 
  }
  
  n_outliers <- sum(abs(resids) > 3, na.rm = TRUE)
  perc_outliers <- (n_outliers / length(resids)) * 100
  
  # D. Spatial Autocorrelation (Moran's I)
  
  df_res <- data.frame(admin2 = df_data$admin2, resid = resids)
  df_agg <- df_res %>%
    group_by(admin2) %>%
    summarise(mean_resid = mean(resid, na.rm = TRUE)) %>%
    arrange(match(admin2, attributes(listw)$region.id)) 
  
  moran_res <- tryCatch({
    moran.test(df_agg$mean_resid, listw, alternative = "two.sided", zero.policy = TRUE)
  }, error = function(e) list(estimate = c("Moran I statistic" = NA), p.value = NA))
  
  return(data.frame(
    Model = model_name,
    WAIC = waic_val,
    DIC = dic_val,
    LogLik = as.numeric(log_lik),
    EDF = df_model, 
    Max_Rhat = rhat,
    Min_ESS = ess,
    Outliers_N = n_outliers,
    Outliers_Pct = perc_outliers,
    Moran_I = as.numeric(moran_res$estimate[1]),
    Moran_P_Val = as.numeric(moran_res$p.value)
  ))
}

results_list <- list()
for(name in names(fitted_models)) {
  if(!is.null(fitted_models[[name]])) {
    results_list[[name]] <- get_all_metrics(fitted_models[[name]], name, df_model, listw_obj)
  }
}

final_table <- do.call(rbind, results_list) %>%
  arrange(WAIC) 

print(final_table)
write.csv(final_table, "model_selection_results.csv", row.names = FALSE)

kbl_latex <- kableExtra::kbl(final_table, format = "latex", booktabs = TRUE, digits = 3,
                             caption = "Comparison of fit statistics, convergence diagnostics, and residual spatial autocorrelation.") %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))
print(kbl_latex)


# ==============================================================================
# 12. QUANTILE RESIDUALS, WORM PLOTS, AND INFLUENCE ANALYSIS
# ==============================================================================

generate_diagnostic_plots <- function(model_list) {
  plot_list <- list()
  
  for(name in names(model_list)) {
    mod <- model_list[[name]]
    if(is.null(mod)) next
    
    p_diagnostic <- tryCatch({
      res <- residuals(mod, type = "quantile")
      res <- res[is.finite(res)]
      
      if(length(res) > 5000) res <- sample(res, 5000)
      
      n <- length(res)
      res_sorted <- sort(res)
      theo_q <- qnorm(ppoints(n))
      
      nsim <- 1000
      sim_mat <- matrix(rnorm(n * nsim), ncol = nsim)
      sim_mat <- apply(sim_mat, 2, sort)
      lower <- apply(sim_mat, 1, quantile, probs = 0.025)
      upper <- apply(sim_mat, 1, quantile, probs = 0.975)
      
      df_qq <- data.frame(Obs = res_sorted, Theo = theo_q, Low = lower, Upp = upper)
      
      # 1. QQ-Plot
      p1 <- ggplot(df_qq, aes(x = Theo, y = Obs)) +
        geom_ribbon(aes(ymin = Low, ymax = Upp), fill = "grey85", alpha = 0.5) +
        geom_point(size = 0.5, alpha = 0.6) +
        geom_abline(color = "#e74c3c", linetype = "dashed") +
        labs(title = paste0(name, ": QQ-Plot"), x = "Theoretical", y = "Quantile Residuals") + 
        theme_minimal()
      
      # 2. Worm Plot (Deviation)
      p2 <- ggplot(df_qq, aes(x = Theo, y = Obs - Theo)) +
        geom_hline(yintercept = 0, color = "#e74c3c") +
        geom_ribbon(aes(ymin = Low - Theo, ymax = Upp - Theo), fill = "grey85", alpha = 0.5) +
        geom_point(size = 0.5, alpha = 0.6) +
        geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "#3498db", se = FALSE, linewidth = 0.5) +
        labs(title = "Worm Plot", x = "Theoretical", y = "Deviation") + 
        ylim(-2, 2) + 
        theme_minimal()
      
      # 3. Influence Plot
      df_inf <- data.frame(Index = 1:n, Resid = res)
      outliers <- df_inf %>% filter(abs(Resid) > 3)
      
      p3 <- ggplot(df_inf, aes(x = Index, y = Resid)) +
        geom_hline(yintercept = c(-3, 3), color = "#e74c3c", linetype = "dotted") +
        geom_point(alpha = 0.4, size = 0.5) +
        geom_point(data = outliers, color = "#e74c3c", size = 1.2) +
        labs(title = "Influence Plot", x = "Index", y = "Residuals") + 
        theme_minimal()
      
      (p1 + p2 + p3)
      
    }, error = function(e) {
      ggplot() + annotate("text", x=0.5, y=0.5, label="Plot Error") + theme_void()
    })
    
    plot_list[[name]] <- p_diagnostic
  }
  return(plot_list)
}

diag_plots <- generate_diagnostic_plots(fitted_models)

if(length(diag_plots) > 0) {
  if(length(diag_plots) <= 4) {
    grid_all <- wrap_plots(diag_plots, ncol = 1)
    ggsave(here("figures","Figure_Diagnostics_All.pdf"), grid_all, width = 14, height = 3.5 * length(diag_plots), limitsize = FALSE)
    print(grid_all)
  } else {
    grid_1 <- wrap_plots(diag_plots[1:4], ncol = 1)
    ggsave(here("figures","Figure_Diagnostics_Part1.pdf"), grid_1, width = 14, height = 14, limitsize = FALSE)
    print(grid_1)
    
    grid_2 <- wrap_plots(diag_plots[5:length(diag_plots)], ncol = 1)
    ggsave(here("figures","Figure_Diagnostics_Part2.pdf"), grid_2, width = 14, height = 3.5 * (length(diag_plots)-4), limitsize = FALSE)
    print(grid_2)
  }
}


# ==============================================================================
# 13. SPATIAL RESIDUAL DIAGNOSTICS TO DETECT UNMODELLED CLUSTERING
# ==============================================================================

generate_spatial_maps <- function(model_list, data, shapefile) {
  map_list <- list()
  
  for (name in names(model_list)) {
    mod <- model_list[[name]]
    if(is.null(mod)) next
    
    res <- tryCatch(residuals(mod, type = "quantile"), error = function(e) NULL)
    
    if(!is.null(res)) {
      data$temp_resid <- res
      df_agg <- data %>% 
        group_by(admin2) %>% 
        summarise(mean_resid = mean(temp_resid, na.rm = TRUE)) %>% 
        ungroup()
      
      df_agg$plot_val <- pmax(pmin(df_agg$mean_resid, 3), -3)
      map_data <- left_join(shapefile, df_agg, by = c("NAME_2"="admin2"))
      
      p_map <- ggplot(map_data) +
        geom_sf(aes(fill = plot_val), color = "black", size = 0.1) +
        scale_fill_gradient2(low = "#2980b9", mid = "white", high = "#c0392b", 
                             midpoint = 0, name = "Avg Resid",
                             limits = c(-3, 3), oob = scales::squish) +
        labs(title = name) +
        theme_void() +
        theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face="bold"))
      
      map_list[[name]] <- p_map
    }
  }
  return(map_list)
}

spatial_maps <- generate_spatial_maps(fitted_models, df_model, shapefile)

if(length(spatial_maps) > 0) {
  final_sp_grid <- wrap_plots(spatial_maps, ncol = 4) + 
    plot_annotation(title = "Spatial Residual Autocorrelation Maps",
                    subtitle = "Checking for unmodelled spatial clustering (Red/Blue clusters).",
                    caption = "Values are district-averaged randomized quantile residuals.")
  
  print(final_sp_grid)
  
  h_val <- ceiling(length(spatial_maps)/4) * 4
  ggsave(here("figures","Figure_Spatial_Residuals_Maps.pdf"), final_sp_grid, width = 14, height = h_val, limitsize = FALSE)
}


# ==============================================================================
# 14. CALIBRATION ASSESSMENT VIA PROBABILITY INTEGRAL TRANSFORM (PIT)
#     Checking the uniformity of PIT values to assess predictive performance.
# ==============================================================================

get_pit_plots_separate <- function(model, label) {
  
  rqr <- tryCatch({
    r <- residuals(model, type = "quantile")
    if(sum(is.finite(r)) < 10) stop("Insufficient finite residuals")
    r[is.finite(r)]
  }, error = function(e) NULL)
  
  if(is.null(rqr)) {
    dummy <- ggplot() + 
      annotate("text", x=0.5, y=0.5, label="PIT Failed", size=3) + 
      labs(title = label) + theme_void() + 
      theme(panel.border = element_rect(colour = "grey", fill=NA, size=1))
    return(list(hist = dummy, ecdf = dummy))
  }
  
  pit <- pnorm(rqr)
  
  p_hist <- ggplot(data.frame(pit=pit), aes(x=pit)) +
    geom_histogram(aes(y=after_stat(density)), bins = 15, 
                   fill="#2c3e50", color="white", alpha=0.8) +
    geom_hline(yintercept = 1, linetype="dashed", color="#c0392b") +
    labs(title = label, x = "PIT", y = "Density") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face="bold", hjust=0.5))
  
  n <- length(pit)
  y <- (1:n)/n
  x <- sort(pit)
  k <- 1.36 / sqrt(n)
  
  df_ecdf <- data.frame(x=x, y=y, lower=pmax(0, y-k), upper=pmin(1, y+k))
  
  p_ecdf <- ggplot(df_ecdf, aes(x=x, y=y)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey80", alpha=0.5) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="#c0392b") +
    geom_line(color="#2980b9", linewidth=0.8) +
    labs(title = label, x = "Theoretical", y = "Empirical") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face="bold", hjust=0.5)) +
    coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1))
  
  return(list(hist = p_hist, ecdf = p_ecdf))
}

list_histograms <- list()
list_ecdfs <- list()

for (model_name in names(fitted_models)) {
  if (is.null(fitted_models[[model_name]])) next
  
  plots <- get_pit_plots_separate(fitted_models[[model_name]], model_name)
  
  list_histograms[[model_name]] <- plots$hist
  list_ecdfs[[model_name]]      <- plots$ecdf
}

if(length(list_histograms) > 0) {
  grid_hist <- wrap_plots(list_histograms, ncol = 4) + 
    plot_annotation(
      title = "",
      subtitle = "Flat bars close to the red line indicate good calibration."
    )
  
  ggsave(here("figures","Figure_PIT_Histograms.pdf"), grid_hist, width = 14, height = 7)
  print(grid_hist)
}

if(length(list_ecdfs) > 0) {
  grid_ecdf <- wrap_plots(list_ecdfs, ncol = 4) + 
    plot_annotation(
      title = "",
      subtitle = "Blue line should remain within the grey confidence bands."
    )
  
  ggsave(here("figures","Figure_PIT_ECDFs.pdf"), grid_ecdf, width = 14, height = 7) 
  print(grid_ecdf)
}


# ==============================================================================
# 15. POSTERIOR ESTIMATES OF FIXED EFFECTS FOR THE SELECTED MODEL (ZIPIG)
#     Extracting and visualizing the coefficients for the best performing model.
# ==============================================================================

final_model_obj <- fitted_models$ZIPIG 

get_fixed_effects <- function(model_entry) {
  mod <- model_entry
  summ <- summary(mod)
  
  process_param <- function(matrix_coef, param_name) {
    if(is.null(matrix_coef)) return(NULL)
    df <- as.data.frame(matrix_coef)
    df$Term <- rownames(df)
    df$Parameter <- param_name
    colnames(df)[which(names(df) == "2.5%")] <- "Lower"
    colnames(df)[which(names(df) == "97.5%")] <- "Upper"
    colnames(df)[which(names(df) == "Mean")] <- "Estimate"
    df <- df %>% filter(!grepl("Intercept", Term))
    return(df[, c("Term", "Estimate", "Lower", "Upper", "Parameter")])
  }
  
  all_fixed <- list()
  for(p in names(summ$model.matrix)) {
    all_fixed[[p]] <- process_param(summ$model.matrix[[p]], p)
  }
  
  return(bind_rows(all_fixed))
}

df_fixed <- get_fixed_effects(final_model_obj)


df_fixed<- df_fixed %>%
  mutate(Term = as.character(Term)) %>%
  mutate(Term = case_when(
    str_detect(Term, "classif_actor1") ~ str_replace(Term, "classif_actor1", "Actor: "),
    
    str_detect(Term, "financing_phase2") ~ "Financing: Transition Phase",
    str_detect(Term, "financing_phase3") ~ "Financing: Regional Phase",
    
    Term == "precip_anomaly" ~ "Precipitation Anomaly",
    Term == "t_index" ~ "Time Trend Index",
    Term == "dist_afungi_lng_s" ~ "Dist. to LNG (Standardized)",
    Term == "alpha" ~ "Alpha (Shape Parameter)",
    
    TRUE ~ Term
  )) %>%
  mutate(Term = str_remove(Term, "\\.{3}\\d+"))


g_forest <- ggplot(df_fixed, aes(x = Estimate, y = reorder(Term, Estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, color = "#2c3e50") +
  geom_point(size = 2.5, color = "#c0392b") +
  facet_wrap(~Parameter, scales = "free_x") +
  labs(title = "",
       subtitle = "Mean and 95% Credible Intervals (Log Scale)",
       x = "Coefficient Estimate", y = NULL) +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "gray90"), 
        axis.text.y = element_text(face = "bold"),
        strip.background = element_rect(fill = "gray95"))

ggsave(here("figures","Figure_Fixed_Effects_Forest.pdf"), g_forest, width = 12, height = 7)


# ==============================================================================
# 16. VISUALIZATION OF NON-LINEAR SMOOTH EFFECTS ON THE LOCATION PARAMETER
#     Plotting the partial effects of penalized splines (P-splines) for key covariates.
# ==============================================================================

get_smooth_effect <- function(model, dataset, term_full, var_name, parameter) {
  grid_vals <- seq(min(dataset[[var_name]], na.rm = TRUE), 
                   max(dataset[[var_name]], na.rm = TRUE), 
                   length.out = 100)
  
  nd <- data.frame(placeholder = grid_vals)
  names(nd) <- var_name 
  
  set.seed(123)
  preds <- predict(model, newdata = nd, model = parameter, term = term_full, 
                   FUN = function(x) { 
                     c(mean = mean(x), 
                       lower = quantile(x, 0.025), 
                       upper = quantile(x, 0.975)) 
                   },
                   intercept = FALSE) 
  
  pred_mat <- as.data.frame(preds)
  cols <- colnames(pred_mat)
  
  df_out <- data.frame(
    Value = grid_vals,
    eff   = pred_mat[, grep("mean", cols)],
    lower = pred_mat[, grep("lower", cols)],
    upper = pred_mat[, grep("upper", cols)],
    Term  = term_full,
    Parameter = parameter
  )
  
  return(df_out)
}

map_terms_mu <- list(
  "s(dist_mueda_Airport_s)"   = "dist_mueda_Airport_s",
  "s(dist_afungi_lng_s)"      = "dist_afungi_lng_s",
  "s(dist_mocimboa_port_s)"   = "dist_mocimboa_port_s"
)

list_mu <- list()

for(term in names(map_terms_mu)) {
  message("Attempting extraction for: ", term)
  
  res <- try(get_smooth_effect(model = final_model_obj, 
                               dataset = df_model, 
                               term_full = term, 
                               var_name = map_terms_mu[[term]], 
                               parameter = "mu"), silent = TRUE)
  
  if(!inherits(res, "try-error")) {
    list_mu[[term]] <- res
  } else {
    message("Failed to extract: ", term)
    print(attr(res, "condition"))
  }
}

df_smooths_all <- bind_rows(list_mu)

if(nrow(df_smooths_all) > 0) {
  labeller_vars <- c(
    "s(dist_mueda_Airport_s)" = "Dist. Mueda (Military)",
    "s(dist_afungi_lng_s)"    = "Dist. Afungi (LNG Gas)",
    "s(dist_mocimboa_port_s)" = "Dist. Mocímboa (Port)"
  )
  
  g_smooths <- ggplot(df_smooths_all, aes(x = Value, y = eff)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2980b9", alpha = 0.3) +
    geom_line(color = "#2980b9", linewidth = 1) +
    facet_wrap(~ Term, scales = "free", labeller = labeller(Term = labeller_vars), ncol = 3) +
    labs(x = "Standardized Distance (Z-Score)", 
         y = "Partial Effect on Log-Intensity",
         caption = "Confidence bands represent 95% Credible Intervals.") +
    theme_bw() +
    theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
  
  print(g_smooths)
} else {
  message("DataFrame is still empty. Check if df_model contains the standardized variables.")
}

ggsave(here("figures", "Figure_NonLinear_Effects_Mu_Final.pdf"), 
       plot = g_smooths, 
       width = 14, 
       height = 5, 
       device = "pdf",
       units = "in",
       dpi = 300
)

# ==============================================================================
# 17. SPATIO-TEMPORAL EVOLUTION OF EXPECTED FATALITY RISK
#     Mapping the posterior expected fatalities across years to visualize the 
#     diffusion of conflict intensity.
# ==============================================================================

model_obj <- fitted_models$ZIPIG
family_name <- model_obj$family$family

message("Calculating Expected Value for Family: ", family_name)
preds <- predict(model_obj, type = "parameter")
mu_vec <- as.numeric(preds$mu)
nu_vec <- as.numeric(preds$nu)

if(any(is.infinite(mu_vec))) {
  max_val <- max(mu_vec[is.finite(mu_vec)], na.rm = TRUE)
  mu_vec[is.infinite(mu_vec)] <- max_val
}

if (family_name %in% c("ZIPIG", "ZINB", "ZANBI")) {
  expected_val <- (1 - nu_vec) * mu_vec
} else {
  expected_val <- mu_vec
}

if (length(expected_val) == nrow(df_model)) {
  df_model$monthly_risk <- expected_val
} else {
  stop("Length mismatch: Predictions do not match data rows.")
}

message("Summary of Monthly Expected Fatalities:")
print(summary(df_model$monthly_risk))


df_spatio <- df_model %>%
  mutate(year = year(date_obj)) %>% 
  group_by(admin2, year) %>%
  summarise(
    total_risk = sum(monthly_risk, na.rm = TRUE), 
    .groups = 'drop'
  )

full_grid <- expand.grid(
  admin2 = unique(shape_subset$NAME_2), 
  year = sort(unique(df_spatio$year))
)

df_spatio_full <- left_join(full_grid, df_spatio, by = c("admin2", "year"))
df_spatio_full$total_risk[is.na(df_spatio_full$total_risk)] <- 0

limite_visual <- quantile(df_spatio_full$total_risk, 0.99, na.rm=TRUE)
if(limite_visual < 5) limite_visual <- max(df_spatio_full$total_risk, na.rm=TRUE)

df_spatio_full$plot_risk <- pmin(df_spatio_full$total_risk, limite_visual)


map_data <- shape_subset %>%
  inner_join(df_spatio_full, by = c("NAME_2" = "admin2"))
g_evolucao <- ggplot(map_data) +
  geom_sf(data = shapefile, fill = "white", color = "gray95") +
  geom_sf(data = shape_subset, fill = "white", color = "gray95") +
  geom_sf(aes(fill = plot_risk), color = "gray80", size = 0.1) +
  geom_sf_text(aes(label = NAME_2), size = 1.5, color = "black", alpha = 0.5, check_overlap = TRUE) +
  
  scale_fill_gradientn(
    colors = c("white", "#feb24c", "#fd8d3c", "#fc4e2a", "#b10026"),
    limits = c(0, limite_visual),
    name = "Expected\nFatalities",
    oob = scales::squish 
  ) +
  
  annotation_scale(
    location = "br",            
    width_hint = 0.3,           
    bar_cols = c("black", "white"),
    text_family = "sans"        
  ) +
  annotation_north_arrow(
    location = "br",            
    which_north = "true",       
    pad_x = unit(0.1, "in"),    
    pad_y = unit(0.2, "in"), 
    height = unit(.8, "cm"),
    width = unit(.8, "cm"),
    style = north_arrow_fancy_orienteering 
  ) +
  
  facet_wrap(~year, ncol = 4) +
  
  labs(
    title = "",
    subtitle = "Annual Aggregated Expected Fatalities (Posterior Mean)",
    caption = paste0("Note: Scale capped at ", round(limite_visual, 1), " fatalities for visualization clarity.")
  ) +
  
  theme_void() +
  theme(
    strip.text = element_text(face = "bold", size = 11, margin = margin(b=5)),
    legend.position = "right",
    legend.title = element_text(size=10, face="bold"),
    plot.title = element_text(face="bold", hjust=0.5),
    plot.subtitle = element_text(hjust=0.5, color="gray40"),
    plot.background = element_rect(fill="white", color=NA) 
  )

print(g_evolucao)


ggsave(here("figures","Figure_SpatioTemporal_Evolution_ExpectedFatalities.pdf"), 
       g_evolucao, width = 14, height = 9, bg = "white")

