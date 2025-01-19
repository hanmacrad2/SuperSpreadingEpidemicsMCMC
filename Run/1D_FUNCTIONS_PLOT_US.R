#PLOT MPOX BY STATE
library(dplyr)

#************************************************
#
# PLOTTING FUNCTIONS
#
#************************************************

format_df_usa_sd_dates <- function(df_mpox){
  
  'Format_df_usa for San Diego dates'
  
  df_result = data.frame()
  list_states = unique(df_mpox$Jurisdiction)
  
  for (i in seq_along(list_states)) {
    
    df_state <- df_mpox[df_mpox$Jurisdiction == list_states[i], ]
    df_state <- df_state[-seq_len(min(3, nrow(df_state))), ] 
    
    df_result <- rbind(df_result, df_state)
  }
  
  return(df_result)
}

#PLOT TOP STATES 
plot_top_states <- function(df_mpox, data_la, data_sd, lwd = 1.5, cex_size = 1.5,
                            legend_location = 'topright', PLOT_ALL = FALSE){
  
  #Plotting params
  list_states = get_list_top_states(df_mpox)
  list_colors = get_list_colors(list_states)
  if(PLOT_ALL){
    title = bquote(paste("Mpox Reported Cases across 56 States & Cities"))
  } else {
    title = bquote(paste("Mpox Reported Cases across Top 10 States & Cities"))
  }
  
  for (i in seq_along(list_states)) {
    
    if(i == 1){
      
      df_state <- df_mpox[df_mpox$Jurisdiction == list_states[i], ]
      
      plot(NA, xlim = c(0, length(df_state$Cases)), 
           ylim = c(0, max(df_state$Cases)),
           type = "n", xlab = "Time", ylab = "Case Count",
           main = title, 
           xaxt = "n", 
           cex.lab = cex_size, cex.axis = cex_size, cex.main = cex_size+0.2)
    } 
    
    df_state <- df_mpox[df_mpox$Jurisdiction == list_states[i], ]
    df_state$timepoint = seq_along(df_state$Cases)
    
    #Plot
    lines(df_state$timepoint, df_state$Cases, type = "o",
          pch = 21, bg = list_colors[i],                    # Use filled circles
          col = list_colors[i], lwd = lwd)

  }

  #X-axis
  time_points = seq(from = 4, to = nrow(df_state), by = 4) 
  df_labels = subset(df_state, timepoint %in% time_points)
  axis(1, at = time_points, labels = df_labels$month_year, las = 1, cex.axis = 0.8)
  
  if(PLOT_ALL){
    list_bottom_states = find_bottom_states(df_mpox)
    plot_bottom_states(list_bottom_states, df_mpox)
  }
  
  #PLOT LA AND SAN DIEGO
  plot_mpox_la_sd(data_la, data_sd)
  
  #Legend
  legend(legend_location, 
         legend = c(list_states, 'San Diego'), 
         fill = c(list_colors, 'black'),
         title = "Top 10 by Cases + San Diego", 
         cex = 0.95, 
         bty = "n")
}

#plot top states
plot_top_states_filtered <- function(df_mpox, df_mpox_filtered, data_la, data_sd, lwd = 1.5, cex_size = 1.5,
                            legend_location = 'topright', PLOT_ALL = FALSE){
  
  #Plotting params
  list_states = get_list_top_states(df_mpox)
  list_colors = get_list_colors(list_states)
  if(PLOT_ALL){
    title = bquote(paste("Mpox Reported Cases across 56 States & Cities"))
  } else {
    title = bquote(paste("Mpox Reported Cases across Top 10 States & Cities"))
  }
  
  for (i in seq_along(list_states)) {
    
    if(i == 1){
      
      df_state <- df_mpox_filtered[df_mpox_filtered$Jurisdiction == list_states[i], ]
      
      plot(NA, xlim = c(0, length(df_state$Cases)), 
           ylim = c(0, max(df_state$Cases)),
           type = "n", xlab = "Time", ylab = "Case Count",
           main = title, 
           xaxt = "n", 
           cex.lab = cex_size, cex.axis = cex_size, cex.main = cex_size+0.2)
    } 
    
    df_state <- df_mpox_filtered[df_mpox_filtered$Jurisdiction == list_states[i], ]
    df_state$timepoint = seq_along(df_state$Cases)
    
    #Plot
    lines(df_state$timepoint, df_state$Cases, type = "o",
          pch = 21, bg = list_colors[i],                    # Use filled circles
          col = list_colors[i], lwd = lwd)
    
  }
  
  #X-axis
  time_points = seq(from = 4, to = nrow(df_state), by = 4) 
  df_labels = subset(df_state, timepoint %in% time_points)
  axis(1, at = time_points, labels = df_labels$month_year, las = 1, cex.axis = 0.8)
  
  if(PLOT_ALL){
    list_bottom_states = find_bottom_states(df_mpox)
    plot_bottom_states(list_bottom_states, df_mpox_filtered)
  }
  
  #PLOT LA AND SAN DIEGO
  plot_mpox_la_sd(data_la, data_sd)
  
  #Legend
  legend(legend_location, 
         legend = c(list_states, 'San Diego'), 
         fill = c(list_colors, 'black'),
         title = "Top 10 by Cases + San Diego", 
         cex = 0.95, 
         bty = "n")
}


get_top_states <- function(df_mpox) {
  
  list_states = find_top_states(df_mpox)
  
  # Filter the dataframe for the selected jurisdictions
  df_mpox_filtered <- df_mpox %>%
    filter(Jurisdiction %in% list_states)
  
  return(df_mpox_filtered)
  
}

get_list_top_states <- function(df_mpox){
  
  list_states = find_top_states(df_mpox)
  
  list_states_final = list()
  list_states_final[c(1:4)] =  list_states[c(1:4)]
  list_states_final[5] = 'LA'
  list_states_final[c(6:10)] =  list_states[c(5:9)]
  
  #print(list_states_final)
  return(list_states_final)
}


get_list_colors <- function(list_states) {
  # Color assignments
  color_map <- list(
    "California" = "darkgreen",
    "New York City" = "red",
    "Texas" = "#8A2BE2",  # Bright Blue Violet
    "Florida" = "#FFFF00",    # Bright Yellow
    "LA" = "blue",          # LA
    "Georgia" = "#FF1493",  # Bright Green  ""
    "Illinois" = "#00FFFF", # Bright Cyan
    "New Jersey" = "#FFA500", # Bright Blue
    "North Carolina" = "#FF6347",
    "Maryland" = "#00FF00", # Bright Orange 
    "San Diego" = "black"   # San Diego
  )
  
  # Generate colors in order of list_states_final
  list_colors <- sapply(list_states, function(state) color_map[[state]])
  return(list_colors)
}

#Find top states
find_top_states <- function(df_mpox, n_top = 9){
 
  top_jurisdictions <- df_mpox %>%
    group_by(Jurisdiction) %>%
    summarise(Total_Cases = sum(Cases, na.rm = TRUE)) %>%
    arrange(desc(Total_Cases)) %>%
    slice_head(n = n_top) # Get the top 10 jurisdictions 
  
  print(top_jurisdictions)
  
  #> sum(data_la$Cases): 2599
  #sum(data_sd$Cases):602
  
  return(top_jurisdictions$Jurisdiction)
}

#Find top states
find_bottom_states <- function(df_mpox, n_top = 9){
  
  n_states = length(unique(df_mpox$Jurisdiction))
  
  bottom_jurisdictions <- df_mpox %>%
    group_by(Jurisdiction) %>%
    summarise(Total_Cases = sum(Cases, na.rm = TRUE)) %>%
    arrange(desc(Total_Cases)) %>%
    slice_tail(n = n_states -n_top) # Get the top 10 jurisdictions 
  
  #(bottom_jurisdictions)
  
  return(bottom_jurisdictions$Jurisdiction)
}

#plot_bottom_states
plot_bottom_states <-function(list_bottom_states, df_mpox, lwd = 1.5){
  
  #Find
  list_colors_extra = sample(rainbow(length(list_bottom_states)))
  
  for (i in seq_along(list_bottom_states)) {
    
    print(paste0('State:', list_bottom_states[i]))
    df_state_bottom <- df_mpox[df_mpox$Jurisdiction == list_bottom_states[i], ]
    df_state_bottom$timepoint = seq_along(df_state_bottom$Cases)
    
    #Plot
    lines(df_state_bottom$timepoint, df_state_bottom$Cases, type = "o",
          pch = 21, bg = list_colors_extra[i],                    # Use filled circles
          col = list_colors_extra[i], lwd = lwd)
    
    #vec_cases = data_top_10_subset$Cases
  }
  
}


#Plot_mpox_la_sd
plot_mpox_la_sd <- function(data_la, data_sd, col_la = 'blue', col_sd = 'black'){
  
  #LA
  data_la$timepoint = seq_along(data_la$Cases)
  lines(data_la$timepoint, data_la$Cases, type = "o",
        pch = 21, bg = col_la,                    
        col = col_la, lwd = 2.5)
  
  #SAN DIEGO
  data_sd$timepoint = seq_along(data_sd$Cases)
  lines(data_sd$timepoint, data_sd$Cases, type = "o",
        pch = 21, bg = col_sd,                    
        col = col_sd, lwd = 4)
}

#**********
#* Apply
#Plot
df_mpox = format_df_usa_sd_dates(df_mpox)

#Filtered Data
#df_mpox_filtered = get_selected_states(df_mpox)

#Plot
plot_top_states(df_mpox, data_la, data_sd)

#Plot all
plot_top_states(df_mpox, data_la, data_sd, PLOT_ALL = TRUE)


#********************************************************
#FILTERED DATA
#********************

# Filter data for dates after January 1, 2023
df_mpox_filtered <- df_mpox %>%
  filter(date > as.Date("2023-01-01"))

df_la_filtered <- data_la %>%
  filter(date > as.Date("2023-01-01"))

df_sd_filtered <- data_sd %>%
  filter(date > as.Date("2023-01-01"))

#Plot
plot_top_states_filtered(df_mpox, df_mpox_filtered, df_la_filtered, df_sd_filtered,
                         legend_location = 'topleft')

#Plot
plot_top_states_filtered(df_mpox, df_mpox_filtered, df_la_filtered, df_sd_filtered,
                         legend_location = 'topleft', PLOT_ALL = TRUE)

#PLot NYC, CALI, SD
df_nyc = subset(df_mpox, Jurisdiction == 'New York City')
df_nyc_filtered <- df_nyc %>%
  filter(date > as.Date("2023-01-01"))
df_nyc_filtered$timepoint = seq_along(df_nyc_filtered$Cases)
