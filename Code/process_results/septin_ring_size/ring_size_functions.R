# ----------------------------------------------------------------------------------------------------
# Functions for computing ring-size
# ----------------------------------------------------------------------------------------------------
# Note either (x1, y1, z1) or (x2, y2, z2) can be vectors.
compute_sphere_distance <- function(x1, y1, z1, x2, y2, z2, R)
{
  dist <- R * acos((x1*x2 + y1*y2 + z1*z2) / R^2)
  return(dist)
}


get_test_points <- function(x, y, z, d_in, d_out, xc, yc, zc)
{
  out <- rep(0, length(x))
  for(i in 1:length(x)){
    dist <- compute_sphere_distance(x[i], y[i], z[i], xc, yc, zc,  1.0)
    if(dist > d_in && dist < d_out){
      out[i] <- 1
    }else{
      out[i] <- 0
    }
  }
  return(out)
}


# Get the minium distance between two nodes, used to figure out where the "hole" begins in a linesearch 
get_min_distace <- function(data){
  min_dist <- Inf
  for(i in 1:(nrow(data))){
    xref <- data$x[i]
    yref <- data$y[i]
    zref <- data$z[i]
    
    xtest <- data$x[-i]
    ytest <- data$y[-i]
    ztest <- data$z[-i]
    
    dist <- min(compute_sphere_distance(xref, yref, zref, xtest, ytest, ztest, 1.0))
    if(is.infinite(min_dist)){
      min_dist <- dist
    }else if(dist > min_dist){
      min_dist <- dist
    }
  }
  return(min_dist)
}

get_max_distace <- function(data, iref){
  
  xref <- data$x[iref]
  yref <- data$y[iref]
  zref <- data$z[iref]
  
  xtest <- data$x[-iref]
  ytest <- data$y[-iref]
  ztest <- data$z[-iref]
  
  dist <- compute_sphere_distance(xref, yref, zref, xtest, ytest, ztest, 1.0)
  distMax <- max(dist)
  distArgMax <- which.max(dist)
  xmax <- xtest[distArgMax]
  ymax <- ytest[distArgMax]
  zmax <- ztest[distArgMax]
  
  return(tibble(max_dist = distMax, x=xmax, y=ymax, z=zmax))
}


get_dist_ring <- function(data)
{
  
  min_dist <- get_min_distace(data)
  eps <- min_dist*0.1
  t <- seq(from = 1e-6, to = 1-1e-6, length.out=10000)
  data_save <- tibble()
  for(i in 1:nrow(data)){
    
    res <- get_max_distace(data, i)
    
    P <- c(res$x[1], res$y[1], res$z[1])
    Q <- c(data$x[i], data$y[i], data$z[i])
    for(j in 1:length(t)){
      t_use <- t[j]
      R <- P + t_use*(Q - P) 
      R <- R / sqrt(sum(R^2))
      dist_tmp <- min(compute_sphere_distance(R[1], R[2], R[3], data$x, data$y, data$z, 1.0))
      if(dist_tmp  > min_dist + eps){
        dist_min <- compute_sphere_distance(R[1], R[2], R[3], data$x[i], data$y[i], data$z[i], 1.0)
        break
      }
    }
    
    data_save <- data_save |> 
      bind_rows(tibble(x=data$x[i], y=data$y[i], z=data$z[i], dmax = res$max_dist[1], dmin = dist_min))
  }
  return(data_save)
}


fibonacci_sphere <- function(samples=1000)
{
  
  x <- rep(0, samples)
  y <- rep(0, samples)
  z <- rep(0, samples)
  phi <- pi * (3. - sqrt(5.))  # golden angle in radians
  
  for(i in 0:(samples-1)){
    y[i] = 1 - (i / (samples - 1)) * 2  # y goes from 1 to -1
    print(sprintf("y[%d] = %.3f", i, y[i]))
    radius = sqrt(1 - y[i] * y[i])  # radius at y
    theta = phi * i  # golden angle increment
    x[i] = cos(theta) * radius
    z[i] = sin(theta) * radius
  }
  return(tibble(x=x, y=y, z=z))
} 


# Find the point c which is a bit further away from pc than b, and lies on the line between pc a b 
find_point_c <- function(pc, b)
{
  
  theta <- sum(pc*b)
  phi <- ifelse(theta > 0.0, 0.02, -0.02)
  n <- c(pc[2]*b[3] - pc[3]*b[2], -(pc[1]*b[3] - pc[3]*b[1]), pc[1]*b[2] - pc[2]*b[1]) # Cross-product to get normal
  
  A1 <- cos(theta+phi); A2 <- 0
  b1 <- b[1]; b2 <- b[2]; b3 <- b[3]
  a1 <- pc[1]; a2 <- pc[2]; a3 <- pc[3]
  n1 <- n[1]; n2 <- n[2]; n3 <- n[3]
  
  c3_1 <- (-A1*(a1*n1*n3 + a2*n2*n3 - a3*n1^2 - a3*n2^2) + (-a1*n2 + a2*n1)*sqrt(-A1^2*n1^2 - A1^2*n2^2 - A1^2*n3^2 + a1^2*n2^2 + a1^2*n3^2 - 2*a1*a2*n1*n2 - 2*a1*a3*n1*n3 + a2^2*n1^2 + a2^2*n3^2 - 2*a2*a3*n2*n3 + a3^2*n1^2 + a3^2*n2^2))/(a1^2*n2^2 + a1^2*n3^2 - 2*a1*a2*n1*n2 - 2*a1*a3*n1*n3 + a2^2*n1^2 + a2^2*n3^2 - 2*a2*a3*n2*n3 + a3^2*n1^2 + a3^2*n2^2)
  c3_2 <- (-A1*(a1*n1*n3 + a2*n2*n3 - a3*n1^2 - a3*n2^2) + (a1*n2 - a2*n1)*sqrt(-A1^2*n1^2 - A1^2*n2^2 - A1^2*n3^2 + a1^2*n2^2 + a1^2*n3^2 - 2*a1*a2*n1*n2 - 2*a1*a3*n1*n3 + a2^2*n1^2 + a2^2*n3^2 - 2*a2*a3*n2*n3 + a3^2*n1^2 + a3^2*n2^2))/(a1^2*n2^2 + a1^2*n3^2 - 2*a1*a2*n1*n2 - 2*a1*a3*n1*n3 + a2^2*n1^2 + a2^2*n3^2 - 2*a2*a3*n2*n3 + a3^2*n1^2 + a3^2*n2^2)
  
  c1_1 <-(A1*n2 + a2*c3_1*n3 - a3*c3_1*n2)/(a1*n2 - a2*n1)
  c1_2 <- (A1*n2 + a2*c3_2*n3 - a3*c3_2*n2)/(a1*n2 - a2*n1)
  
  c2_1 <- (-A1*n1 - a1*c3_1*n3 + a3*c3_1*n1)/(a1*n2 - a2*n1)
  c2_2 <-(-A1*n1 - a1*c3_2*n3 + a3*c3_2*n1)/(a1*n2 - a2*n1)
  
  c_1 <- c(c1_1, c2_1, c3_1)
  c_2 <- c(c1_2, c2_2, c3_2)
  
  dist1 <- compute_sphere_distance(b1, b2, b3, c1_1, c2_1, c3_1, 1.0)  
  dist2 <- compute_sphere_distance(b1, b2, b3, c1_2, c2_2, c3_2, 1.0)  
  
  if(dist1 < dist2){
    return(c_1)
  }else{
    return(c_2)
  }
}


cross_prod <- function (x, y)
{
  return(c(x[2]*y[3] - x[3]*y[2], -(x[1]*y[3] - x[3]*y[1]), x[1]*y[2] - x[2]*y[1]))
}


get_inner_outer_ring_rad <- function(data_ring)
{
  
  # Compute approximated ring-centre 
  pc <- c(mean(data_ring$x), mean(data_ring$y), mean(data_ring$z))
  pc <- pc / sqrt(sum(pc^2))
  
  # Find point furhest away from pc in the ring (point b)
  which_point_furthest_away <- which.max(compute_sphere_distance(pc[1], pc[2], pc[3], data_ring$x, data_ring$y, data_ring$z, 1.0)) 
  b <- c(data_ring$x[which_point_furthest_away], data_ring$y[which_point_furthest_away], data_ring$z[which_point_furthest_away])
  
  # Given the c-point we can now parameterice a circle around pc 
  c_point <- find_point_c(pc, b)
  
  theta_list = seq(0, 2*pi-0.02, by=0.02)
  v_rot <- lapply(theta_list, function(theta){
    r_point <- c_point*cos(theta) + cross_prod(pc, c_point)*sin(theta) + pc*sum(c_point*pc)*(1 - cos(theta))
    return(tibble(x=r_point[1], y = r_point[2], z = r_point[3]))}) |> 
    bind_rows()
  
  min_distance <- get_min_distace(data_ring)
  data_ring$is_inner_point <- F
  data_ring$is_outer_point <- F
  t <- seq(from = 1e-4, to = 1-1e-4, by = 1e-3)
  
  for(i in 1:nrow(v_rot)){
    
    check_inner <- T
    p_end <- as.numeric(v_rot[i, 1:3])
    for(j in 1:length(t)){
      
      t_use <- t[j]
      R_tmp <-  pc + t_use*(p_end - pc) 
      R <- R_tmp / sqrt(sum(R_tmp^2))
      dist_all_point <- compute_sphere_distance(R[1], R[2], R[3], data_ring$x, data_ring$y, data_ring$z, 1.0)
      
      if(check_inner == T && min(dist_all_point) < min_distance){
        i_min <- which.min(dist_all_point)
        data_ring$is_inner_point[i_min] <- T
        check_inner <- F
      }
      
      if(check_inner == F && min(dist_all_point) > min_distance*2.0){
        i_min <- which.min(dist_all_point)
        data_ring$is_outer_point[i_min] <- T
        break
      }
    }
  }
  
  # Compute inner and outer distances 
  inner_point <- data_ring |> filter(is_inner_point == T)
  d_inner_data <- tibble(dist = compute_sphere_distance(pc[1], pc[2], pc[3], inner_point$x, inner_point$y, inner_point$z, 1.0), 
                         label = "Inner_points")
  outer_point <- data_ring |> filter(is_outer_point == T)
  d_outer_data <- tibble(dist = compute_sphere_distance(pc[1], pc[2], pc[3], outer_point$x, outer_point$y, outer_point$z, 1.0), 
                         label = "Outer_points")
  distance_data <- bind_rows(d_inner_data, d_outer_data)
  
  return(list(dist_data=distance_data, ring_data=data_ring))
}


compute_ring_size = function(data_points, scale_dist=2.5)
{
  
  dist_outer <- data_points |> filter(is_outer_point == TRUE)
  dist_inner <- data_points |> filter(is_inner_point == TRUE)
  data_ret = tibble()
  for(i in 1:nrow(dist_inner)){
    
    max_dist_point_i = get_max_distace(dist_inner, i) 
    point_a = c(max_dist_point_i$x[1], max_dist_point_i$y[1], max_dist_point_i$z[1])
    point_b = c(dist_inner$x[i], dist_inner$y[i], dist_inner$z[i])
    c_point = find_point_c(point_a, point_b)
    
    # Find the outer point closest to c 
    dist_to_c = compute_sphere_distance(c_point[1], c_point[2], c_point[3], dist_outer$x, dist_outer$y, dist_outer$z, 1.0)
    i_min = which.min(dist_to_c)
    d_thickness = compute_sphere_distance(dist_inner$x[i], dist_inner$y[i], dist_inner$z[i], dist_outer$x[i_min], dist_outer$y[i_min], dist_outer$z[i_min], 1.0)
    
    # Rotate c_point
    p_centre = (point_a + 0.5*(point_b - point_a)) 
    p_centre = p_centre / sqrt(sum(p_centre^2))
    c_point_rot = c_point*cos(pi) + cross_prod(p_centre, c_point)*sin(pi) + p_centre*sum(c_point*p_centre)*(1 - cos(pi))
    i_min2 = which.min(compute_sphere_distance(c_point_rot[1], c_point_rot[2], c_point_rot[3], dist_outer$x, dist_outer$y, dist_outer$z, 1.0))
    d_max = compute_sphere_distance(dist_outer$x[i_min2], dist_outer$y[i_min2], dist_outer$z[i_min2], dist_outer$x[i_min], dist_outer$y[i_min], dist_outer$z[i_min], 1.0)
    
    data_ret = bind_rows(data_ret, tibble(x=dist_inner$x[i], y=dist_inner$y[i], z=dist_inner$z[i], 
                                          d_min = max_dist_point_i$max_dist*scale_dist, d_thick = d_thickness*scale_dist, 
                                          d_max = d_max*scale_dist))
  }
  
  return(data_ret)
}


process_no_cable <- function(dir_data, index_list, tag, filter_t=NULL)
{
  
  data <- read_csv(file.path(dir_data, "Distribution_data.csv"), col_types = cols())
  
  data_ret <- tibble()
  for(j in 1:length(index_list)){
    
    data_dist <- data |> filter(index == index_list[j])
    print(sprintf("Max time = %.3f", max(data_dist$t, na.rm = T)))
    if(!is.null(filter_t)){
      data_dist <- data_dist |> filter(t < filter_t)
    }
    print(sprintf("Max time post = %.3f", max(data_dist$t, na.rm = T)))
    top_time <- max(data_dist$t)
    
    if(top_time > 1000){
      next
    }
    
    data_dist_ <- data_dist |> 
      filter(t == top_time) |> 
      group_by(x, y, z) |> 
      summarise(P = mean(P))
    maxP <- max(data_dist_$P)
    data_dist_ <- data_dist_ |> filter(P > maxP * 0.2)
    res <- get_inner_outer_ring_rad(data_dist_) 
    data_tmp <- compute_ring_size(res[[2]]) |> 
      mutate(tag = tag, 
             index = index_list[j])
    data_ret <- bind_rows(data_ret, data_tmp)
  }
  return(data_ret)
}   


process_ring_size <- function(tag, dir_data1, dir_data2, case_list, index_list=c(1, 2, 3, 4, 5), filter_t=NULL)
{
  
  # Read the data 
  for(i in 1:length(case_list)){
    case <- case_list[i]
    
    dir_experiment <- file.path(dir_data1, case)
    data_dist1_full <- read_csv(file.path(dir_experiment, "Distribution_data.csv"), col_types = cols())
    dir_experiment <- file.path(dir_data2, case)
    data_dist2_full <- read_csv(file.path(dir_experiment, "Distribution_data.csv"), col_types = cols())
    
    if(!is.null(filter_t)){
      data_dist1_full <- data_dist1_full |> filter(t < filter_t)
      data_dist2_full <- data_dist2_full |> filter(t < filter_t)
    }
    
    data_sum <- tibble()
    
    for(j in 1:length(index_list)){
      
      data_dist1 <- data_dist1_full |> filter(index == index_list[j])
      top_times1 <- sort(unique(data_dist1$t), decreasing = T)[1]
      data_dist2 <- data_dist2_full |> filter(index == index_list[j])
      top_times2 <- sort(unique(data_dist2$t), decreasing = T)[1]
      
      min_top_time <- min(c(top_times1, top_times2))
      min_top_time <- 200
      top_times1 <- unique(data_dist1$t)[which.min(abs(unique(data_dist1$t) - min_top_time))]
      top_times2 <- unique(data_dist2$t)[which.min(abs(unique(data_dist2$t) - min_top_time))]
      
      data_dist1 <- data_dist1 |> 
        filter(t %in% top_times1) |> 
        group_by(x, y, z) |> 
        summarise(P = mean(P))
      maxP <- max(data_dist1$P)
      data_p1 <- data_dist1 |> filter(P > maxP * 0.5)
      res1 <- get_inner_outer_ring_rad(data_p1) 
      dist1 <- compute_ring_size(res1[[2]]) |> mutate(exocytosis = "Not_wide") 
      
      data_dist2 <- data_dist2 |> 
        filter(t %in% top_times2) |> 
        group_by(x, y, z) |> 
        summarise(P = mean(P))
      maxP <- max(data_dist2$P)
      data_p2 <- data_dist2 |> filter(P > maxP * 0.5)
      res2 <- get_inner_outer_ring_rad(data_p2) 
      dist2 <- compute_ring_size(res2[[2]]) |> mutate(exocytosis = "Wide") 
      
      data_plot = bind_rows(dist1, dist2) |> 
        mutate(d_mean = (d_min+d_max)*0.5)
      
      data_sum <- bind_rows(data_sum, data_plot |> mutate(index = index_list[j]))
      
      p1 <- ggplot(data_plot, aes(d_min, fill=exocytosis, color=exocytosis)) + 
        geom_density(alpha=0.5) + 
        labs(x = "Inner ring radius", y = "density", title = "Distribution of inner ring radius") +
        scale_color_manual(values = cbPalette[c(4, 6)], name = "Exocytosis") + 
        scale_fill_manual(values = cbPalette[c(4, 6)], name = "Exocytosis") + 
        my_theme + theme(legend.position = "bottom")
      
      p2 <- ggplot(data_plot, aes(d_max, fill=exocytosis, color=exocytosis)) + 
        geom_density(alpha=0.5) + 
        labs(x = "Outer ring radius", y = "density", title = "Distribution of inner ring radius") +
        scale_color_manual(values = cbPalette[c(4, 6)], name = "Exocytosis") + 
        scale_fill_manual(values = cbPalette[c(4, 6)], name = "Exocytosis") + 
        my_theme + theme(legend.position = "bottom")
      
      p3 <- ggplot(data_plot, aes(d_mean, fill=exocytosis, color=exocytosis)) + 
        geom_density(alpha=0.5) +
        labs(x = "(dmin+dmax)*0.5", y = "density", title = "Distribution of mean outer and inner ring radius") +
        scale_color_manual(values = cbPalette[c(4, 6)], name = "Exocytosis") + 
        scale_fill_manual(values = cbPalette[c(4, 6)], name = "Exocytosis") + 
        my_theme + theme(legend.position = "bottom")
      
      dir_save <- file.path("..", "..", "..", "Results", "Septin_ring_v1", tag, case)
      if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
      dir_save_ <- file.path(dir_save, index_list[j])
      if(!dir.exists(dir_save_)) dir.create(dir_save_)
      
      ggsave(file.path(dir_save_, "Inner.png"), p1, dpi=300, width = BASE_WIDTH, height = BASE_HEIGHT)
      ggsave(file.path(dir_save_, "Outer.png"), p2, dpi=300, width = BASE_WIDTH, height = BASE_HEIGHT)
      ggsave(file.path(dir_save_, "Mean.png"), p3, dpi=300, width = BASE_WIDTH, height = BASE_HEIGHT)
    }
    data_plot_sum <- data_sum |> 
      group_by(index, exocytosis) |> 
      summarise(median_inner = mean(d_min), 
                median_outer = mean(d_max), 
                median_mean = mean(d_mean))
    
    p1_sum <- ggplot(data_plot_sum, aes(index, median_inner)) + 
      geom_point(aes(color = exocytosis), size=4.0) + 
      scale_x_continuous(breaks = seq(min(index_list), by = 2, to = max(index_list))) + 
      scale_color_manual(values = cbPalette[-1]) + 
      my_theme
    p2_sum <- ggplot(data_plot_sum, aes(index, median_mean)) + 
      geom_point(aes(color = exocytosis), size=4.0) + 
      scale_x_continuous(breaks = seq(min(index_list), by = 2, to = max(index_list))) + 
      scale_color_manual(values = cbPalette[-1]) + 
      my_theme
    p3_sum <- ggplot(data_plot_sum, aes(index, median_outer)) + 
      geom_point(aes(color = exocytosis), size=4.0) + 
      scale_x_continuous(breaks = seq(min(index_list), by = 2, to = max(index_list))) + 
      scale_color_manual(values = cbPalette[-1]) + 
      my_theme
    ggsave(file.path(dir_save, "Sum_Inner.svg"), p1_sum, dpi=300, width = BASE_WIDTH, height = BASE_HEIGHT)
    ggsave(file.path(dir_save, "Sum_Mean.svg"), p2_sum, dpi=300, width = BASE_WIDTH, height = BASE_HEIGHT)
    ggsave(file.path(dir_save, "Sum_Outer.svg"), p3_sum, dpi=300, width = BASE_WIDTH, height = BASE_HEIGHT)
  }
  return(data_plot_sum)
}
