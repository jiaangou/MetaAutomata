#' @title Neighborhood

#' @description Returns neighborhood coordinates given focal cell and neighborhood size

#' @param x,y coordinates of focal cell
#' @param df dataframe describing the grid, used to get the boundary coordinates
#' @param nh_size size of the neighborhood to be returned
#' @param neighborhood the type of neighborhood to use; vonNeumann by default, Moore otherwise
#' @param torus defined how the boundary condition should be treated; if TRUE, then the grid is toroidal

#' @import dplyr
#' @export

neighborhood <- function(x, y, df, nh_size = 1, neighborhood = "vonNeumann",torus = TRUE, origin = FALSE){

  #Neighborhood centered at the origin (0, 0)
  offsets <- expand.grid(x = -nh_size:nh_size, y = -nh_size:nh_size)%>% #This is equivalent to the Moore neighborhood
    filter(x != 0 | y != 0)#remove origin

  #Slice off corners to get von Neumann neighborhood
  if(neighborhood == 'vonNeumann'){
    offsets <- offsets%>%
      filter(abs(x) + abs(y) <= nh_size)
  }

  #Adjust neighborhood coordinates to center at focal x and y
  x_coords <- offsets$x + x
  y_coords <- offsets$y + y

  #Get boundary conditions from landscape df
  x_max <- max(df$x)
  y_max <- max(df$y)

  #Toroidal function - correct indices that are out of bounds
  toroidal <- function(max_i, coords){

    #indices of values that are out of bounds
    smaller_i <- which(coords<1)  # < 1
    larger_i <- which(coords>max_i)  # > max

    #correct those indices
    coords[smaller_i] <- max_i + coords[smaller_i]
    coords[larger_i] <-  coords[larger_i] - max_i

    return(coords)

  }

  #If torus, then correct indices that are out of bounds
  if(torus == TRUE){

    #Correct out of bound indices --------------
    if(any(x_coords< 1| x_coords>x_max)){
      x_coords <- toroidal(max_i = x_max, coords = x_coords)
    }

    if(any(y_coords < 1 | y_coords>y_max)){
      y_coords <- toroidal(max_i = y_max, coords = y_coords)
    }

  }

  #Get neighborhood data
  coords <- data.frame(x = x_coords, y = y_coords)

  #If origin is TRUE, put the original coordinates back in
  if(origin == TRUE){
    coords <- coords%>%
      bind_rows(data.frame(x = x, y = y))

  }
  # nh <- coords%>%
  #   left_join(df, by = c('x', 'y'))%>%
  #   relocate(ID, .before = 'x')

  return(coords)

}
