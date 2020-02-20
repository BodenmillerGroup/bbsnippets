MultiLayerColourVector <- function(layer1, layer2 = NULL, layer3 = NULL,
                                   col_layer1 = NULL, force_colours = FALSE){
  # Check inputs
  if(!is.character(layer1)){
    stop("'layer1' not a character vector.")
  }
  
  if(!is.null(layer2)){
    if(!is.list(layer2)){
      stop("'layer2' not a list.")
    }
    if(is.null(names(layer2))){
      stop("Please supply a named list as 'layer2'")
    }
    if(!all(unlist(lapply(layer2, is.character)))){
      stop("'layer2' contains other than character entries.")
    }
    if(!all(names(layer2) %in% layer1)){
      stop("Some/all names of 'layer2' cannot be found in 'layer1'")
    }
  }
  
  if(!is.null(layer3)){
    if(!is.list(layer3)){
      stop("'layer2' not a list.")
    }
    if(is.null(names(layer3))){
      stop("Please supply a named list as 'layer3'")
    }
    if(!all(unlist(lapply(layer3, is.character)))){
      stop("'layer3' contains other than character entries.")
    }
    if(!all(names(layer3) %in% as.character(unlist(layer2)))){
      stop("Some/all names of 'layer3' cannot be found in 'layer2'")
    }
  }
  
  if(!is.null(col_layer1)){
    if(!is.character(col_layer1) || is.null(names(col_layer1))){
      stop("Colours for layer 1 need to be a named character vector")
    }
    if(!all(names(col_layer1) %in% layer1)){
      stop("Some of the classes specified in 'col_layer1' are not in 'layer1'")
    }
    if(!all(as.character(col_layer1) %in% c("green", "blue", "purple",
                                            "red", "orange", "yellow", "brown"))){
      stop("Only the following colours are supported:\n", 
           "green, blue, purple, red, orange, yellow, brown")
    }
  }
  
  if(!is.logical(force_colours)){
    stop("'force_colours' needs to be a logical.")
  }
  
  # I will store the colour vectors in a list
  out <- list(layer1 = NULL, layer2 = NULL, layer3 = NULL)
  
  # Select colours for main classes
  # Here, I select 'nicer' colours for the main classes
  # Standard colour that are supplied by the user are replaced
  # Colours: green, yellow, blue, grey, purple, red, orange
  if(length(layer1) > 7){
    stop("Only 7 classes allowed in first layer")
  }
  
  cur_colours <- c("#33A02C", "#1F78B4", "#6A3D9A", 
                   "#E31A1C", "#FF7F00", "yellow2", "#B15928")
  names(cur_colours) <- c("green", "blue", "purple",
                          "red", "orange", "yellow", "brown")
  
  if(is.null(col_layer1)){
    layer1_vec <- cur_colours[1:length(layer1)]
    names(layer1_vec) <- layer1
  } else if(!force_colours){
    layer1_vec <- cur_colours[col_layer1]
    names(layer1_vec) <-  names(col_layer1)
  } else{
    layer1_vec <- col_layer1
  }
  
  # Layer 2
  # Green
  green_col <- c("springgreen1", "springgreen3", "greenyellow", 
                 "palegreen1", "palegreen3", 
                 "green3", "green4", "darkgreen")
  
  # Blue
  blue_col <- c("skyblue", "lightblue", "turquoise", "turquoise3", 
                "deepskyblue", "blue1", "deepskyblue4", "blue4")
  
  # Purple
  purple_col <- c("plum1", "plum3", "violet", "magenta3", 
                  "darkorchid1", "mediumorchid3", "blueviolet", "darkorchid4")
  
  # Red
  red_col <- c("indianred1", "lightcoral", "pink2", "indianred3", "red1", 
                  "firebrick2", "red3", "red4")
  
  # Orange
  orange_col <- c("lightsalmon", "sienna1", "salmon", "tan1", "darkorange", 
                  "tan3", "orange2", "tan4")
  
  # Yellow
  yellow_col <- c("yellow1", "lightgoldenrod1", "khaki", "gold1", 
                  "goldenrod2", "gold3", "yellow3", "goldenrod4")
  
  # Brown
  brown_col <- c("coral","coral2", "tomato2",
                  "peru", "sienna3", "brown3", "sienna", "brown4")
  
  
}