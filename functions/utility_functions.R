# Utility functions

# load all R scripts in a specified directory
sourceDir <- function(directory) {
  script_files <- list.files(directory, pattern = "\\.R$", full.names = TRUE)
  invisible(lapply(script_files, source))
}

# create function code template
write_function <- function(fun, dest = NULL) {
  if(is.null(dest)){
    dest = getwd()
  }
  fun_name <- deparse(substitute(fun))
  params <- as.list(formals(fun))
  param_names <- names(params)
  
  # Create the function code
  code <- sprintf("#' %s\n#'\n", fun_name)
  for (param in param_names) {
    code <- paste0(code, sprintf("#' @param %s Description of argument\n", param))
  }
  code <- paste0(code, sprintf("#'\n#' @return Description of the return value\n"))
  code <- paste0(code, sprintf("#'\n#' @export\n"))
  code <- paste0(code, sprintf("%s <- %s\n", fun_name, paste(deparse(fun), collapse = "\n ")))
  
  # browser()
  # if(GPT_description){
  #   if(require("gptstudio", quietly = TRUE)) {
  #     library(gptstudio)
  #     library(gpttools)
  #     code
  #   } else {
  #     warning("'gptstudio' package is not installed. Code will remain unannotated.")
  #   }
  # }
  
  # Save the function code to a file
  file_path <- file.path(dest, paste0(fun_name, ".R"))
  
  writeLines(code, file_path)
  
  message(sprintf("Function code saved to: %s", file_path)) 
  
  # Open the code file
  file.edit(file_path)
}

## Custom install & load function
# TODO: check availability of all packages (installed, cran, github), 
# if not all installed, list their availability, and give numbered options:
# 1. load installed packages
# 2. install only available on cran
# 3. attempt to install CRAN and github 
# eg:
# These packages have more recent versions available.
# It is recommended to update all of them.
# Which would you like to update?
#   
#   1: All                          
# 2: CRAN packages only           
# 3: None                         
# 4: dplyr (1.1.2 -> 1.1.3) [CRAN]

# see example here:
# https://github.com/ropensci/rnaturalearth/blob/704772792ecafb0f29b9f68089fdab80ccf4cb43/R/install-rnaturalearthdata.R#L40
# and base R package update where it gives you multiple options including "cran only"
# consider something that searches github, perhaps based on:
# https://github.com/hoxo-m/githubinstall
# https://cloud.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html
load_lib <- function(...) {
  # convert passed arguments in ... to vector of character
  packages <- as.character(substitute(c(...)))[-1]  # Exclude the 'c'
  
  # function to load libraries
  load <- function(packages){ 
    for (pkg in packages) {
      require(pkg, character.only = TRUE, quietly = TRUE)
    }
  }
  
  # Identify installed packages
  installed <- sapply(packages, function(pkg) {
    tryCatch({
      find.package(pkg)
      TRUE
    }, error = function(e) {
      FALSE
    })
  })
  
  # Install and load missing packages
  if (sum(!installed) > 0) {
    # ask whether to install missing packages?
    cat(sprintf("%s %s %s not installed.\nAttempt to install packages from Cran? (y/n)\n", 
                ifelse(length(packages[!installed]) > 1,
                       "Packages ", "The package "),
                paste(packages[!installed], collapse = ", "), 
                ifelse(length(packages[!installed]) > 1,
                       " are ", " is ")
    ))
    choice <- readline()
    
    if (tolower(choice) == "y") {
      # Identify packages in CRAN
      cran <- packages[!installed] %in% rownames(available.packages())
      
      # Install and load missing packages
      for (pkg in packages[!installed][cran]) {
        install.packages(pkg)
        load(pkg)
      }
      
      if(sum(!cran) > 0){
        # ask whether to try to install remaining packages from github
        cat(sprintf("%s %s %s not available on CRAN.\nAttempt to install from GitHub using githubinstall? (y/n) \n", 
                    ifelse(length(packages[!installed][!cran]) > 1,
                           "Packages ", "The package "),
                    paste(packages[!installed][!cran], collapse = ", "), 
                    ifelse(length(packages[!installed][!cran]) > 1,
                           " are ", " is ")
        ))
        choice <- readline()
        
        # attempt to install from github:
        if (tolower(choice) == "y") {
          githubinstall::gh_install_packages(packages[!installed][!cran])
        }
      }
    } else {
      cat("Loading installed packages.\n")
      load(packages[installed])
      warning(sprintf("The following package%s remain%s uninstalled:\n%s",
                      ifelse(length(packages[!installed]) > 1,"s",""),
                      ifelse(length(packages[!installed]) > 1,"s",""),
                      paste(packages[!installed], collapse = ", ")))
    }
  } else {
    load(packages)
  }
}

# # install/load packages
# sapply(packages, function(.x){
#   if (!require(.x, character.only = TRUE, quietly = TRUE)) {
#     install.packages(.x, repos = "http://cran.us.r-project.org")
#   }
#   require(.x, character.only = TRUE)
# })
# }

# not in (opposite of %in%)
'%!in%' <- function(x,y)!('%in%'(x,y))

# save generic plot
sv <- function(x, filename = "plot.png", path = ".", width = 500, height = 500, ...) {
  # capture x argument call
  expr <- substitute(x)
  # ensure path exists
  make_path(path)
  # combine file name
  filename <- file.path(path, filename)
  # save plot
  png(filename, width, height, ...)
  if(isRplot(x, expr)) { 
    out <- capture.output(print(x))
  } else if (inherits(x, "ggplot")) {
    print(x)  # for ggplot objects
  } else {
    # Attempt to use generic plot function for other types
    plot(x, ...)
  }

  dev.off()
  # # # return console output if not just "NULL"
  # if(out[[1]] != "NULL") cat("Plot ruturned:\n", out, "\n")
  # print path to file (without remote home directory)
  cat("Saved file to:\n", gsub(gsub("\\$HOME", "", system2("echo", args = "$HOME", stdout = TRUE)), "", 
                              getwd()), 
      gsub("^\\.", "", filename), "\n", sep = "")
  # View by running browseURL("Z:/R/beetles./plot.png")
}

isRplot <- function(x, expr = NULL) {
  # Allow passing an already substituted expression for nested calls
  if (is.null(expr)) {
    expr <- substitute(x)
  }
  
  # Recursive helper function to search for plot or hist calls
  containsPlotCall <- function(exp) {
    if (is.call(exp)) {
      if (as.character(exp[[1]]) %in% c("plot", "hist")) {
        return(TRUE)
      }
      for (i in seq_along(exp)) {
        if (containsPlotCall(exp[[i]])) {
          return(TRUE)
        }
      }
    } else if (is.list(exp) || is.expression(exp)) {
      for (item in exp) {
        if (containsPlotCall(item)) {
          return(TRUE)
        }
      }
    }
    return(FALSE)
  }
  
  containsPlotCall(expr)
}

# make missing directory
make_path <- function(...){
  path <- file.path(...)
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# update display number
update_display <- function(DISPLAY = NULL) {
  if (is.null(DISPLAY)) {
    DISPLAY <- Sys.getenv("DISPLAY")
  }
  
  if (DISPLAY == "") stop("No DISPLAY number defined")
  
  # Extract the numeric part of the DISPLAY string
  display_number <- as.numeric(sub("^.*:([0-9]+)\\..*$", "\\1", DISPLAY))
  
  if (is.na(display_number)) stop("Invalid DISPLAY number")
  
  # Update the DISPLAY environment variable
  Sys.setenv(DISPLAY = sprintf("localhost:%.1f", display_number))
  
  # Evaluate SSH and DISPLAY environment variables from tmux
  system("eval $(tmux showenv -s | grep -E '^(SSH|DISPLAY)')")
  
  # Echo the DISPLAY value
  system("echo $DISPLAY")
  
  # Try to close any existing graphics devices
  tryCatch({
    dev.off()
  }, error = function(e) {})
  
  # Open a new X11 device
  x11()
  plot.new()
  text(sprintf("Display number set to: %.1f", display_number), x = 0.5, y = 0.5)
}

# combine vector of character strings into proper grammar list
combine_strings <- function(strings, collapse = ", ", final = "and") {
  n <- length(strings)
  
  if (n == 1) {
    return(strings)
  } else if (n == 2) {
    return(paste(strings, collapse = paste0(" ", final, " ")))
  } else {
    return(paste0(paste(strings[1:(n-1)], collapse = collapse), 
                  collapse, final, " ", strings[n]))
  }
}
