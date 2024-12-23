.onAttach <- function(libname, pkgname){
  packageStartupMessage(c("The multvardiv package contains the same tools as mggd, ",
                          "and similar tools for other multivariate distributions.\n\n",
                          "mggd will no longer be maintained. You may wish to install ",
                          "multvardiv instead."))
}
