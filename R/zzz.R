.onLoad <- function(libname, pkgname) {
    library.dynam("RSpectra", pkgname, libname);
}

.onUnload <- function(libpath) {
    library.dynam.unload("RSpectra", libpath);
}
