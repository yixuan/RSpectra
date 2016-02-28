.onLoad <- function(libname, pkgname) {
    library.dynam("rARPACK", pkgname, libname);
}

.onUnload <- function(libpath) {
    library.dynam.unload("rARPACK", libpath);
}
