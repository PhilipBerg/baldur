% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compile_model.R
\name{compile_model}
\alias{compile_model}
\title{Wrapper for compiling Stan model}
\usage{
compile_model()
}
\value{
A complied Stan model
}
\description{
As of now, Rstan models complied with packages cannot sample when ran in a multidplyr backend.
As a temporary solution, this function compiles the Stan model and can be used to run Baldur on several cores.
}
