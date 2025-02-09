% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcano.R
\name{makevolcano}
\alias{makevolcano}
\title{Creates Volcanos}
\usage{
makevolcano(
  data = NULL,
  DEFolder = "./data/",
  ylab = "FDR",
  ycol = "FDR",
  plotLocation = "./images/",
  level = 0,
  pCutoff = 0.05,
  FCutoff = 1,
  xpadding = 1.5,
  ypadding = 5,
  col = c("grey30", "forestgreen", "royalblue", "red2")
)
}
\arguments{
\item{data}{input data (list of \code{data.frames}s)}

\item{DEFolder}{folder for DE data if \code{df} is \code{NULL}}

\item{plotLocation}{folder to plot images}

\item{level}{level(s) to search at}

\item{pCutoff}{cutoff adjusted p for volcano}

\item{FCutoff}{cutoff log2fold change for volcano}
}
\value{
nothing
}
\description{
Creates Volcanos
}
