
# global options for R
options(
  htmltools.dir.version = FALSE, 
  formatR.indent = 2, width = 55, 
  digits = 2,scipen=999,tinytex.verbose = TRUE,
  knitr.kable.NA = '',
  echo=FALSE, warning=FALSE, message=FALSE,comment="")

# global options for knitr
knitr::opts_chunk$set(fig.align='center',echo = FALSE,message = FALSE,
                      warning = FALSE, comment="",
                      fig.width=14, fig.height=7,
                      out.height=450,  # do not set both height and width
                      dev = 'svglite',retina =1.5 # improve resolution
                      ) 

# global options for DT
options(htmltools.preserve.raw = FALSE)
options(DT.options = list(dom ="t" ,  # pure table with no search blank
                          columnDefs = list(
                            list(className = "dt-center", targets = "_all"), # align center
                            list(visible=FALSE,targets=0) # hide index column
                            )
                          )
        )

# global options for servr pkg

options(servr.interval = 0.5) # control time to refresh the preview
options(servr.daemon = TRUE) # unlock thread when infinite moon render

