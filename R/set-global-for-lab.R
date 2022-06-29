
# global options for R
options(
  htmltools.dir.version = FALSE, 
  formatR.indent = 2, width = 55, 
  digits = 2,scipen=999,tinytex.verbose = TRUE,
  knitr.kable.NA = '',
  echo=FALSE, warning=FALSE, message=FALSE,comment="")

# global options for knitr
knitr::opts_chunk$set(fig.align='center',
                      echo = TRUE,message = FALSE,
                      warning = FALSE, 
                      comment="",
                      attr.source='.numberLines',
                      class.source = "foldable",
                      fig.width=17, fig.height=7,
                      out.height = 600,
                      dev = 'svglite',fig.retina = 1 # improve resolution
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
