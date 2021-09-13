# A list of todo items off the top of our heads

- vignette, fully fledged to show off the nice method :)
- roxygen with markdown, for ease of styling?

- description fields to be proper filled
- biocviews
- license
- namespace imports and so - handled at best in centralized place?
- readme Rmd? this would help if we have "some code dependent" content

- CI integration - GH actions, quite easy to setup and very handy!

- unit testing 

- data -> enough with the small ones, or need EHub object?

- harmonize filenames/function names

- Working off R CMD check + BiocCheck

- cleanup of notes from gg-based calls?

- dedicated examples

- pkgdown?
- code of conduct?
- NEWS.md
- message vs print

- some refs for the manuscript and to manuscripts

- have verbose parameter to control the amount of messages printed out?

#
# Further notes from brainstorming
#

- return() in bemkl_supervised...R
- try() can be removed (e.g. in compare_immune_response/ predict_immune_response.R)
- consistent variable/attribute style 
- verify if attributes are really used in functions (see compute_CCpair.R --> something like "compute_cell_fraction()")
- validate the input objects (create checks ...)
- compute_gold_standard.R -->  more general form (dynamic function)
- compute.yers_expIS etc. are all the same --> use intern function e.g. "compute_signature_generics" 
--> outsource the simmilar functionallity into a separate file/function for the genes as "meta-level" (e.g. compute.IMPRES.R)
- compute.TIDE.R etc. --> find more robust way for calling python script (+ change paste0() to file.path())
- change F to FALSE and T to TRUE
- ipsmap.R --> parameterize hard coded threshold
- predict_with_bemkl.R & predict_with_multi...R --> K to parameter/ probably optimize the for loops
- general: consistent spacing and style in all functions
- predict_with_multi... --> replace do.call(c) with unlist()


#
# Arsenij's notes:
#
- compute_cell_fraction.R l.41 old_cellnames are not used further
- make one function for all remove.genes.ICB_proxies

- error occurs, when trying to run compute_ICB_genes: Error in "rownames<-". 
