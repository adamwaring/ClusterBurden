## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7
)

## -----------------------------------------------------------------------------
library(ClusterBurden)

x = annotate_with_gnomad_freqs(dataset = data.table(symbol="MYH7", chr="14", pos=23882985, ref="G", alt="A"))
print(x)

x = annotate_with_gnomad_freqs(dataset = data.table(symbol="MYH7", txconsq="c.5704G>C"))
print(x)

## -----------------------------------------------------------------------------
get_residue_position("p.R502W")

x = data.table(pconsq = c("p.Leu1903Leu", "p.Lys1173_Ala1174delinsThr"))
x[,protein_position := get_residue_position(pconsq)]
x

## -----------------------------------------------------------------------------
x = data.table(ref = c("AGGATGG", "G"), alt = c("A", "GCACACA"))

x[,n_res:=mapply(function(x, y) max(nchar(x)%/%3, nchar(y)%/%3), ref, alt)]
x

## -----------------------------------------------------------------------------
x = data.table(chr = "14", pos = 23882071:23882081, over_10=1)

format_coverage(x)

