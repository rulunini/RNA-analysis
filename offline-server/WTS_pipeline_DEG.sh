#!/bin/bash

Rscript ./GEO_analysis.r -i ./ -c sample_counts.csv -m sample_meta.csv -o ./geo_output
Rscript ./2_batch_effect.r -i ./geo_output -o ./batch_output
Rscript ./3_gsva.r -i ./batch_output -c sample_counts.csv -m sample_meta.csv -o ./gsva_output
Rscript ./4_deg_analysis.r -i ./gsva_output -o ./deg_output
