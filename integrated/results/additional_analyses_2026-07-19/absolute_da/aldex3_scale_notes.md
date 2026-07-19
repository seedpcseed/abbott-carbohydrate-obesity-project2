Scale model: sample.sm
s.mu = log2(gene_copies_per_ul)
s.var grid = 0.05, 0.25, 1.0 (provisional; facility qPCR technical variance unavailable)
Formula: ~ carbohydrate * timepoint + obesity_group + (1 | donor_id)
method = lme4; nsample = 128
Counts = round(genus_relative * asv_read_depth); genera with >75% zeros amalgamated to 'other'
