# Copy this file to config/plot_settings.R and complete it with your own evidence.
# The completed local file is excluded from version control. No example study
# observations, sample decisions, or measured values are supplied here.
# Run Stage 14 from the repository root after preparing the required input files.

plot_settings <- list(
  qc = list(
    # Choose one of groupA_rep1..3 or groupB_rep1..3 for the amber highlight.
    focus_sample = NA_character_,
    # State your evidence-based action, uncertainty and any follow-up needed.
    complexity_decision = NA_character_,
    # Describe the actual featureCounts strandedness setting used for your input.
    assignment_strandedness = NA_character_,
    problem_screen = data.frame(
      problem = c("Poor-quality reads", "Low coverage", "Adaptor contamination",
                  "PCR duplicates", "GC-content bias", "Sample mislabelling"),
      status = rep("", 6),
      evidence = rep("", 6),
      decision = rep("", 6),
      # For each row choose candidate, context, or not_supported to set its colour.
      status_group = rep("", 6),
      stringsAsFactors = FALSE
    ),
    read_length = list(
      sample_id = NA_character_,  # One of the six generic sample IDs above.
      mate = NA_character_,       # R1 or R2.
      # Counts for the two read-length classes in this selected mate.
      short_reads = NA_real_,
      long_reads = NA_real_,
      short_length_bp = NA_integer_,
      long_length_bp = NA_integer_,
      # Use the actual denominator checked; do not infer complete pairing.
      paired_matches = NA_real_,
      paired_total = NA_real_,
      missing_mates = NA_real_,
      malformed_records = NA_real_,
      mismatched_mates = NA_real_,
      # Keep each panel label brief; use \n for intended line breaks.
      downstream_performance = NA_character_,
      downstream_context = NA_character_,
      decision = NA_character_,
      decision_reason = NA_character_,
      unresolved_notes = NA_character_,
      limitations = NA_character_
    )
  )
)
