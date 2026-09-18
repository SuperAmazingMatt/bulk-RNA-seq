#!/usr/bin/env Rscript

# Check paired FASTQ integrity, summarise mate-length combinations, and plot QC.
#
# This script reads gzip-compressed FASTQ files in chunks. It never modifies the
# input files and does not hold an entire FASTQ file in memory.
#
# Example environment-module setup:
#   module load R/4.3.3
#
# Inspect the files that would be analysed (no FASTQ scan and no output files):
#   Rscript scripts/check_and_plot_pair_lengths.R --stage both --samples all --dry-run
#
# Run the full raw + trimmed-paired check:
#   Rscript scripts/check_and_plot_pair_lengths.R --stage both --samples all

options(stringsAsFactors = FALSE)

abort <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

usage <- function(status = 0L) {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript check_and_plot_pair_lengths.R --stage raw|trimmed|both [options]\n\n",
      "Required:\n",
      "  --stage VALUE          raw, trimmed, or both\n\n",
      "Options:\n",
      "  --samples VALUE        all (default) or comma-separated sample names\n",
      "  --chunk-size N         FASTQ records per chunk (default: 50000)\n",
      "  --project-root PATH    Project root; normally inferred from scripts/\n",
      "  --output-dir PATH      Results directory; normally created automatically\n",
      "  --export-failures      Save only true pairing failures as .fastq.gz\n",
      "  --dry-run              Show matched input files without scanning them\n",
      "  --help                 Show this help\n\n",
      "Definitions:\n",
      "  Unequal R1/R2 lengths are valid QC observations, not pairing failures.\n",
      "  A true failure is a malformed record, mismatched primary ID, wrong mate\n",
      "  tag, missing mate, or truncated FASTQ record.\n"
    )
  )
  quit(save = "no", status = status)
}

parse_cli <- function(args) {
  out <- list(
    stage = NULL,
    samples = "all",
    chunk_size = 50000L,
    project_root = NULL,
    output_dir = NULL,
    export_failures = FALSE,
    dry_run = FALSE
  )

  value_options <- c(
    "--stage", "--samples", "--chunk-size", "--project-root", "--output-dir"
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg == "--help") {
      usage(0L)
    } else if (arg == "--export-failures") {
      out$export_failures <- TRUE
      i <- i + 1L
      next
    } else if (arg == "--dry-run") {
      out$dry_run <- TRUE
      i <- i + 1L
      next
    }

    if (grepl("^--[^=]+=", arg)) {
      key <- sub("=.*$", "", arg)
      value <- sub("^[^=]+=", "", arg)
    } else {
      key <- arg
      if (!(key %in% value_options)) {
        abort("Unknown option: %s", key)
      }
      if (i == length(args)) {
        abort("Option %s requires a value", key)
      }
      value <- args[[i + 1L]]
      i <- i + 1L
    }

    if (!(key %in% value_options)) {
      abort("Unknown option: %s", key)
    }

    if (key == "--stage") out$stage <- tolower(value)
    if (key == "--samples") out$samples <- value
    if (key == "--chunk-size") {
      parsed <- suppressWarnings(as.integer(value))
      if (is.na(parsed) || parsed < 1L) {
        abort("--chunk-size must be a positive integer")
      }
      out$chunk_size <- parsed
    }
    if (key == "--project-root") out$project_root <- value
    if (key == "--output-dir") out$output_dir <- value

    i <- i + 1L
  }

  if (is.null(out$stage)) {
    abort("--stage is required (raw, trimmed, or both)")
  }
  if (!(out$stage %in% c("raw", "trimmed", "both"))) {
    abort("Invalid --stage value: %s", out$stage)
  }
  out
}

script_path <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- full_args[grepl("^--file=", full_args)]
  if (!length(file_arg)) return(NA_character_)
  sub("^--file=", "", file_arg[[1L]])
}

infer_project_root <- function(user_root = NULL) {
  if (!is.null(user_root)) {
    root <- normalizePath(path.expand(user_root), mustWork = TRUE)
  } else {
    this_script <- script_path()
    if (is.na(this_script)) {
      root <- normalizePath(getwd(), mustWork = TRUE)
    } else {
      script_dir <- dirname(normalizePath(this_script, mustWork = TRUE))
      root <- if (basename(script_dir) == "scripts") dirname(script_dir) else getwd()
      root <- normalizePath(root, mustWork = TRUE)
    }
  }

  if (!dir.exists(file.path(root, "data"))) {
    abort("Raw-data directory not found: %s", file.path(root, "data"))
  }
  if (!dir.exists(file.path(root, "results"))) {
    abort("Results directory not found: %s", file.path(root, "results"))
  }
  root
}

resolve_output_path <- function(path, project_root) {
  if (is.null(path)) {
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    return(file.path(
      project_root, "results", "pairing_check",
      sprintf("R_%s_%s", stamp, Sys.getpid())
    ))
  }

  expanded <- path.expand(path)
  is_absolute <- grepl("^/", expanded) || grepl("^[A-Za-z]:[/\\\\]", expanded)
  if (!is_absolute) expanded <- file.path(project_root, expanded)
  normalizePath(expanded, mustWork = FALSE)
}

discover_stage <- function(stage, project_root) {
  if (stage == "raw") {
    input_dir <- file.path(project_root, "data")
    r1_pattern <- "_L[0-9]+_R1_001\\.fastq\\.gz$"
    r1_files <- list.files(input_dir, pattern = r1_pattern, full.names = TRUE)
    r2_files <- sub("_R1_001\\.fastq\\.gz$", "_R2_001.fastq.gz", r1_files)
    samples <- sub(r1_pattern, "", basename(r1_files))
  } else if (stage == "trimmed") {
    input_dir <- file.path(project_root, "results", "02_trimmomatic_results")
    if (!dir.exists(input_dir)) {
      abort("Trimmed-data directory not found: %s", input_dir)
    }
    r1_pattern <- "_R1\\.trim_pe\\.fastq\\.gz$"
    r1_files <- list.files(input_dir, pattern = r1_pattern, full.names = TRUE)
    r2_files <- sub("_R1\\.trim_pe\\.fastq\\.gz$", "_R2.trim_pe.fastq.gz", r1_files)
    samples <- sub(r1_pattern, "", basename(r1_files))
  } else {
    abort("Internal error: unsupported stage %s", stage)
  }

  if (!length(r1_files)) {
    abort("No %s R1 files found in %s", stage, input_dir)
  }
  missing_r2 <- !file.exists(r2_files)
  if (any(missing_r2)) {
    abort(
      "Missing %s R2 counterpart(s): %s",
      stage,
      paste(basename(r2_files[missing_r2]), collapse = ", ")
    )
  }
  if (anyDuplicated(samples)) {
    abort("Duplicate normalized sample names detected for the %s stage", stage)
  }

  data.frame(
    stage = stage,
    sample = samples,
    r1_file = normalizePath(r1_files, mustWork = TRUE),
    r2_file = normalizePath(r2_files, mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

discover_pairs <- function(stage_option, sample_option, project_root) {
  stages <- if (stage_option == "both") c("raw", "trimmed") else stage_option
  pairs <- do.call(rbind, lapply(stages, discover_stage, project_root = project_root))

  if (tolower(sample_option) != "all") {
    wanted <- trimws(strsplit(sample_option, ",", fixed = TRUE)[[1L]])
    wanted <- wanted[nzchar(wanted)]
    if (!length(wanted)) abort("--samples did not contain a sample name")

    for (stage in stages) {
      present <- pairs$sample[pairs$stage == stage]
      missing <- setdiff(wanted, present)
      if (length(missing)) {
        abort(
          "Sample(s) not found at %s stage: %s",
          stage,
          paste(missing, collapse = ", ")
        )
      }
    }
    pairs <- pairs[pairs$sample %in% wanted, , drop = FALSE]
  }

  pairs <- pairs[order(match(pairs$stage, c("raw", "trimmed")), pairs$sample), ]
  rownames(pairs) <- NULL
  pairs
}

normalize_primary_id <- function(header) {
  first_token <- sub("[[:space:]].*$", "", sub("^@", "", header))
  sub("/[12]$", "", first_token)
}

extract_mate_tag <- function(header) {
  clean <- sub("^@", "", header)
  first_token <- sub("[[:space:]].*$", "", clean)
  has_second <- grepl("[[:space:]]", clean)
  second_token <- rep("", length(clean))
  second_token[has_second] <- sub(
    "[[:space:]].*$", "",
    sub("^[^[:space:]]+[[:space:]]+", "", clean[has_second])
  )

  tag <- rep(NA_integer_, length(clean))

  first_has_tag <- grepl("/[12]$", first_token)
  tag[first_has_tag] <- as.integer(sub("^.*/([12])$", "\\1", first_token[first_has_tag]))

  pending <- is.na(tag)
  second_is_casava <- pending & grepl("^[12]:", second_token)
  tag[second_is_casava] <- as.integer(substr(second_token[second_is_casava], 1L, 1L))

  pending <- is.na(tag)
  second_has_slash_tag <- pending & grepl("/[12]$", second_token)
  tag[second_has_slash_tag] <- as.integer(
    sub("^.*/([12])$", "\\1", second_token[second_has_slash_tag])
  )
  tag
}

parse_fastq_block <- function(lines) {
  n_records <- length(lines) %/% 4L
  trailing_lines <- length(lines) %% 4L

  if (!n_records) {
    return(list(
      lines = lines,
      n = 0L,
      trailing_lines = trailing_lines,
      header = character(),
      sequence = character(),
      plus = character(),
      quality = character(),
      primary_id = character(),
      mate_tag = integer(),
      length = integer(),
      valid = logical()
    ))
  }

  starts <- seq.int(1L, by = 4L, length.out = n_records)
  header <- lines[starts]
  sequence <- lines[starts + 1L]
  plus <- lines[starts + 2L]
  quality <- lines[starts + 3L]
  primary_id <- normalize_primary_id(header)
  sequence_length <- nchar(sequence, type = "bytes")
  quality_length <- nchar(quality, type = "bytes")

  valid <- startsWith(header, "@") &
    nzchar(primary_id) &
    startsWith(plus, "+") &
    sequence_length > 0L &
    sequence_length == quality_length

  list(
    lines = lines,
    n = n_records,
    trailing_lines = trailing_lines,
    header = header,
    sequence = sequence,
    plus = plus,
    quality = quality,
    primary_id = primary_id,
    mate_tag = extract_mate_tag(header),
    length = sequence_length,
    valid = valid
  )
}

record_lines <- function(parsed, indices) {
  if (!length(indices)) return(character())
  unlist(lapply(indices, function(i) {
    start <- 4L * (i - 1L) + 1L
    parsed$lines[start:(start + 3L)]
  }), use.names = FALSE)
}

add_reason <- function(reason, condition, label) {
  indices <- which(condition)
  if (!length(indices)) return(reason)
  reason[indices] <- ifelse(
    nzchar(reason[indices]),
    paste0(reason[indices], ";", label),
    label
  )
  reason
}

add_length_counts <- function(count_env, r1_length, r2_length) {
  if (!length(r1_length)) return(invisible(NULL))
  tab <- table(paste(r1_length, r2_length, sep = ","))
  for (key in names(tab)) {
    old <- if (exists(key, envir = count_env, inherits = FALSE)) {
      get(key, envir = count_env, inherits = FALSE)
    } else {
      0
    }
    assign(key, old + as.numeric(tab[[key]]), envir = count_env)
  }
  invisible(NULL)
}

length_counts_frame <- function(count_env, stage, sample, denominator) {
  keys <- ls(envir = count_env, all.names = TRUE)
  if (!length(keys)) {
    return(data.frame(
      stage = character(), sample = character(),
      r1_length = integer(), r2_length = integer(),
      length_difference = integer(), absolute_difference = integer(),
      pair_count = numeric(), percent_of_valid_matched_pairs = numeric()
    ))
  }

  split_keys <- strsplit(keys, ",", fixed = TRUE)
  r1_length <- as.integer(vapply(split_keys, `[[`, character(1L), 1L))
  r2_length <- as.integer(vapply(split_keys, `[[`, character(1L), 2L))
  counts <- as.numeric(unlist(mget(keys, envir = count_env), use.names = FALSE))
  difference <- r1_length - r2_length

  data.frame(
    stage = stage,
    sample = sample,
    r1_length = r1_length,
    r2_length = r2_length,
    length_difference = difference,
    absolute_difference = abs(difference),
    pair_count = counts,
    percent_of_valid_matched_pairs = if (denominator > 0) {
      100 * counts / denominator
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

scan_pair <- function(
    stage, sample, r1_file, r2_file, chunk_size,
    output_dir, export_failures = FALSE, preview_limit = 1000L) {

  message(sprintf("[%s / %s] scanning paired FASTQ files", stage, sample))
  started <- Sys.time()

  r1_connection <- gzfile(r1_file, open = "rt")
  r2_connection <- gzfile(r2_file, open = "rt")
  on.exit(try(close(r1_connection), silent = TRUE), add = TRUE)
  on.exit(try(close(r2_connection), silent = TRUE), add = TRUE)

  length_env <- new.env(hash = TRUE, parent = emptyenv())
  preview <- list()
  preview_n <- 0L

  safe_stem <- gsub("[^A-Za-z0-9_.-]", "_", paste(stage, sample, sep = "__"))
  failure_dir <- file.path(output_dir, "failed_reads")
  failure_manifest <- NULL
  failed_r1_connection <- NULL
  failed_r2_connection <- NULL
  exported_r1 <- 0
  exported_r2 <- 0

  if (export_failures) {
    dir.create(failure_dir, recursive = TRUE, showWarnings = FALSE)
    failure_manifest <- file(
      file.path(failure_dir, paste0(safe_stem, ".failure_manifest.tsv")),
      open = "wt"
    )
    writeLines(
      "stage\tsample\tposition\treason\tr1_header\tr2_header",
      failure_manifest
    )
    on.exit(try(close(failure_manifest), silent = TRUE), add = TRUE)

    failed_r1_connection <- gzfile(
      file.path(failure_dir, paste0(safe_stem, ".failed_R1.fastq.gz")),
      open = "wt"
    )
    failed_r2_connection <- gzfile(
      file.path(failure_dir, paste0(safe_stem, ".failed_R2.fastq.gz")),
      open = "wt"
    )
    on.exit(try(close(failed_r1_connection), silent = TRUE), add = TRUE)
    on.exit(try(close(failed_r2_connection), silent = TRUE), add = TRUE)
  }

  open_failed_r1 <- function() {
    if (is.null(failed_r1_connection)) abort("Internal error: R1 failure output is closed")
    failed_r1_connection
  }

  open_failed_r2 <- function() {
    if (is.null(failed_r2_connection)) abort("Internal error: R2 failure output is closed")
    failed_r2_connection
  }

  add_preview <- function(position, reason, r1_header, r2_header) {
    room <- preview_limit - preview_n
    if (room <= 0L || !length(position)) return(invisible(NULL))
    take <- seq_len(min(room, length(position)))
    preview[[length(preview) + 1L]] <<- data.frame(
      stage = stage,
      sample = sample,
      position = position[take],
      reason = reason[take],
      r1_header = gsub("[\\t\\r\\n]", " ", r1_header[take]),
      r2_header = gsub("[\\t\\r\\n]", " ", r2_header[take]),
      stringsAsFactors = FALSE
    )
    preview_n <<- preview_n + length(take)
    invisible(NULL)
  }

  write_manifest <- function(position, reason, r1_header, r2_header) {
    if (!export_failures || !length(position)) return(invisible(NULL))
    clean <- function(x) gsub("[\\t\\r\\n]", " ", x)
    rows <- paste(
      stage, sample, position, clean(reason), clean(r1_header), clean(r2_header),
      sep = "\t"
    )
    writeLines(rows, failure_manifest)
    invisible(NULL)
  }

  totals <- list(
    r1_records = 0, r2_records = 0, compared_positions = 0,
    valid_both_records = 0, valid_matched_pairs = 0,
    matching_primary_ids = 0, id_mismatches = 0,
    malformed_r1 = 0, malformed_r2 = 0,
    missing_r1 = 0, missing_r2 = 0,
    trailing_lines_r1 = 0, trailing_lines_r2 = 0,
    mate_tags_detected_r1 = 0, mate_tags_detected_r2 = 0,
    wrong_mate_tags_r1 = 0, wrong_mate_tags_r2 = 0,
    equal_length_pairs = 0, unequal_length_pairs = 0
  )

  position_offset <- 0
  next_progress <- 5000000

  repeat {
    r1_lines <- readLines(r1_connection, n = 4L * chunk_size, warn = FALSE)
    r2_lines <- readLines(r2_connection, n = 4L * chunk_size, warn = FALSE)
    if (!length(r1_lines) && !length(r2_lines)) break

    p1 <- parse_fastq_block(r1_lines)
    p2 <- parse_fastq_block(r2_lines)
    n_common <- min(p1$n, p2$n)
    n_positions <- max(p1$n, p2$n)
    positions <- position_offset + seq_len(n_positions)

    totals$r1_records <- totals$r1_records + p1$n
    totals$r2_records <- totals$r2_records + p2$n
    totals$compared_positions <- totals$compared_positions + n_common
    totals$malformed_r1 <- totals$malformed_r1 + sum(!p1$valid)
    totals$malformed_r2 <- totals$malformed_r2 + sum(!p2$valid)
    totals$trailing_lines_r1 <- totals$trailing_lines_r1 + p1$trailing_lines
    totals$trailing_lines_r2 <- totals$trailing_lines_r2 + p2$trailing_lines
    totals$mate_tags_detected_r1 <- totals$mate_tags_detected_r1 + sum(!is.na(p1$mate_tag))
    totals$mate_tags_detected_r2 <- totals$mate_tags_detected_r2 + sum(!is.na(p2$mate_tag))

    wrong_all_r1 <- !is.na(p1$mate_tag) & p1$mate_tag != 1L
    wrong_all_r2 <- !is.na(p2$mate_tag) & p2$mate_tag != 2L
    totals$wrong_mate_tags_r1 <- totals$wrong_mate_tags_r1 + sum(wrong_all_r1)
    totals$wrong_mate_tags_r2 <- totals$wrong_mate_tags_r2 + sum(wrong_all_r2)

    if (n_common > 0L) {
      common <- seq_len(n_common)
      valid1 <- p1$valid[common]
      valid2 <- p2$valid[common]
      valid_both <- valid1 & valid2
      id_match <- p1$primary_id[common] == p2$primary_id[common]
      wrong1 <- wrong_all_r1[common]
      wrong2 <- wrong_all_r2[common]
      good <- valid_both & id_match & !wrong1 & !wrong2

      totals$valid_both_records <- totals$valid_both_records + sum(valid_both)
      totals$matching_primary_ids <- totals$matching_primary_ids + sum(valid_both & id_match)
      totals$id_mismatches <- totals$id_mismatches + sum(valid_both & !id_match)
      totals$valid_matched_pairs <- totals$valid_matched_pairs + sum(good)

      if (any(good)) {
        good_r1_length <- p1$length[common][good]
        good_r2_length <- p2$length[common][good]
        same_length <- good_r1_length == good_r2_length
        totals$equal_length_pairs <- totals$equal_length_pairs + sum(same_length)
        totals$unequal_length_pairs <- totals$unequal_length_pairs + sum(!same_length)
        add_length_counts(length_env, good_r1_length, good_r2_length)
      }

      reason <- rep("", n_common)
      reason <- add_reason(reason, !valid1, "malformed_R1")
      reason <- add_reason(reason, !valid2, "malformed_R2")
      reason <- add_reason(reason, valid_both & !id_match, "primary_ID_mismatch")
      reason <- add_reason(reason, wrong1, "wrong_R1_mate_tag")
      reason <- add_reason(reason, wrong2, "wrong_R2_mate_tag")
      failed <- nzchar(reason)

      if (any(failed)) {
        failed_indices <- which(failed)
        add_preview(
          positions[failed_indices], reason[failed_indices],
          p1$header[failed_indices], p2$header[failed_indices]
        )
        write_manifest(
          positions[failed_indices], reason[failed_indices],
          p1$header[failed_indices], p2$header[failed_indices]
        )

        exportable <- which(valid_both & (!id_match | wrong1 | wrong2))
        if (export_failures && length(exportable)) {
          writeLines(record_lines(p1, exportable), open_failed_r1())
          writeLines(record_lines(p2, exportable), open_failed_r2())
          exported_r1 <- exported_r1 + length(exportable)
          exported_r2 <- exported_r2 + length(exportable)
        }
      }
    }

    if (p1$n > n_common) {
      extra <- (n_common + 1L):p1$n
      totals$missing_r2 <- totals$missing_r2 + length(extra)
      extra_reason <- rep("missing_R2", length(extra))
      extra_headers <- p1$header[extra]
      add_preview(
        positions[extra], extra_reason, extra_headers,
        rep("", length(extra))
      )
      write_manifest(
        positions[extra], extra_reason, extra_headers,
        rep("", length(extra))
      )
      exportable <- extra[p1$valid[extra]]
      if (export_failures && length(exportable)) {
        writeLines(record_lines(p1, exportable), open_failed_r1())
        exported_r1 <- exported_r1 + length(exportable)
      }
    }

    if (p2$n > n_common) {
      extra <- (n_common + 1L):p2$n
      totals$missing_r1 <- totals$missing_r1 + length(extra)
      extra_reason <- rep("missing_R1", length(extra))
      extra_headers <- p2$header[extra]
      add_preview(
        positions[extra], extra_reason,
        rep("", length(extra)), extra_headers
      )
      write_manifest(
        positions[extra], extra_reason,
        rep("", length(extra)), extra_headers
      )
      exportable <- extra[p2$valid[extra]]
      if (export_failures && length(exportable)) {
        writeLines(record_lines(p2, exportable), open_failed_r2())
        exported_r2 <- exported_r2 + length(exportable)
      }
    }

    position_offset <- position_offset + n_positions
    if (position_offset >= next_progress) {
      message(sprintf(
        "[%s / %s] %s positions checked",
        stage, sample,
        format(position_offset, big.mark = ",", scientific = FALSE, trim = TRUE)
      ))
      while (next_progress <= position_offset) next_progress <- next_progress + 5000000
    }
  }

  true_failure <- totals$id_mismatches > 0 ||
    totals$malformed_r1 > 0 || totals$malformed_r2 > 0 ||
    totals$missing_r1 > 0 || totals$missing_r2 > 0 ||
    totals$trailing_lines_r1 > 0 || totals$trailing_lines_r2 > 0 ||
    totals$wrong_mate_tags_r1 > 0 || totals$wrong_mate_tags_r2 > 0

  status <- if (true_failure) {
    "FAIL"
  } else if (totals$unequal_length_pairs > 0) {
    "PASS_WITH_LENGTH_DIFFERENCES"
  } else {
    "PASS"
  }

  elapsed_minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  unequal_percent <- if (totals$valid_matched_pairs > 0) {
    100 * totals$unequal_length_pairs / totals$valid_matched_pairs
  } else {
    NA_real_
  }

  summary <- data.frame(
    stage = stage,
    sample = sample,
    status = status,
    r1_records = totals$r1_records,
    r2_records = totals$r2_records,
    compared_positions = totals$compared_positions,
    valid_both_records = totals$valid_both_records,
    matching_primary_ids = totals$matching_primary_ids,
    id_mismatches = totals$id_mismatches,
    valid_matched_pairs = totals$valid_matched_pairs,
    equal_length_pairs = totals$equal_length_pairs,
    unequal_length_pairs = totals$unequal_length_pairs,
    unequal_length_percent = unequal_percent,
    malformed_r1 = totals$malformed_r1,
    malformed_r2 = totals$malformed_r2,
    missing_r1 = totals$missing_r1,
    missing_r2 = totals$missing_r2,
    trailing_lines_r1 = totals$trailing_lines_r1,
    trailing_lines_r2 = totals$trailing_lines_r2,
    mate_tags_detected_r1 = totals$mate_tags_detected_r1,
    mate_tags_detected_r2 = totals$mate_tags_detected_r2,
    wrong_mate_tags_r1 = totals$wrong_mate_tags_r1,
    wrong_mate_tags_r2 = totals$wrong_mate_tags_r2,
    exported_failed_r1_records = exported_r1,
    exported_failed_r2_records = exported_r2,
    elapsed_minutes = elapsed_minutes,
    r1_file = r1_file,
    r2_file = r2_file,
    stringsAsFactors = FALSE
  )

  message(sprintf(
    "[%s / %s] %s; valid matched pairs=%s; unequal lengths=%s (%.4f%%)",
    stage, sample, status,
    format(totals$valid_matched_pairs, big.mark = ",", scientific = FALSE),
    format(totals$unequal_length_pairs, big.mark = ",", scientific = FALSE),
    ifelse(is.na(unequal_percent), 0, unequal_percent)
  ))

  list(
    summary = summary,
    length_counts = length_counts_frame(
      length_env, stage, sample, totals$valid_matched_pairs
    ),
    preview = if (length(preview)) do.call(rbind, preview) else NULL
  )
}

write_tsv <- function(data, path) {
  write.table(
    data, file = path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
}

save_qc_plots <- function(summary, length_counts, output_dir) {
  sample_levels <- sort(unique(summary$sample))
  stage_levels <- intersect(c("raw", "trimmed"), unique(summary$stage))

  plot_summary <- summary
  plot_summary$sample <- factor(plot_summary$sample, levels = sample_levels)
  plot_summary$stage <- factor(plot_summary$stage, levels = stage_levels)

  dodge <- ggplot2::position_dodge(width = 0.8)
  p_percentage <- ggplot2::ggplot(
    plot_summary,
    ggplot2::aes(x = sample, y = unequal_length_percent, fill = stage)
  ) +
    ggplot2::geom_col(position = dodge, width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f%%", unequal_length_percent)),
      position = dodge, vjust = -0.35, size = 3
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) sprintf("%.3f%%", x),
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      title = "Unequal mate lengths among valid matched read pairs",
      x = NULL,
      y = "Pairs with R1 length different from R2 (%)",
      fill = "Stage",
      caption = "Unequal lengths are QC observations, not pairing failures."
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank()
    )

  heat <- length_counts
  heat$sample <- factor(heat$sample, levels = sample_levels)
  heat$stage <- factor(heat$stage, levels = stage_levels)
  heat$panel <- interaction(heat$stage, heat$sample, sep = " | ", drop = TRUE)

  p_heatmap <- ggplot2::ggplot(
    heat,
    ggplot2::aes(
      x = r1_length, y = r2_length,
      fill = percent_of_valid_matched_pairs
    )
  ) +
    ggplot2::geom_tile(width = 1, height = 1) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      colour = "#B2182B", linewidth = 0.35
    ) +
    ggplot2::facet_wrap(~panel) +
    ggplot2::scale_fill_gradient(
      low = "#DEEBF7", high = "#08519C", trans = "log10",
      name = "% of valid\nmatched pairs"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Paired-read length combinations",
      x = "R1 length (bp)",
      y = "R2 length (bp)",
      caption = "The dashed diagonal marks equal R1 and R2 lengths; fill uses a log10 scale."
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 8)
    )

  difference_counts <- aggregate(
    pair_count ~ stage + sample + length_difference,
    data = length_counts,
    FUN = sum
  )
  denominators <- summary[, c("stage", "sample", "valid_matched_pairs")]
  difference_counts <- merge(
    difference_counts, denominators,
    by = c("stage", "sample"), all.x = TRUE
  )
  difference_counts$percentage <-
    100 * difference_counts$pair_count / difference_counts$valid_matched_pairs
  difference_counts$sample <- factor(difference_counts$sample, levels = sample_levels)
  difference_counts$stage <- factor(difference_counts$stage, levels = stage_levels)

  p_difference <- ggplot2::ggplot(
    difference_counts,
    ggplot2::aes(x = length_difference, y = percentage, fill = stage)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~sample, scales = "free_y") +
    ggplot2::scale_y_continuous(
      trans = "log10",
      labels = function(x) sprintf("%.4g%%", x)
    ) +
    ggplot2::labs(
      title = "Distribution of mate-length differences",
      x = "R1 length - R2 length (bp)",
      y = "Valid matched pairs (%, log10 scale)",
      fill = "Stage"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))

  save_one <- function(plot, stem, width, height) {
    ggplot2::ggsave(
      file.path(output_dir, paste0(stem, ".png")),
      plot = plot, width = width, height = height,
      units = "in", dpi = 300, bg = "white"
    )
    ggplot2::ggsave(
      file.path(output_dir, paste0(stem, ".pdf")),
      plot = plot, width = width, height = height,
      units = "in", device = "pdf", bg = "white"
    )
  }

  save_one(p_percentage, "unequal_pair_percentage", 9, 5.5)
  save_one(p_heatmap, "pair_length_heatmap", 11, 8)
  save_one(p_difference, "pair_length_difference_distribution", 11, 7)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) usage(0L)
  options <- parse_cli(args)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    abort("The ggplot2 package is required but is not available")
  }

  project_root <- infer_project_root(options$project_root)
  pairs <- discover_pairs(options$stage, options$samples, project_root)

  if (options$dry_run) {
    cat(sprintf("Project root: %s\n\n", project_root))
    print(pairs, row.names = FALSE, right = FALSE)
    cat(sprintf("\nDry run complete: %d paired file set(s); no FASTQ data were read.\n", nrow(pairs)))
    return(invisible(NULL))
  }

  output_dir <- resolve_output_path(options$output_dir, project_root)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    abort("Output directory already exists and is not empty: %s", output_dir)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) abort("Could not create output directory: %s", output_dir)

  settings <- data.frame(
    setting = c(
      "project_root", "stage", "samples", "chunk_size",
      "export_failures", "R_version", "ggplot2_version", "started"
    ),
    value = c(
      project_root, options$stage, options$samples, options$chunk_size,
      options$export_failures, R.version.string,
      as.character(utils::packageVersion("ggplot2")),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
  )
  write_tsv(settings, file.path(output_dir, "run_settings.tsv"))
  write_tsv(pairs, file.path(output_dir, "input_file_pairs.tsv"))

  results <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    results[[i]] <- scan_pair(
      stage = pairs$stage[[i]],
      sample = pairs$sample[[i]],
      r1_file = pairs$r1_file[[i]],
      r2_file = pairs$r2_file[[i]],
      chunk_size = options$chunk_size,
      output_dir = output_dir,
      export_failures = options$export_failures
    )

    current_summary <- do.call(rbind, lapply(results[seq_len(i)], `[[`, "summary"))
    current_counts <- do.call(rbind, lapply(results[seq_len(i)], `[[`, "length_counts"))
    write_tsv(current_summary, file.path(output_dir, "pairing_summary.tsv"))
    write_tsv(current_counts, file.path(output_dir, "pair_length_counts.tsv"))
  }

  summary <- do.call(rbind, lapply(results, `[[`, "summary"))
  length_counts <- do.call(rbind, lapply(results, `[[`, "length_counts"))
  preview_parts <- lapply(results, `[[`, "preview")
  preview_parts <- preview_parts[!vapply(preview_parts, is.null, logical(1L))]

  unequal_summary <- summary[, c(
    "stage", "sample", "status", "valid_matched_pairs",
    "equal_length_pairs", "unequal_length_pairs", "unequal_length_percent"
  )]
  write_tsv(unequal_summary, file.path(output_dir, "unequal_length_summary.tsv"))

  if (length(preview_parts)) {
    write_tsv(
      do.call(rbind, preview_parts),
      file.path(output_dir, "pairing_failures_preview.tsv")
    )
  }

  if (nrow(length_counts)) {
    save_qc_plots(summary, length_counts, output_dir)
  } else {
    warning("No valid matched pairs were available for plotting")
  }

  cat(sprintf("\nCompleted. Results written to:\n%s\n", output_dir))
}

tryCatch(
  main(),
  error = function(error) {
    message("ERROR: ", conditionMessage(error))
    quit(save = "no", status = 1L)
  }
)

