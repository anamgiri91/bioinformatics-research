# ================================================================
# Task 13
# Create Sample-Level Flag Summary Excel Workbook
#
# Input:
#   Task13_Normal_selfref_flags.csv
#   Task13_Normal_extref_flags.csv
#   Task13_Tumor_selfref_flags.csv
#   Task13_Tumor_extref_flags.csv
#
# Output:
#   Task13_SampleFlagSummary.xlsx
#
# For each of the 53 samples:
#   - number of -1 flags
#   - number of  0 flags
#   - number of +1 flags
#   - percentage of -1 flags
#   - percentage of  0 flags
#   - percentage of +1 flags
#   - total flagged CpGs (-1 + +1)
#   - total flagged percentage
# ================================================================


# ------------------------------------------------
# 1. Load package
# ------------------------------------------------

if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop(
        "openxlsx is not installed.\n",
        "Run install.packages('openxlsx') first."
    )
}

library(openxlsx)


# ------------------------------------------------
# 2. Directories
# ------------------------------------------------

result_dir <- "/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results"

output_file <- file.path(
    result_dir,
    "Task13_SampleFlagSummary.xlsx"
)


# ------------------------------------------------
# 3. Input files
# ------------------------------------------------

normal_self_file <- file.path(
    result_dir,
    "Task13_Normal_selfref_flags.csv"
)

normal_ext_file <- file.path(
    result_dir,
    "Task13_Normal_extref_flags.csv"
)

tumor_self_file <- file.path(
    result_dir,
    "Task13_Tumor_selfref_flags.csv"
)

tumor_ext_file <- file.path(
    result_dir,
    "Task13_Tumor_extref_flags.csv"
)


# ------------------------------------------------
# 4. Function to process one Task 13 CSV
# ------------------------------------------------

get_sample_summary <- function(file, dataset_name) {

    cat("\n========================================\n")
    cat("Processing:", dataset_name, "\n")
    cat("File:", basename(file), "\n")
    cat("========================================\n")

    if (!file.exists(file)) {
        stop("File does not exist: ", file)
    }

    # Read CSV
    dat <- read.csv(
        file,
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    cat(
        "Dimensions:",
        nrow(dat),
        "x",
        ncol(dat),
        "\n"
    )

    # Expected dimensions
    if (ncol(dat) != 65) {
        stop(
            "Expected 65 columns but found ",
            ncol(dat),
            " in ",
            basename(file)
        )
    }

    # ------------------------------------------------
    # Sample columns
    #
    # Columns 1-12 = metadata and CpG-level summaries
    # Columns 13-65 = 53 sample flags
    # ------------------------------------------------

    sample_flags <- dat[, 13:65]

    sample_names <- names(sample_flags)

    cat(
        "Number of samples:",
        length(sample_names),
        "\n"
    )

    # Convert flags to numeric
    sample_flags[] <- lapply(
        sample_flags,
        function(x) as.numeric(as.character(x))
    )


    # ------------------------------------------------
    # Count each flag for each sample
    # ------------------------------------------------

    count_neg1 <- colSums(
        sample_flags == -1,
        na.rm = TRUE
    )

    count_zero <- colSums(
        sample_flags == 0,
        na.rm = TRUE
    )

    count_pos1 <- colSums(
        sample_flags == 1,
        na.rm = TRUE
    )


    # ------------------------------------------------
    # Number of usable CpGs
    #
    # Important for external reference:
    # unmatched CpGs are NA.
    # ------------------------------------------------

    usable_cpgs <- colSums(
        !is.na(sample_flags)
    )


    # ------------------------------------------------
    # Percentages
    # ------------------------------------------------

    pct_neg1 <- 100 * count_neg1 / usable_cpgs

    pct_zero <- 100 * count_zero / usable_cpgs

    pct_pos1 <- 100 * count_pos1 / usable_cpgs


    # ------------------------------------------------
    # Total flagged
    #
    # A flag is considered an outlier if:
    #     flag == -1 OR flag == +1
    # ------------------------------------------------

    total_flagged <- count_neg1 + count_pos1

    total_flagged_pct <- (
        100 * total_flagged / usable_cpgs
    )


    # ------------------------------------------------
    # Create summary table
    # ------------------------------------------------

    result <- data.frame(

        Sample = sample_names,

        Dataset = dataset_name,

        Usable_CpGs = usable_cpgs,

        Flag_neg1_Count = count_neg1,

        Flag_0_Count = count_zero,

        Flag_pos1_Count = count_pos1,

        Flag_neg1_Pct = round(pct_neg1, 3),

        Flag_0_Pct = round(pct_zero, 3),

        Flag_pos1_Pct = round(pct_pos1, 3),

        Total_Flagged = total_flagged,

        Total_Flagged_Pct =
            round(total_flagged_pct, 3),

        stringsAsFactors = FALSE
    )

    return(result)
}


# ------------------------------------------------
# 5. Process all four datasets
# ------------------------------------------------

normal_self <- get_sample_summary(
    normal_self_file,
    "Normal - Self Reference"
)

normal_ext <- get_sample_summary(
    normal_ext_file,
    "Normal - External Reference"
)

tumor_self <- get_sample_summary(
    tumor_self_file,
    "Tumor - Self Reference"
)

tumor_ext <- get_sample_summary(
    tumor_ext_file,
    "Tumor - External Reference"
)


# ------------------------------------------------
# 6. Print sample summaries
# ------------------------------------------------

cat("\n\n========== NORMAL SELF ==========\n")
print(normal_self)

cat("\n\n========== NORMAL EXTERNAL ==========\n")
print(normal_ext)

cat("\n\n========== TUMOR SELF ==========\n")
print(tumor_self)

cat("\n\n========== TUMOR EXTERNAL ==========\n")
print(tumor_ext)


# ------------------------------------------------
# 7. Create combined long-format table
# ------------------------------------------------

all_results <- rbind(
    normal_self,
    normal_ext,
    tumor_self,
    tumor_ext
)


# ------------------------------------------------
# 8. Sort samples naturally
# ------------------------------------------------

sort_samples <- function(x) {

    sample_num <- as.numeric(
        gsub("[^0-9]", "", x$Sample)
    )

    x[
        order(sample_num),
        ,
        drop = FALSE
    ]
}


normal_self <- sort_samples(normal_self)
normal_ext  <- sort_samples(normal_ext)

tumor_self <- sort_samples(tumor_self)
tumor_ext  <- sort_samples(tumor_ext)

all_results <- rbind(
    normal_self,
    normal_ext,
    tumor_self,
    tumor_ext
)


# ------------------------------------------------
# 9. Create wide comparison table
# ------------------------------------------------

make_wide <- function(data, prefix) {

    output <- data[, c(
        "Sample",
        "Flag_neg1_Count",
        "Flag_0_Count",
        "Flag_pos1_Count",
        "Flag_neg1_Pct",
        "Flag_0_Pct",
        "Flag_pos1_Pct",
        "Total_Flagged",
        "Total_Flagged_Pct"
    )]

    names(output)[-1] <- paste0(
        prefix,
        "_",
        names(output)[-1]
    )

    output
}


normal_wide <- merge(
    make_wide(normal_self, "Self"),
    make_wide(normal_ext, "External"),
    by = "Sample"
)

tumor_wide <- merge(
    make_wide(tumor_self, "Self"),
    make_wide(tumor_ext, "External"),
    by = "Sample"
)


# Sort N1-N53 / T1-T53
normal_wide <- normal_wide[
    order(
        as.numeric(
            gsub("[^0-9]", "", normal_wide$Sample)
        )
    ),
]

tumor_wide <- tumor_wide[
    order(
        as.numeric(
            gsub("[^0-9]", "", tumor_wide$Sample)
        )
    ),
]


# ------------------------------------------------
# 10. Create full comparison table
# ------------------------------------------------

all_wide <- merge(
    normal_wide,
    tumor_wide,
    by.x = "Sample",
    by.y = "Sample",
    all = TRUE
)


# ------------------------------------------------
# 11. Create workbook
# ------------------------------------------------

wb <- createWorkbook()


# ------------------------------------------------
# 12. Styles
# ------------------------------------------------

title_style <- createStyle(
    fontSize = 14,
    textDecoration = "bold"
)

header_style <- createStyle(
    fontSize = 11,
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = "Bottom"
)

highlight_style <- createStyle(
    fgFill = "#FFF2CC",
    textDecoration = "bold"
)


# ------------------------------------------------
# 13. Function to add worksheet
# ------------------------------------------------

add_table_sheet <- function(
    wb,
    sheet_name,
    data,
    title = NULL
) {

    addWorksheet(
        wb,
        sheet_name
    )

    if (!is.null(title)) {

        writeData(
            wb,
            sheet_name,
            title,
            startRow = 1,
            startCol = 1
        )

        addStyle(
            wb,
            sheet_name,
            title_style,
            rows = 1,
            cols = 1
        )

        start_row <- 3

    } else {

        start_row <- 1

    }

    writeData(
        wb,
        sheet_name,
        data,
        startRow = start_row,
        startCol = 1
    )

    addStyle(
        wb,
        sheet_name,
        header_style,
        rows = start_row,
        cols = 1:ncol(data),
        gridExpand = TRUE
    )

    setColWidths(
        wb,
        sheet_name,
        cols = 1:ncol(data),
        widths = "auto"
    )

    freezePane(
        wb,
        sheet_name,
        firstActiveRow = start_row + 1,
        firstActiveCol = 2
    )

    # Highlight N14-N17 or T14-T17
    standout <- c(
        "N14", "N15", "N16", "N17",
        "T14", "T15", "T16", "T17"
    )

    rows_to_highlight <- which(
        data$Sample %in% standout
    )

    if (length(rows_to_highlight) > 0) {

        addStyle(
            wb,
            sheet_name,
            highlight_style,
            rows = rows_to_highlight + start_row,
            cols = 1:ncol(data),
            gridExpand = TRUE,
            stack = TRUE
        )
    }
}


# ------------------------------------------------
# 14. Summary sheet
# ------------------------------------------------

add_table_sheet(
    wb,
    "Summary",
    all_wide,
    "Task 13 - Sample-Level Flag Summary"
)


# ------------------------------------------------
# 15. Normal sheet
# ------------------------------------------------

add_table_sheet(
    wb,
    "Normal",
    normal_wide,
    "Normal - Self vs External Reference"
)


# ------------------------------------------------
# 16. Tumor sheet
# ------------------------------------------------

add_table_sheet(
    wb,
    "Tumor",
    tumor_wide,
    "Tumor - Self vs External Reference"
)


# ------------------------------------------------
# 17. Individual analysis sheets
# ------------------------------------------------

add_table_sheet(
    wb,
    "Normal_Self",
    normal_self,
    "Normal - Self Reference"
)

add_table_sheet(
    wb,
    "Normal_Ext",
    normal_ext,
    "Normal - External Reference"
)

add_table_sheet(
    wb,
    "Tumor_Self",
    tumor_self,
    "Tumor - Self Reference"
)

add_table_sheet(
    wb,
    "Tumor_Ext",
    tumor_ext,
    "Tumor - External Reference"
)


# ------------------------------------------------
# 18. Save workbook
# ------------------------------------------------

saveWorkbook(
    wb,
    output_file,
    overwrite = TRUE
)


# ------------------------------------------------
# 19. Final message
# ------------------------------------------------

cat("\n\n")
cat("====================================================\n")
cat("TASK 13 SAMPLE SUMMARY COMPLETE\n")
cat("====================================================\n")
cat("Output:\n")
cat(output_file, "\n\n")

cat("Sheets created:\n")
cat("  Summary\n")
cat("  Normal\n")
cat("  Tumor\n")
cat("  Normal_Self\n")
cat("  Normal_Ext\n")
cat("  Tumor_Self\n")
cat("  Tumor_Ext\n")
cat("\n====================================================\n")
