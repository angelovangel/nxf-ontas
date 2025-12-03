#! /usr/bin/env Rscript

library(vroom)
library(dplyr)
library(reactable)
library(htmlwidgets)
library(htmltools) # Added library for easily adding HTML elements like a title


# takes all tsv files in current dir 

files <- list.files(path = ".", pattern = "*.tsv", full.names = T)
df <- vroom(files, delim = "\t", id = 'file') %>%
  mutate(file = basename(file)) %>%
  mutate(file = tools::file_path_sans_ext(file))


summary_df <- df %>%
  group_by(file) %>%
  summarise(
    Total_Bases_Ontarget = sum(bases, na.rm = T),
    Avg_Target_CovX = ifelse(sum(region_len) > 0, Total_Bases_Ontarget / sum(region_len, na.rm = T), 0)
  )
    
# Define column specifications for formatting and alignment
columns_def <- list(
    file = colDef(
        name = "Sample", 
        align = "left",
        # Ensure the column doesn't get treated as a number
        format = colFormat(digits = 0),
    ),
    Total_Bases_Ontarget = colDef(
        name = "Total bases on target",
        align = "right",
        # Format as integer with comma separator
        format = colFormat(separators = TRUE, digits = 0),
        style = list(fontFamily = "monospace")
    ),
    Avg_Target_CovX = colDef(
        name = "Avg. coverage (X)",
        align = "right",
        # Format as number with two decimal places
        format = colFormat(digits = 2),
        style = list(fontFamily = "monospace")
    )
)

# Define a modern, clean theme
modern_theme <- reactableTheme(
  # General styles
  color = "hsl(0, 0%, 15%)", # Dark grey text
  backgroundColor = "hsl(0, 0%, 100%)", # White background
  borderColor = "hsl(0, 0%, 90%)", # Very light grey border
  # Header styles
  headerStyle = list(
    background = "hsl(0, 0%, 98%)", # Almost white header
    fontWeight = 600, # Medium weight
    fontSize = "14px",
    textTransform = "uppercase",
    # Add subtle border to separate header
    borderBottom = "2px solid hsl(0, 0%, 90%)"
  ),
  # Row highlight style - using a subtle, modern blue/teal
  highlightColor = "rgba(0, 150, 136, 0.08)", # Subtle teal highlight (Original: #009688)
  # Input and button styles
  inputStyle = list(
    borderColor = "hsl(0, 0%, 80%)",
    borderRadius = "4px"
  ),
  # Pagination button styles
  pageButtonHoverStyle = list(backgroundColor = "hsl(0, 0%, 90%)"),
  pageButtonActiveStyle = list(backgroundColor = "hsl(0, 0%, 80%)")
)


# Create the reactable table object with the new theme and 50% width
react_table <- reactable(
    summary_df,
    columns = columns_def,
    # Make the table searchable and sortable by default
    searchable = TRUE,
    filterable = FALSE, 
    highlight = TRUE,
    defaultPageSize = 10,
    theme = modern_theme, 
    width = "75%"
)


# --- WRAP TABLE WITH A TITLE USING HTMLTOOLS::TAGLIST ---
# Define the title element
report_title <- tags$h2(
    "Target Coverage Summary Report", 
    style = "font-family: Arial, sans-serif; 
             color: hsl(0, 0%, 10%); 
             margin-bottom: 20px;"
)

# Combine the title and the reactable into a single list
report_content <- tagList(
    report_title,
    react_table
)

# --- Save the reactable widget as a self-contained HTML file ---
output_file <- "coverage_summary_report.html"

# Use the combined 'report_content' and set the 'title' argument
htmlwidgets::saveWidget(
    react_table, 
    file = output_file, 
    selfcontained = TRUE,
    title = "Target Coverage Report" # Sets the browser tab title
)

cat(paste0("Interactive HTML report saved to: ", output_file, "\n"))