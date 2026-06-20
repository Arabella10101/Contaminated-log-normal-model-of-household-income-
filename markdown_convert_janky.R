library(reticulate)

# 1. Ensure the correct virtual environment is used
use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)

# 2. Initialize the converter
mid <- import("markitdown")
md_converter <- mid$MarkItDown()

# 3. Define your directories
source_dir <- "C:/Users/arabe/Documents/Research_Project/research papers/to convert"
target_dir <- "C:/Users/arabe/Documents/Research_Project/research papers/markdown files"

# Create the target directory if it doesn't exist yet
if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}

# 4. Get a list of all PDF files in the source folder
pdf_files <- list.files(path = source_dir, pattern = "\\.pdf$", full.names = TRUE, ignore.case = TRUE)

# 5. Loop through each file and convert
if (length(pdf_files) == 0) {
  cat("No PDF files found in the source directory.\n")
} else {
  for (pdf_path in pdf_files) {
    # Extract the base file name (e.g., "paper1.pdf")
    file_name <- basename(pdf_path)
    
    # Change the extension to .md (e.g., "paper1.md")
    md_name <- sub("\\.pdf$", ".md", file_name, ignore.case = TRUE)
    
    # Create the full output path
    output_md_path <- file.path(target_dir, md_name)
    
    cat("Converting:", file_name, "...\n")
    
    # Wrap in tryCatch so one bad PDF doesn't crash the whole loop
    tryCatch({
      # Run conversion
      result <- md_converter$convert(pdf_path)
      
      # Save the markdown file
      writeLines(result$text_content, con = output_md_path, useBytes = TRUE)
      
      cat("Successfully saved to:", md_name, "\n\n")
    }, error = function(e) {
      cat("❌ Error converting", file_name, ":", conditionMessage(e), "\n\n")
    })
  }
  cat("All done! Check your markdown folder.")
}