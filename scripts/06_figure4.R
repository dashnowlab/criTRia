library(grid)

draw_evidence_table <- function() {
  rows <- data.frame(
    section = rep(
      c("Genetic/bioinformatic\nevidence", "Clinical evidence", "Other"),
      times = c(11, 1, 2)
    ),
    information = c(
      "Inheritance pattern",
      "Genomic coordinates",
      "Genomic feature\n(e.g., exon, intron, UTR)",
      "Size of repeat (case)",
      "Size of repeat (control)",
      "Proposed allelic thresholds",
      "Pathogenic",
      "Alternatives",
      "Flanking",
      "Interruptions",
      "Methodology of discovery",
      "Case history",
      "Ancestry",
      "Standardized terminology (such\nas HPO terms)"
    ),
    essential = c(
      TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
      FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE
    ),
    motif = c(rep(FALSE, 6), rep(TRUE, 4), rep(FALSE, 4)),
    height = c(20, 20, 32, rep(20, 10), 30)
  )
  
  x <- c(4, 137, 221, 305, 389, 473)
  y <- c(6, 24, 24 + cumsum(rows$height))
  
  colors <- c(
    header = "#84BBC3",
    information = "#B8DDE7",
    status = "#D3E8EE",
    border = "#000000",
    status_grid = "#9DA4A6"
  )
  
  pushViewport(viewport(xscale = c(0, 478), yscale = c(0, 333)))
  on.exit(popViewport())
  
  # Convert top-down image coordinates to grid's bottom-up coordinates.
  y_pos <- function(value) unit(333 - value, "native")
  
  cell <- function(left, right, top, bottom, fill, label = "",
                   bold = FALSE, border = colors[["border"]]) {
    center_x <- unit(mean(c(left, right)), "native")
    center_y <- y_pos(mean(c(top, bottom)))
    
    grid.rect(
      x = center_x, y = center_y,
      width = unit(right - left, "native"),
      height = unit(bottom - top, "native"),
      gp = gpar(fill = fill, col = border, lwd = 0.6)
    )
    
    if (nzchar(label)) {
      grid.text(
        label, x = center_x, y = center_y,
        gp = gpar(
          fontfamily = "sans", fontsize = 8.2, lineheight = 0.95,
          fontface = if (bold) "bold" else "plain"
        )
      )
    }
  }
  
  # Information spans the main label and motif subcolumn.
  header_edges <- x[c(1, 2, 4, 5, 6)]
  headings <- c("", "Information", "Essential", "If relevant")
  
  for (i in seq_along(headings)) {
    cell(
      header_edges[i], header_edges[i + 1], y[1], y[2],
      colors[["header"]], headings[i], bold = TRUE
    )
  }
  
  # Merge the category cells across consecutive rows.
  sections <- rle(rows$section)
  section_end <- cumsum(sections$lengths)
  section_start <- section_end - sections$lengths + 1L
  
  for (i in seq_along(sections$values)) {
    cell(
      x[1], x[2], y[section_start[i] + 1], y[section_end[i] + 2],
      colors[["header"]], sections$values[i], bold = TRUE
    )
  }
  
  motif_rows <- which(rows$motif)
  cell(
    x[2], x[3], y[min(motif_rows) + 1], y[max(motif_rows) + 2],
    colors[["information"]], "Motif"
  )
  
  for (i in seq_len(nrow(rows))) {
    top <- y[i + 1]
    bottom <- y[i + 2]
    label_left <- if (rows$motif[i]) x[3] else x[2]
    
    cell(
      label_left, x[4], top, bottom,
      colors[["information"]], rows$information[i]
    )
    
    cell(x[4], x[5], top, bottom,
         colors[["status"]], border = colors[["status_grid"]])
    cell(x[5], x[6], top, bottom,
         colors[["status"]], border = colors[["status_grid"]])
    
    # Draw the checkmark without relying on a font glyph.
    mark_col <- if (rows$essential[i]) 4 else 5
    mark_x <- mean(x[c(mark_col, mark_col + 1)])
    mark_y <- mean(c(top, bottom))
    
    grid.lines(
      x = unit(mark_x + c(-5, -2, 6), "native"),
      y = y_pos(mark_y + c(1, 4, -4)),
      gp = gpar(col = colors[["border"]], lwd = 0.75)
    )
  }
  
  # Extend the section dividers across the status columns.
  for (edge in y[section_start + 1]) {
    grid.lines(
      x = unit(x[c(4, 6)], "native"),
      y = y_pos(rep(edge, 2)),
      gp = gpar(col = colors[["border"]], lwd = 0.6)
    )
  }
  
  cell(x[1], x[6], y[1], tail(y, 1), fill = NA)
}

# Preview in RStudio.
if (interactive()) {
  grid.newpage()
  draw_evidence_table()
}

# Export a PDF at the reference aspect ratio.
pdf(
  "paper/fig4_minimum_evidence.pdf",
  width = 1434 / 300,
  height = 999 / 300,
  bg = "white"
)
grid.newpage()
draw_evidence_table()
invisible(dev.off())