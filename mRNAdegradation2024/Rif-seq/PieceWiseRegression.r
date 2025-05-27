
# .r script by Seunghyeon Kim.
# This .r file converts linear mRNA level to log mRNA level and generates graphs containing log mRNA level and segmented fit into a subdirectory.
# Put .r script and .xlsx files in the same folder and run.
# Check Example.xlsx and change it to your data. You also have to check Sampling time. (Line 18)
# Visit https://cran.r-project.org/package=segmented for segmented package of R.
# Last updated: May 27, 2025.

# Load packages
library(readxl)
library(openxlsx)
library(segmented)

# Create a folder to save plots
plot_dir <- "plots/"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# Sampling time [minutes]
time <- c(0, 0.5, 1, 2, 3, 4, 5, 6, 8, 10)

# Open .xlsx file
file_path <- "Example.xlsx" # Example.xlsx contains linear scale mRNA level
wb <- loadWorkbook(file_path)
data <- read_excel(file_path, sheet = "Sheet1", col_names = TRUE)
name <- data[[1]]
data2 <- data[, -1]
data3 <- t(data2)
data3 <- log(data3 / data3[1])
data3[is.na(data3) | is.nan(data3) | is.infinite(data3)] <- NA

all_outputs <- list() # An empty list

for (i in 1:ncol(data3)) {
  data4 <- data.frame(time, data5 = data3[, i])
  data4 <- na.omit(data4)
  model <- lm(data5 ~ time, data = data4) # lm = linear fitting model
  lm.coef <- summary(model)$coef  # fitting results
  davies_test <- davies.test(model, seg.Z = ~ time, k = 20) # davies test
  davies_test_pvalue <- davies_test$p.value # p value

  gene_name <- ifelse(is.na(name[i]) | name[i] == "", paste0("Sample_", i), name[i])

  if (davies_test_pvalue < 0.05) {
    Piecewise <- TRUE
    lm_segment <- segmented(model, seg.Z = ~ time, npsi = 1)

    breakpoint <- lm_segment$psi[1, "Est."]
    slopes <- slope(lm_segment)
    slope1 <- slopes$time[1, "Est."]
    slope2 <- slopes$time[2, "Est."]
    intercept <- intercept(lm_segment)
    intercepts1 <- intercept$time[1, "Est."]
    intercepts2 <- intercept$time[2, "Est."]

    # Save graphs to .pdf files, -1 means 1 breakpoint
    pdf(file.path(plot_dir, paste0(gene_name, "-1.pdf")))
    par(pin = c(3, 2.5))
    plot(time, data4$data5, col = "blue", pch = 16, cex = 1.5,
         xlab = "Time after adding rifampicin [min]", ylab = "Log mRNA level [AU]", main = paste(gene_name, " - Segmented"))
    plot(lm_segment, col = "red", lwd = 2, add = TRUE)
    dev.off()

  } else {
    Piecewise <- FALSE
    breakpoint <- 0  # Linear fit
    slope1 <- lm.coef[2, 1]
    slope2 <- NA  # There's no slope2 since it's linear fit
    intercepts1 <- lm.coef[1, 1]
    intercepts2 <- NA  # There's no intercepts2 since it's linear fit

    # Save graphs to .pdf files, -0 means 0 breakpoints (linear fit)
    pdf(file.path(plot_dir, paste0(gene_name, "-0.pdf")))
    par(pin = c(3, 2.5))
    plot(time, data4$data5, col = "blue", pch = 16, cex = 1.5,
         xlab = "Time after adding rifampicin [min]", ylab = "Log mRNA level [AU]", main = paste(gene_name, " - Linear"))
    abline(model, col = "red", lwd = 2)
    dev.off()
  }

  # Save fitting results to the list
  all_output <- list(davies_test_pvalue, Piecewise, breakpoint, c(slope1, slope2), intercepts1, intercepts2)
  all_outputs[[i]] <- all_output
}

# Dataframe
result_df <- data.frame(
  davies_test_pvalue = sapply(all_outputs, function(x) x[[1]]),  # p value from davis test
  Piecewise = sapply(all_outputs, function(x) x[[2]]),           # Piecewise
  breakpoint = sapply(all_outputs, function(x) x[[3]]),          # Breakpoint
  slope1 = sapply(all_outputs, function(x) ifelse(is.na(x[[4]][1]), NA, x[[4]][1])),
  slope2 = sapply(all_outputs, function(x) ifelse(is.na(x[[4]][2]), NA, x[[4]][2])),
  intercept1 = sapply(all_outputs, function(x) ifelse(is.na(x[[5]]), NA, x[[5]])),
  intercept2 = sapply(all_outputs, function(x) ifelse(is.na(x[[6]]), NA, x[[6]]))
)

Gene_name <- data[[1]]
final_df <- cbind(Gene_name, result_df)
addWorksheet(wb, "Regression Results") # Add 'Regression Results' sheet in the original .xlsx file
writeData(wb, "Regression Results", final_df)
saveWorkbook(wb, file_path, overwrite = TRUE)
