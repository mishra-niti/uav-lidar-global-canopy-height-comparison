# UAV LiDAR and Global Canopy Height Product Comparison

An R workflow for comparing global canopy-height rasters with a high-resolution UAV LiDAR canopy height model (CHM). The current analysis supports two products, labeled **GCH 2020** and **GFCH 2019**, and evaluates each product against the 95th percentile of valid UAV CHM cells within each global-product pixel.

The workflow calculates pixel-level errors, summary accuracy statistics, and individual and combined error-distribution figures. Input rasters and analysis outputs are intentionally excluded from this repository.

## Analysis overview

For each global canopy-height product, the script:

1. reads and screens the UAV CHM and global-product raster;
2. aligns coordinate reference systems when necessary;
3. identifies the valid area of overlap;
4. assigns UAV CHM cells to global-product pixels;
5. calculates UAV canopy-height P95 and the number of valid UAV cells in each product pixel;
6. retains product pixels containing at least 10 valid UAV cells by default;
7. calculates bias, absolute error, MAE, RMSE, error quantiles, and Pearson correlation; and
8. writes comparison tables, summary tables, histograms, and R session information.

The principal error metric is:

\[
\text{error} = \text{global product height} - \text{UAV CHM P95}
\]

- A **negative** error indicates that the global product underestimates canopy height relative to UAV LiDAR P95.
- A **positive** error indicates that the global product overestimates canopy height relative to UAV LiDAR P95.

## Repository contents

| Path | Purpose |
| --- | --- |
| `multi_product_compare.R` | Complete comparison and figure-generation workflow |
| `data/README.md` | Input-data requirements and naming guidance |
| `outputs/` | Default destination for generated tables and figures |
| `CITATION.cff` | Citation metadata recognized by GitHub |
| `LICENSE` | MIT license for the source code |
| `.gitignore` | Prevents local rasters and generated results from being committed |

## Requirements

- R (a recent version is recommended)
- R packages: `terra`, `dplyr`, `ggplot2`, and `readr`

Install the required packages once from an R console:

```r
install.packages(c("terra", "dplyr", "ggplot2", "readr"))
```

## Input data

The default configuration expects these files:

```text
data/uav_chm.tif
data/gch_2020.tif
data/gfch_2019.tif
```

Different file names or absolute paths can be specified in section 2 of `multi_product_compare.R`.

All input rasters should:

- contain canopy height in metres;
- have a defined coordinate reference system;
- overlap spatially;
- use `NA`/NoData for invalid cells; and
- be clipped, or otherwise limited, to a practical area surrounding the UAV coverage.

The script converts UAV CHM values below 0 m or above 80 m to `NA`. For the global products, values at or below 0 m and above 80 m are converted to `NA`. Modify these thresholds if they are not appropriate for the study system. Treating zero as invalid means that treeless product pixels are excluded.

The data are not distributed with this code. Users must obtain the relevant products independently and cite the precise datasets and versions they use.

## Run the workflow

1. Download or clone this repository.
2. Place the three input rasters in `data/`, or edit their paths in the script.
3. Open a terminal in the repository root and run:

```bash
Rscript multi_product_compare.R
```

The script can also be opened and sourced from RStudio, with the repository root set as the working directory.

## Outputs

By default, all results are written to `outputs/`.

| Output | Description |
| --- | --- |
| `GCH_2020_vs_UAV_P95_pixel_comparison.csv` | Pixel-level GCH–UAV comparison |
| `GFCH_2019_vs_UAV_P95_pixel_comparison.csv` | Pixel-level GFCH–UAV comparison |
| `*_vs_UAV_P95_summary.csv` | Per-product summary statistics |
| `GCH_GFCH_vs_UAV_P95_combined_summary.csv` | Combined product summary |
| `GCH_GFCH_vs_UAV_P95_combined_pixel_comparison.csv` | Combined pixel-level table |
| `*_minus_UAV_P95_histogram.png/.pdf` | Per-product error distributions |
| `GCH_GFCH_minus_UAV_P95_combined_histogram.png/.pdf` | Faceted comparison figure |
| `session_info.txt` | R version, platform, and loaded-package information |

The pixel-level tables include the product-pixel identifier, product height, UAV P95, number of valid UAV cells, signed error, and absolute error.

## Statistical definitions

- **Mean bias:** mean of `product height - UAV P95`.
- **Median bias:** median of `product height - UAV P95`.
- **MAE:** mean absolute error.
- **RMSE:** root mean square error.
- **Pearson r:** linear correlation between product height and UAV P95.
- **Error quantiles:** 5th, 25th, 75th, and 95th percentiles of signed error.

## Important considerations

- The reference is UAV **P95**, not the maximum UAV canopy height.
- The default requirement of 10 valid UAV cells per product pixel is a quality-control choice and may need sensitivity testing.
- Edge pixels may contain only partial UAV coverage; inspect `uav_n` or apply a stricter coverage rule when appropriate.
- If coordinate systems differ, the script uses nearest-neighbor reprojection to avoid smoothing height values. For a formal validation study, document all alignment and resampling decisions and consider preprocessing the product rasters to a common analysis extent.
- Product uncertainty, acquisition-date differences, geolocation error, canopy change, and differing canopy-height definitions can all contribute to observed disagreement.
- Record the full name, source, version, spatial resolution, acquisition year, and citation for every input product in the associated manuscript or data record.

## Citation

Use GitHub's **Cite this repository** option, which reads `CITATION.cff`. The global canopy-height datasets must be cited separately using their authoritative publications or data records.

## License

The source code is released under the [MIT License](LICENSE). This license does not apply to any UAV LiDAR data, global canopy-height products, or third-party datasets used with the workflow.

