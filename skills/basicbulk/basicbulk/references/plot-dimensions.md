# Plot Dimensions

Use these default R Markdown chunk dimensions when generating basicBulk reports. `N` means number of contrasts.

| Plot | fig.width | fig.height | Function args | Notes |
|------|-----------|------------|---------------|-------|
| `PCAgrid` | 9 | 7 | | 12x9 for many groups |
| `DEGbar` | 4 | 3 | | Increase width for long group names |
| `basicVolcano` | N<=2: 6, N<=4: 9, N>=5: 12 | 6 for single, 4 for faceted | | Width scales with contrasts |
| `basicMA` | N<=2: 6, N<=4: 9, N>=5: 12 | 6 for single, 4 for faceted | | Same scaling as volcano |
| `basicHeatmap` | 7 | 6 | | Per-contrast heatmaps can be narrower in width (6) |
| `basicUpset` | N<=3: 9, N>=4: 12 | N<=3: 6, N>=4: 8 | | Overlaps section, UpSet tab; scales with contrast count |
| Overlap Expression (`genePlot`) | 9 | 6 | `ncol = 3` | Overlaps section; top shared genes as a dot grid, scale height with gene count |
| `gg4way` | 5.5 | 6 | | Fixed; one plot per tab |
| `corTriplot` | 4 | 2.75 | | LogFC Correlations Overview; compact by default — keep it small, do not enlarge at low/moderate N (this is the default that renders right). Resize up (toward ~9×8.5) only when a genuinely high contrast count makes the cells unreadable, or split by control instead |
| `gg4way` combined grid | 12 | scales with rows | | Optional small-N grid (`N<=4`); height ~4 per row of panels |
| `topGSEAplot` | N<=3: 6, N>=4: 9 | N<=3: 4, N>=4: 6 | | Both scale with contrasts |
| `GSEAplotList` | 6 | 3 | | Per-plot in tabset pills |
| Pathway heatmap (`basicHeatmap`) | 7 | 6 | `pathways = TRUE` | Pathways section, Heatmap tab |
| Gene set scores (`genePlot`) | 9 | 5 | | Gene Set Scores section; scale height with pathway count |
| `genePlot` | 9 | 6 | `nrow = 2` | Use `nrow = 2` for multi-gene grids; scale height with gene count |
| `geneHistogram` | 9 | 7 | | Per-gene histogram plots |

## High-cardinality categorical legends

Past ~20 categories the `ggsci` d3 palette runs out of distinct colours **and** a default legend overflows the figure edges. These are two separate problems — fix both:

- **Colour:** switch from `ggsci::scale_*_d3()` to `pals::glasbey()` (up to ~32 distinct colours), e.g. `ggplot2::scale_color_manual(values = pals::glasbey(nlevels(f)) %>% purrr::set_names(levels(f)))`.
- **Layout:** a palette swap alone still overflows — force a fixed wrap and shrink the keys so a wide legend fits: `+ ggplot2::guides(<aes> = ggplot2::guide_legend(ncol = N, theme = ggplot2::theme(legend.byrow = TRUE))) + ggplot2::theme(legend.text = ggplot2::element_text(size = 7), legend.key.size = grid::unit(9, "pt"))`, and give the figure enough width/height for the legend strip.

Both are needed past ~20 categories (e.g. 20+ crossed groups, or a dozen long treatment names).

Use `fig.show="hide"` for chunks that produce plots as side effects, including `augere.de::runVoom()`, `augere.gsea::runContrast()` and `augere.gsea::runPrecomputed()`.

Set `include=FALSE` for pure computation chunks and for the heatmap build chunk (it emits `png()`/`pdf()` + ggsave "Saving … image" device messages that `fig.show="hide"` does not suppress; `include=FALSE` still runs the chunk, keeping the built objects available to the print chunk).
