Work in progress: tests and a single-command pipeline are being added

## Week 1: Learning the FCS Format on a Single File

Before analysing a real dataset, I wanted to understand what an FCS file actually contains and how to handle it in Python. I used one example file from an Attune cytometer, from the [`tlnagy/fcsexamples`](https://github.com/tlnagy/fcsexamples) repository (originally from FlowRepository), and read it with [FlowKit](https://github.com/whitews/FlowKit).

### What's in the file

The file holds 22,188 events, recorded in 2010. Each event is a row, and each column is a channel:

- `Event` and `Time`: a counter and an acquisition timestamp, not measurements of the cell. I leave both out of the analysis.
- `FSC` and `SSC`: forward scatter (size) and side scatter (granularity).
- Six fluorescence channels: `BL1–BL3` from the blue laser and `VL1–VL3` from the violet laser.
- Each channel in three forms: area (`-A`), height (`-H`) and width (`-W`).

Two things were missing, and both shaped what the file could be used for:

- **No marker names.** The `$PnS` labels, where the stain for each channel is normally recorded, were empty. So there was no way to know what any fluorescence channel measured, and therefore nothing biological to gate.
- **No compensation matrix.** I wrote a small check that searches the file's metadata for any spillover or compensation key. It found none, so compensation wasn't possible for this file.

The metadata also varied from what I expected: the `cyt` field, which normally names the instrument, held the acquisition software version instead. A tool that handles files from different instruments can't assume any metadata field exists, so my code reads these fields without relying on them.

### Seeing the data

A plain scatter plot of FSC-A against SSC-A squashed most of the events into one corner: a few extreme events stretched the axes to around 10 million, while most cells sat far below that. Switching to a density plot (a 2D histogram) on its own didn't help, for the same reason.

Two changes made the data readable: limiting the plotted range to the 0.5th–99.5th percentiles, and using a logarithmic colour scale. That revealed a dense main population and a sparser group of larger, more granular events.

![FSC-A vs SSC-A, density plot with log colour scale](flow-cytomtery-gating/images/attune_fsc_ssc_density.png)

### A surprise in the fluorescence channels

Plotting BL1-A against VL1-A showed almost every event on a single tight diagonal line, grouped into distinct, evenly spaced clusters. In a stained cell sample, two channels from different lasers shouldn't track each other this closely. This pattern is typical of multi-peak calibration beads, which are made at a series of fixed brightness levels and fluoresce across many channels at once. Together with the tight scatter profile and the missing marker names, this suggests the file is a bead sample used for instrument quality control rather than cells. The file came without a sample description, so I can't confirm it.

### Transformation

On a linear scale, the low-intensity clusters were crammed together near zero. Applying an arcsinh transformation (FlowKit's `AsinhTransform`, with the top of scale set to the file's measurement range) spread them out into about eight clearly separated clusters, roughly evenly spaced along the diagonal. Even spacing on a log-like scale means each cluster is brighter than the last by a similar multiple, which is how multi-peak beads are designed, and which strengthens the bead interpretation.

![BL1-A vs VL1-A before and after arcsinh transformation](flow-cytomtery-gating/images/attune_bl1_vl1_raw_vs_asinh.png)

Two things I learned the hard way:

- **Transformations are for fluorescence channels.** Scatter channels are normally kept on a linear scale. FlowKit reflects this: `apply_transform` leaves scatter channels untouched by default.
- **`apply_transform` changes the sample, not an existing DataFrame.** My first transformed plot looked identical to the raw one because I had plotted a DataFrame created *before* the transform. The transformed values have to be requested afterwards with `as_dataframe(source="xform")`.

### What Week 1 established

By the end of the week, I could read any FCS file, inspect its channels and metadata, check for a compensation matrix, plot scatter and density views, and transform fluorescence data. The file itself had no biology to gate, which made the next step clear: I needed a real, documented cell dataset with known markers and, ideally, expert gating to compare against.

## Week 2: Gating a Real Dataset

### The data

I used the haematopoietic stem cell transplant (HSCT) dataset from the FlowCAP-I challenge, available on FlowRepository (FR-FCM-ZZY6): 30 FCS files, each one sample. According to the dataset description, the samples are mouse bone marrow from transplant experiments in Connie Eaves's laboratory at the Terry Fox Laboratory, BC Cancer. Before staining, erythroid precursors were removed by magnetic depletion of Ter119-labelled cells.

Each file has six channels, all as heights: `FSC-H`, `SSC-H` and `FL1-H` to `FL4-H`. As with the Attune file, the files don't record which marker is in which channel, and there's no compensation matrix or instrument metadata. I found the panel in the documentation of the R package [`ks`](https://cran.r-project.org/package=ks), which redistributes part of this dataset:

| Channel | Dye | Marker | What it tells you |
|---|---|---|---|
| FL1-H | FITC | CD45.1 | Donor-derived white blood cells |
| FL2-H | PE | Mac1 (CD11b) | Myeloid cells |
| FL3-H | PI | Live/dead | Dead cells (PI only enters damaged cells) |
| FL4-H | APC | CD45.2 | Recipient-derived white blood cells |

This mapping came from a secondary source, so I treated it as an assumption until Week 3, when I confirmed it value by value against the `ks` data.

All fluorescence values run from 0 to 1023, which means the instrument had already log-scaled them. Applying a further transformation would distort the data, so I gated the raw values.

### The biology

This is a mouse transplant experiment. Donor and recipient strains carry different versions of the CD45 protein, CD45.1 and CD45.2, so the two can be told apart by staining. Weeks after the transplant, the fraction of donor cells in the recipient's bone marrow shows how well the transplanted stem cells engrafted. Donor myeloid cells are particularly informative: myeloid cells live only days, so donor myeloid cells found long after the transplant must have been freshly made by donor stem cells.

### Choosing the gates

Each gate answers one biological question, and the population I keep is decided by the biology, not by where most events sit. I chose cut-offs by eye on file 001, then checked them on other files. All histograms use 256 bins over the full 0–1023 range with a log-scaled y-axis, so that small populations stay visible and files can be compared directly.

**Cells versus debris (FSC-H vs SSC-H).** There were no events below about 235 on FSC-H. The instrument was most likely set to record only events above a forward-scatter threshold, so most debris was never recorded, and I didn't apply a debris gate.

**Live versus dead (FL3-H, PI).** The PI histogram had no clean valley: a long, smoothly falling tail up to about 600–650, then a sparse region, then a small cluster near 1000. Because PE from the Mac1 stain can leak into the PI detector, I plotted FL2-H against FL3-H. There was no diagonal smear, so the tail isn't spillover from Mac1. I set the cut-off at 700, at the start of the sparse region. To check how much the choice mattered, I counted live cells in file 001 at several cut-offs:

| PI cut-off | Live cells | Change from 700 |
|---|---|---|
| 600 | 9,721 | −133 (1.3% of events) |
| 700 | 9,854 | — |
| 800 | 9,886 | +32 (0.3%) |
| 900 | 9,901 | +47 (0.5%) |

Above 700, the exact value barely matters. Below it, the cut-off starts cutting into the tail of live cells, which is why I didn't go lower.

**Recipient versus the rest (FL4-H, CD45.2).** Two clearly separated groups, with an empty band around 500–600. Cut-off: 550. The recipient population sat in the same place in every normally stained file, so one cut-off works across the batch.

**Donor versus CD45-negative (FL1-H, CD45.1).** Among cells low on CD45.2, some were bright on CD45.1 and some dim. A cell low on both CD45 versions isn't a donor white blood cell; it's most likely a cell carrying little or no CD45, such as erythroid cells that survived the Ter119 depletion. So I split this group in two rather than treating it all as donor. In files 022, 025 and 028, both groups sat in consistent places, but the valley between them shifted, from about 400 to about 500, depending on which group was larger. A cut-off of 450 sat in the sparse region in all three.

**Myeloid versus non-myeloid (FL2-H, Mac1).** A clear valley between about 400 and 550. Cut-off: 500. Mac1-negative cells are a mixture (B cells and their precursors, stem and progenitor cells, immature myeloid cells and others), so I describe them as Mac1-negative rather than lymphoid.

Two features appear in every histogram: tall spikes at exactly 0 (events below the detector's floor) and small piles at 1023 (events above its ceiling). Both are real events, not missing data.

The final gating tree:

```
All events
└── Live (FL3-H < 700)
    ├── Recipient (FL4-H ≥ 550)                       → Mac1+ / Mac1−
    ├── Donor (FL4-H < 550 and FL1-H ≥ 450)            → Mac1+ / Mac1−
    └── CD45-negative (FL4-H < 550 and FL1-H < 450)
```

### What I got wrong first

My first attempt at choosing cut-offs went wrong in instructive ways. I drew a tight box around the densest cluster on FSC/SSC, which would have thrown away real granular cells. I kept the low Mac1 peak, which is the non-myeloid population, the opposite of what I wanted. And I plotted with percentile-clipped axes, which had been useful for the Attune file but here hid the bright donor cells. Plotting the full 0–1023 range with more bins and a log-scaled y-axis, and deciding each gate from its biological question, fixed all three.

### Results for file 001

| Population | Cells | % of parent |
|---|---|---|
| Live (of 9,936 events) | 9,854 | 99.2% |
| Recipient (of live) | 9,461 | 96.0% |
| Donor (of live) | 4 | 0.04% |
| CD45-negative (of live) | 389 | 3.9% |
| Mac1+ within recipient | 710 | 7.5% |
| Mac1+ within donor | 1 | (4 cells, too few for a percentage) |
| Donor share of all myeloid cells | 1 of 711 | 0.14% |

File 001 has essentially no donor engraftment. Early on, I saw a small group of cells low on CD45.2 and was tempted to call them donor cells. Splitting that group by CD45.1 showed that 389 of its 393 cells are dim on both CD45 versions, so they aren't donor white blood cells at all. Only 4 cells have the donor profile. This is why I kept the cautious label "CD45.2-low" until CD45.1 could confirm what those cells were.

### Batch results

All 30 files were gated with the same cut-offs, and the results saved to `results/hsct_gating_summary.csv`.

**The last four files turned out to be controls.** Files 027–030 didn't look like transplant samples. Checking all four fluorescence channels in each, a clear pattern emerged:

| File | Channels with a positive population | Most likely |
|---|---|---|
| 027 | CD45.2 only | Single-stain control for CD45.2 (APC) |
| 028 | CD45.1 only | Single-stain control for CD45.1 (FITC) |
| 029 | a small Mac1 group only | Possibly a single-stain control for Mac1 (PE), weakly stained |
| 030 | none | Unstained control |

All four had a normal PI pattern, which fits PI being added to every tube. This is the standard set of compensation controls: an unstained tube plus one tube per fluorescent stain, usually acquired together at the end of a run. The dataset description doesn't mention controls, so these identifications are inferences from the data, but four files fitting one pattern makes the interpretation fairly strong. I've marked them with a `likely_control` column in the results table and left them out of the summaries below.

**The controls also provided a check on the gates.** Each single-stain control shows where cells that are negative for a given stain end up:

| Channel | Negative cells end around | My cut-off |
|---|---|---|
| FL1 (CD45.1), from 027 | ~350 | 450 |
| FL2 (Mac1), from 027 and 030 | ~300–350 | 500 |
| FL4 (CD45.2), from 028 | ~290 | 550 |

Every cut-off sits above where the negative population ends, which is the same check a lab would make when setting gates from controls.

Across the remaining 26 transplant samples:

| Measure | Median | Range |
|---|---|---|
| Live (% of events) | 98.8% | 97.2–99.2% |
| Recipient (% of live) | 91.8% | 16.0–99.6% |
| Donor (% of live) | 0.03% | 0–82.4% |
| CD45-negative (% of live) | 2.8% | 0.4–20.2% |
| Mac1+ within recipient | 47.0% | 2.7–67.2% |
| Donor share of myeloid cells | 0.02% | 0–94.7% |

I report medians and ranges rather than means because the values are very uneven: only 10 of the 26 samples have more than 1% donor cells, while most have almost none. Donor contribution ranges from essentially absent, as in file 001, to the large majority of cells, as in file 022.

The share of Mac1-positive cells among recipient cells varied from 2.7% to 67.2%. The lowest values come from samples dominated by donor cells, such as 005 and 006. That fits the biology: once donor stem cells have taken over myeloid production, few recipient cells are myeloid. File 001's 7.5% also sits within this range, so it's a normal sample rather than a sign that the channel mapping was wrong.

Percentages from very small groups are unreliable. Where a group has only a handful of cells, like file 001's 4 donor cells, the count matters more than the percentage.

### Why the cut-offs weren't automated

I explored finding each file's cut-off automatically, by smoothing each histogram and locating the valley between its two main peaks, and decided against it. In a real lab, gates are set using control samples, which show where "negative" ends for each stain. When I started, the dataset appeared to have no controls, so an automated cut-off would have had nothing to be calibrated against. I used fixed gates chosen by eye instead, applied identically to every sample, which is standard practice within one experiment.

The controls I later found in files 027–030 confirmed those gates: every cut-off sits above the negative population in the relevant control. With documented controls, the next step would be to set the cut-offs from the controls directly rather than by eye.

## Week 3: Validating the Gates Against Expert Labels

### Finding the labels

FlowCAP-I scored automated methods against expert manual gating, so expert labels exist for this dataset, but they weren't among the files on FlowRepository. I found them in the R package [`ks`](https://cran.r-project.org/package=ks), which includes four of the HSCT samples as an example dataset: 39,128 events from mice 5, 6, 9 and 12, with four fluorescence values and an expert label for every event. I exported the data from R to a CSV file:

```r
install.packages("ks")
library(ks)
data(hsct)
write.csv(hsct, "data/hsct_ks_labels.csv", row.names = FALSE)
```

### Checking the labels line up with my files

The labels are only useful if each one can be attached to the right event. For each of the four mice, I compared the `ks` data with the matching FCS file (005, 006, 009 and 012): the number of events matched, and every fluorescence value matched, in the same order. That means the labels can be attached to my events by position.

It also settled a question left open since Week 2. The `ks` columns are named by marker (`FITC.CD45.1`, `PE.Ly65Mac1`, `PI.LiveDead`, `APC.CD45.2`), and their values are identical to FL1-H, FL2-H, FL3-H and FL4-H. So the channel mapping I had assumed is confirmed.

### What the expert labels mean

The labels are numbered 0–4, but the numbers are just codes, and the documentation doesn't say which population each one is. I worked it out from the data, by profiling each class and by cross-tabulating it against my own groups.

The first cross-table, with all four files pooled, made no sense: recipient cells appeared to be split between two expert classes, and one class seemed to contain both recipient and donor cells. Making a separate table for each file resolved it. Within each file, each of my groups fell almost entirely into one expert class, but the numbers meant different populations in different files. Class 1, for example, is recipient Mac1− cells in files 005 and 009, but donor Mac1+ cells in 006 and recipient Mac1+ cells in 012. The labels were assigned file by file, which works for FlowCAP's purpose of scoring algorithms but means they can't be pooled.

Across all four files, the experts used the same populations I did:

- **recipient and donor cells, each split into Mac1-positive and Mac1-negative**
- **class 0, a catch-all:** it contains essentially all the cells I classed as dead, most of my CD45-negative cells, and a few borderline cells from every other group

### Scoring the agreement

For each expert class except class 0, I matched it to the group of mine with the largest overlap in that file, then calculated:

- **precision:** of the cells in my group, the share the experts also put in that class
- **recall:** of the cells in the experts' class, the share I found
- **F-measure:** a combined score from 0 to 1, the same measure FlowCAP used

| File | Expert class | My group | Precision | Recall | F-measure |
|---|---|---|---|---|---|
| 005 | 1 | recipient Mac1− | 0.969 | 1.000 | 0.984 |
| 005 | 2 | donor Mac1− | 0.976 | 0.995 | 0.986 |
| 005 | 3 | donor Mac1+ | 0.957 | 0.957 | 0.957 |
| 006 | 1 | donor Mac1+ | 0.937 | 0.985 | 0.960 |
| 006 | 2 | recipient Mac1− | 0.957 | 1.000 | 0.978 |
| 006 | 3 | donor Mac1− | 0.974 | 0.996 | 0.985 |
| 009 | 1 | recipient Mac1− | 0.975 | 0.998 | 0.986 |
| 009 | 2 | recipient Mac1+ | 0.981 | 0.993 | 0.987 |
| 009 | 3 | donor Mac1+ | 0.995 | 0.996 | 0.996 |
| 009 | 4 | donor Mac1− | 0.983 | 0.995 | 0.989 |
| 012 | 1 | recipient Mac1+ | 0.978 | 0.999 | 0.988 |
| 012 | 2 | recipient Mac1− | 0.971 | 0.998 | 0.984 |
| 012 | 3 | donor Mac1− | 0.966 | 1.000 | 0.983 |
| 012 | 4 | donor Mac1+ | 0.979 | 0.996 | 0.987 |

In file 005, expert class 4 contains a single cell. A score from one cell means nothing, so I excluded it.

**Across the 14 scored populations, the F-measure ranges from 0.957 to 0.996, with a median of 0.985.** Fixed gates chosen by eye, without controls, reproduce expert manual gating very closely.

### Where my gates and the experts disagree

Recall is almost always higher than precision: I find nearly every cell the experts placed in a population, but my groups also include a few cells the experts left in class 0. These are most likely borderline cells near the gate edges, which the experts chose not to assign. My gates are slightly more inclusive than theirs rather than wrong.

The one consistent disagreement in the other direction is in file 005, where 54 cells the experts called donor fell into my CD45-negative group. They were donor by the experts' judgement but just below my CD45.1 cut-off of 450, so for this file the cut-off may be slightly high. I set all the cut-offs before looking at the expert labels and didn't change them afterwards, so these scores are a fair test of gates chosen independently.

### Limitations

- Expert labels were available for only 4 of the 26 transplant samples, so the validation covers those four. The other 22 were gated with the same cut-offs.
- The meaning of each expert class was inferred from the data rather than documented.
- Expert labels are a reference, not absolute truth: different analysts gate the same file slightly differently, which is part of why FlowCAP existed.