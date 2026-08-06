# Querying Cell Painting pixels as SQL, with duckdb_zarr

A small, real dataset plus two SQL scripts, for kicking the tires on
[`duckdb_zarr`](https://duckdb.org/community_extensions/extensions/duckdb_zarr).

The point of the demo: an n-dimensional array is a function from coordinates to
values, and its relational projection is one row per cell.
That view is useless for compute (convolution, segmentation, registration) but
it is exactly right for the part of imaging work that is really database work -
select, filter, join, reduce.
In Cell Painting the axes are *semantic*: `dim_0` is a perturbation and `dim_1`
is a stain, so the interesting questions all cross from pixels into plate
metadata.

## Run it

Needs DuckDB **1.5.1 or newer** - `duckdb_zarr` is published for 1.5.1, 1.5.2 and
1.5.3 (linux/macOS/Windows), and nothing earlier.

```sh
duckdb -c ".read demo.sql"          # local files in this repo
duckdb -c ".read demo_remote.sql"   # same thing over raw.githubusercontent.com
```

## The data

`jump_rab30_mini/` - two real RAB30 sites from JUMP
[cpg0016](https://github.com/jump-cellpainting/datasets), the ones used by the
JUMP hub's
[`14_display_perturbation_images`](https://broadinstitute.github.io/jump_hub/howto/notebooks/14_display_perturbation_images.html)
notebook:

| | modality | source | plate | well | site |
|-|-|-|-|-|-|
| `dim_0 = 0` | ORF | `source_4` | `BR00123947` | G01 | 1 |
| `dim_0 = 1` | CRISPR | `source_13` | `CP-CC9-R1-05` | L07 | 0 |

- `images.zarr` - `(site, channel, y, x)` = `(2, 5, 384, 384)` `uint16`, gzip.
  384x384 center crop at **native resolution**, so pixel statistics (including
  saturation) are real rather than an artifact of downsampling.
- `sites.csv`, `channels.csv` - the lookup tables that give the integer
  dimension indices their meaning.
- `make_bundle.py` - regenerates the whole thing from JUMP. Run `./make_bundle.py`.

2.2 MB total. Cell Painting Gallery images are CC0.

## What the queries show

`demo.sql` Q1 replaces the notebook's per-image `np.percentile(img, 99.5)` -
the call that picks `vmax` for display - with one `GROUP BY` over every site and
channel at once:

```
┌──────────┬───────────┬──────────────┬─────────┬─────────┬────────┬───────────┬────────┬─────────┐
│ modality │  Source   │    Plate     │  Well   │ channel │  mean  │ vmax_99_5 │  max   │ clipped │
├──────────┼───────────┼──────────────┼─────────┼─────────┼────────┼───────────┼────────┼─────────┤
│ CRISPR   │ source_13 │ CP-CC9-R1-05 │ L07     │ AGP     │  326.6 │       693 │   1774 │       0 │
│ CRISPR   │ source_13 │ CP-CC9-R1-05 │ L07     │ DNA     │  184.1 │       640 │   1000 │       0 │
│ CRISPR   │ source_13 │ CP-CC9-R1-05 │ L07     │ ER      │  805.9 │      2235 │   3527 │       0 │
│ CRISPR   │ source_13 │ CP-CC9-R1-05 │ L07     │ Mito    │ 2095.3 │      4673 │   6992 │       0 │
│ CRISPR   │ source_13 │ CP-CC9-R1-05 │ L07     │ RNA     │ 1196.4 │      3214 │   4362 │       0 │
│ ORF      │ source_4  │ BR00123947   │ G01     │ AGP     │ 4039.3 │     17062 │  33165 │       0 │
│ ORF      │ source_4  │ BR00123947   │ G01     │ DNA     │ 2140.3 │     26996 │  39759 │       0 │
│ ORF      │ source_4  │ BR00123947   │ G01     │ ER      │ 4143.7 │     31301 │  50892 │       0 │
│ ORF      │ source_4  │ BR00123947   │ G01     │ Mito    │ 3056.6 │     17388 │  31646 │       0 │
│ ORF      │ source_4  │ BR00123947   │ G01     │ RNA     │ 5006.5 │     41262 │  65535 │      12 │
└──────────┴───────────┴──────────────┴─────────┴─────────┴────────┴───────────┴────────┴─────────┘
```

An order-of-magnitude intensity difference between sources and 12 saturated
pixels, in one statement. Q2 is `argmax` + `unravel_index` as a window function.
Q3 is spatial binning as a `GROUP BY`.

## Notes for the maintainers

Everything below was hit while building this. Versions: DuckDB 1.5.1 and 1.5.2,
`duckdb_zarr` from community extensions, `zarr-python` 3.x, linux_amd64.

**1. The Cell Painting Gallery's OME-Zarr is unreadable.**
This is the big one. `cpg0004-lincs` is the only Cell Painting Gallery dataset
converted to OME-Zarr, and every variant of it (`images_zarr`,
`images_zarr_050`, `images_zarr_withdownscale8`) is written by `bioformats2raw`
as Zarr v2 with `blosc`/`lz4`, `>u2`, chunks `[1,1,1,1024,1024]`. Four errors
stack up, each behind the last:

```
remote, as-is                      -> Remote Zarr v2 stores currently require consolidated metadata (.zmetadata)
local + hand-written .zmetadata    -> Zarr metadata key "shuffle" is not a string
  + shuffle as a string            -> Blosc codec metadata is missing integer typesize
  + typesize                       -> Only Blosc with cname=zstd is currently supported
```

The terminal blocker is the codec. `blosc`/`lz4` support would make the whole
gallery queryable in place. The two metadata errors look like v2 blosc config
being parsed with v3 field conventions - v2 spells `shuffle` as an int and has
no `typesize`. And `bioformats2raw` output has no `.zmetadata` and never will,
so remote v2 needs to work without it.

**2. `zstd` is rejected, and it is `zarr-python`'s default.**
A store written the obvious way (`zarr.create_array(...)` with no `compressors=`)
fails with `does not yet support Zarr v3 codec: zstd`. `blosc` fails too.
`gzip` and uncompressed work - hence `GzipCodec()` in `make_bundle.py`.

**3. Remote and local take different arguments.**
Group discovery is not implemented for remote v3, so the same store needs two
different call shapes:

```sql
zarr_cells('jump_rab30_mini/images.zarr', 'images')         -- local
zarr_cells('https://.../images.zarr/images', '')            -- remote
```

**4. Filtering on a joined column breaks the scan.**

```sql
-- Invalid Input Error: Unsupported pushed filter type for zarr_cells
SELECT ... FROM zarr_cells(...) z JOIN chans c ON c.idx = z.dim_1 WHERE c.channel = 'Mito'
```

DuckDB pushes a semi-join filter into the scan and it is rejected. Filtering on
`dim_1 = 3` directly is fine, so the workaround is a `MATERIALIZED` CTE - see
`demo.sql` Q2. This is the roughest edge for the join-heavy use case, since
joining dimension indices to metadata is the main reason to want this at all.

**5. `zarr` vs `duckdb_zarr`.**
The [`zarr`](https://duckdb.org/community_extensions/extensions/zarr) extension
(xqlsystems) has a docs page but no published binaries for any platform or
version - every `community-extensions.duckdb.org` URL 404s. This demo uses
[`duckdb_zarr`](https://github.com/WayScience/duckdb_zarr) (WayScience), which
is published for 1.5.1 and up.
