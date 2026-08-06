-- Cell Painting QC with no Python in it.
--
--   curl -LO <url>/jump_rab30_mini.zip && unzip -q jump_rab30_mini.zip
--   duckdb -c ".read demo.sql"          -- needs DuckDB 1.5.1 or newer
--
-- Bundle: two RAB30 sites (ORF + CRISPR) from JUMP cpg0016, 5 channels,
-- 384x384 center crop at native resolution. See make_bundle.py inside the zip.

INSTALL duckdb_zarr FROM community;
LOAD duckdb_zarr;

-- zarr_cells projects the (site, channel, y, x) array to rows of
-- (dim_0, dim_1, dim_2, dim_3, value) -- bare integer indices. The joins to
-- sites.csv and channels.csv are the only thing that make the array mean anything.
CREATE VIEW px AS
  SELECT s.modality, s.Source, s.Plate, s.Well, c.channel,
         z.dim_2 AS y, z.dim_3 AS x, z.value
  FROM zarr_cells('jump_rab30_mini/images.zarr', 'images') z
  JOIN 'jump_rab30_mini/sites.csv'    s ON s.idx = z.dim_0
  JOIN 'jump_rab30_mini/channels.csv' c ON c.idx = z.dim_1;

.print '== Q1: np.percentile(img, 99.5) from the notebook, every site x channel at once =='
SELECT modality, Source, Plate, Well, channel,
       round(avg(value), 1) AS mean,
       quantile_cont(value, 0.995)::INT AS vmax_99_5,
       max(value) AS max,
       count(*) FILTER (value = 65535) AS clipped
FROM px GROUP BY ALL ORDER BY modality, channel;

.print ''
.print '== Q2: brightest Mito pixel per site -- argmax + unravel_index, as a query =='
-- MATERIALIZED because filtering on a joined column trips a pushdown bug;
-- filter on dim_N instead and the scan reads only the chunks it needs.
WITH mito AS MATERIALIZED (
  SELECT dim_0, dim_3 AS x, dim_2 AS y, value
  FROM zarr_cells('jump_rab30_mini/images.zarr', 'images') WHERE dim_1 = 3
)
SELECT s.modality, s.Plate, m.x, m.y, m.value
FROM mito m JOIN 'jump_rab30_mini/sites.csv' s ON s.idx = m.dim_0
QUALIFY row_number() OVER (PARTITION BY s.modality ORDER BY m.value DESC) = 1;

.print ''
.print '== Q3: 32x32 binning of the CRISPR DNA channel, top-left 128px =='
SELECT dim_3 // 32 AS bx, dim_2 // 32 AS biny, round(avg(value)) AS m
FROM zarr_cells('jump_rab30_mini/images.zarr', 'images')
WHERE dim_0 = 1 AND dim_1 = 1 AND dim_2 < 128 AND dim_3 < 128
GROUP BY ALL ORDER BY biny, bx LIMIT 8;
