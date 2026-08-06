-- Cell Painting QC over HTTPS, with no Python and no download step.
--
--   duckdb -c ".read demo_remote.sql"      -- needs DuckDB 1.5.1
--
-- Data: two RAB30 sites (ORF + CRISPR) from JUMP cpg0016, 5 channels,
-- 384x384 center crop at native resolution. See make_bundle.py in the repo.

INSTALL duckdb_zarr FROM community;
LOAD duckdb_zarr;

SET VARIABLE base = 'https://raw.githubusercontent.com/USER/REPO/main/jump_rab30_mini';

-- Remote v3 stores need the array node itself, not the group:
--   .../images.zarr/images  with array_path ''   (not .../images.zarr with 'images')
SET VARIABLE arr = getvariable('base') || '/images.zarr/images';

CREATE VIEW px AS
  SELECT s.modality, s.Source, s.Plate, s.Well, c.channel,
         z.dim_2 AS y, z.dim_3 AS x, z.value
  FROM zarr_cells(getvariable('arr'), '') z
  JOIN read_csv(getvariable('base') || '/sites.csv')    s ON s.idx = z.dim_0
  JOIN read_csv(getvariable('base') || '/channels.csv') c ON c.idx = z.dim_1;

.print '== Q1: the notebook''s np.percentile(img, 99.5) vmax, every site x channel at once =='
SELECT modality, Source, Plate, Well, channel,
       round(avg(value), 1) AS mean,
       quantile_cont(value, 0.995)::INT AS vmax_99_5,
       max(value) AS max,
       count(*) FILTER (value = 65535) AS clipped
FROM px GROUP BY ALL ORDER BY modality, channel;

.print ''
.print '== Q2: brightest Mito pixel per site -- argmax + unravel_index, as a query =='
-- MATERIALIZED because filtering on a joined column trips a pushdown bug;
-- filter on dim_N instead and the scan fetches only the chunks it needs.
WITH mito AS MATERIALIZED (
  SELECT dim_0, dim_3 AS x, dim_2 AS y, value
  FROM zarr_cells(getvariable('arr'), '') WHERE dim_1 = 3
)
SELECT s.modality, s.Plate, m.x, m.y, m.value
FROM mito m JOIN read_csv(getvariable('base') || '/sites.csv') s ON s.idx = m.dim_0
QUALIFY row_number() OVER (PARTITION BY s.modality ORDER BY m.value DESC) = 1;
