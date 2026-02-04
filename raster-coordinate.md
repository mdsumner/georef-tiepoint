# Raster-to-Model Coordinate Transformations: A Historical and Technical Analysis

## Status: DRAFT FOR REVIEW AND EXPANSION
**Author:** Michael Sumner (AAD)  
**Reviewer checklist:** GDAL devs, xarray/CF community  
**Last updated:** February 2026

---

## Executive Summary

This document explores the historical development and current landscape of methods for associating raster image coordinates with geospatial reference systems. It traces the lineages of tiepoints, affine geotransforms, GCPs, RPCs, geolocation arrays, and coordinate variable approaches—and considers how these divergent paths affect modern formats like GeoZarr.

The central thesis: **GeoTIFF ossified around the simplest affine case while carrying vestigial complexity; Zarr/CF has the structural flexibility to handle all cases but the conventions are immature. Neither format has elegantly unified the spectrum from simple affine through curvilinear grids. GDAL has become the practical Rosetta Stone that papers over these divergent lineages.**

---

## 1. GeoTIFF Tiepoints: Historical Context

### 1.1 Origins (1994-1995)

GeoTIFF emerged from practical necessity. By the early 1990s, satellite imagery from SPOT and Landsat generated vast volumes of raster data, but TIFF files required external "world files" or sidecar metadata for geospatial context. SPOT Image Corporation proposed an initial GeoTIFF structure in September 1994, limited to UTM.

**Key figures:**
- Dr. Niles Ritter (NASA JPL) - primary author and coordinator
- Mike Ruth (SPOT Image) - key contributor
- Ed Grissom (Intergraph Corporation) - initial leadership
- Roger Lott (EPSG) - projection parametrization consultation

The GeoTIFF Working Group (formed 1994) rapidly expanded to 140+ subscribers. Version 1.0 was released November 1995 (Ritter & Ruth, 1995).

**Source:** OGC GeoTIFF 1.1 Standard; Wikipedia; spec history documents

### 1.2 The Tiepoint Mechanism

The spec defines three coordinate transformation tags:

1. **ModelTiepointTag** (33922): Stores raster→model tiepoint pairs as `(I,J,K,X,Y,Z)` tuples
2. **ModelPixelScaleTag** (33550): Defines pixel spacing `(ScaleX, ScaleY, ScaleZ)`
3. **ModelTransformationTag** (34264): Full 4×4 affine transformation matrix

**The crucial relationship:**

```
Single tiepoint + ModelPixelScaleTag ≡ 6-parameter affine geotransform
```

The math is unambiguous: if only one tiepoint `(I,J,K,X,Y,Z)` and scale `(Sx, Sy, Sz)` are specified:

```
| X |   | Sx   0   0   Tx |   | I |
| Y | = | 0  -Sy   0   Ty | × | J |
| Z |   | 0    0  Sz   Tz |   | K |
| 1 |   | 0    0   0    1 |   | 1 |

where: Tx = X - I*Sx, Ty = Y + J*Sy, Tz = Z - K*Sz
```

### 1.3 Multiple Tiepoints: "Not Baseline GeoTIFF"

**This is explicitly stated in the spec:**

> "Case 5: The raster data cannot be fit into the model space with a simple affine transformation (rubber-sheeting required). Use only the ModelTiepoint tag, and specify as many tiepoints as your application requires. **Note, however, that this is not a Baseline GeoTIFF implementation, and should not be used for interchange**; it is recommended that the image be geometrically rectified first, and put into a standard projected coordinate system."

**Source:** GeoTIFF spec 2.6.3 "Cookbook for Defining Transformations"

The spec also notes: "For orthorectification or mosaicking applications a large number of tiepoints may be specified on a mesh over the raster image. **However, the definition of associated grid interpolation methods is not in the scope of the current GeoTIFF spec.**"

Even more striking, the spec explicitly refuses to commit to any interpolation semantics even for simple cases:

> "However, tiepoints are only to be considered exact at the points specified; **thus defining such a set of bounding tiepoints does not imply that the model space locations of the interior of the image may be exactly computed by a linear interpolation of these tiepoints.**"

This is a deliberate design choice: the format stores knowledge (these pixels correspond to these coordinates), but the application applies judgment (how to interpolate between them).

Even more striking, the spec explicitly refuses to commit to ANY interpolation semantics, even for simple cases:

> "However, tiepoints are only to be considered exact at the points specified; **thus defining such a set of bounding tiepoints does not imply that the model space locations of the interior of the image may be exactly computed by a linear interpolation of these tiepoints.**"

This is remarkable: even if you provide 3 corner tiepoints, you *cannot assume* the interior is linearly interpolated. The spec stores points but refuses to define what happens between them.

**Assessment:** Multiple tiepoints were *mentioned* as theoretically possible, but:
- No interpolation method was standardized
- Explicitly called "not Baseline"
- Explicitly discouraged for interchange
- The OGC 1.1 standardization (2019) repeated this guidance
- The spec explicitly punts on interpolation semantics entirely

**⚠️ NEEDS VERIFICATION:** Were multiple tiepoints ever seriously used anywhere? What did GDAL do when encountering them historically? 

### 1.4 The Design Philosophy: Carry vs Warp

The GeoTIFF spec's refusal to define interpolation methods reveals a deliberate architectural boundary:

- **The format carries knowledge**: "these pixels map to these coordinates"
- **The application applies judgment**: polynomial order, TPS, bilinear, etc.
- **For interchange, just rectify first**: the spec literally recommends this

This maps directly to GDAL's fundamental distinction between `RasterIO` and `Warp`:

| Operation | What it does | Transformation handling |
|-----------|--------------|------------------------|
| **RasterIO** | Reads/writes pixels | Carries geotransform/CRS metadata through unchanged |
| **Warp** | Reprojects imagery | Seeks out transformation info (geotransform, GCPs, RPCs, geolocation arrays), applies user-specified method |

RasterIO doesn't care about rectification—it reads pixels and passes metadata along. Warp is where the transformation *method* gets chosen: `-order 1`, `-tps`, `-rpc`, etc.

GeoTIFF was designed for the RasterIO world: store the georeferencing information, let the application decide what to do with it. The spec explicitly avoids becoming a Warp specification.

### 1.5 The 3D Vestige

The tiepoint format `(I,J,K,X,Y,Z)` carries a Z dimension "in anticipation of future support for 3D digital elevation models and vertical coordinate systems." In practice:

- K and Z are almost always set to zero
- ScaleZ is "primarily used to map the pixel value of a digital elevation model into the correct Z-scale, and so for most other purposes this value should be zero"

**Assessment:** This appears to be forward-looking design that never materialized as intended.

---

## 2. The GCP/RPC/Geolocation Array Divergence

### 2.1 Ground Control Points (GCPs)

GCPs emerged as a separate mechanism in GDAL, distinct from GeoTIFF tiepoints:

- Stored as dataset-level metadata
- Used with polynomial transformations (order 1-3) or Thin Plate Spline
- Require warping to produce georeferenced output

**GDAL's approach:** `gdalwarp` with `-order N` or `-tps` flags, or `gdal_translate` to attach GCPs, followed by warping.

**Relationship to tiepoints:** Both are conceptually "point pairs mapping pixel to georef coordinates," but:
- GeoTIFF tiepoints were intended for *intrinsic* georeferencing
- GCPs are for *post-hoc* correction/registration

**⚠️ NEEDS VERIFICATION:** Were GeoTIFF tiepoints and GDAL GCPs ever intended to converge? Or were they always seen as separate domains? 

### 2.2 Rational Polynomial Coefficients (RPCs)

RPCs became the satellite vendor's solution for the sensor-model problem. GDAL RFC 22 (circa 2006-2007?) added formal RPC support:

- Stored as metadata in the "RPC" domain
- Describe the physical relationship between image coordinates and ground coordinates via rational polynomials
- Support orthorectification with DEM

**Key products using RPCs:** GeoEye, DigitalGlobe, SPOT (some products)

**Source:** GDAL RFC 22: RPC Georeferencing

### 2.3 Geolocation Arrays

GDAL RFC 4 added geolocation array support for swath and curvilinear data:

- Common in AVHRR, Envisat, HDF4/5, netCDF products
- Full 2D lat/lon arrays per pixel (or subsampled)
- Stored in the "GEOLOCATION" domain

**Key quote from RFC 4:** "It is common in AVHRR, Envisat, HDF and netCDF data products to distribute geolocation for raw or projected data in this manner, and current approaches to representing this as very large numbers of GCPs, or greatly subsampling the geolocation information to provide more reasonable numbers of GCPs are inadequate for many applications."

**Assessment:** Geolocation arrays solved a real problem that neither GeoTIFF tiepoints nor sparse GCPs could handle: high-resolution curvilinear grids from satellite swath data.

---

## 3. GDAL as the Rosetta Stone

### 3.1 The Unified Transformer API

GDAL's `GDALCreateGenImgProjTransformer2()` presents a unified view across:

- `GEOTRANSFORM` - Simple affine
- `GCP_POLYNOMIAL` - GCP-based polynomial (order 1-3)
- `GCP_TPS` - GCP-based Thin Plate Spline
- `GCP_HOMOGRAPHY` - GCP-based homography (added later)
- `GEOLOC_ARRAY` - Full geolocation arrays
- `RPC` - Rational Polynomial Coefficients

**The `METHOD` option allows forcing a specific approach when multiple are available.**

### 3.2 The Precedence Problem

By default:
1. Geotransform takes precedence if available
2. GCPs used if `bGCPUseOK=TRUE` (which `gdalwarp` always sets)
3. Geolocation arrays only used if nothing else available

This creates practical issues for data with multiple georeferencing mechanisms.

**⚠️ NEEDS VERIFICATION:** Is this precedence behavior intentional design or historical accident? 

---

## 4. CF Conventions and the Coordinate Array Lineage

### 4.1 Different Origins, Different Mental Models

CF conventions emerged from climate/ocean modeling where:

- Lat/lon coordinate arrays were the natural primitive
- Regular grids used 1D coordinate variables (COARDS style)
- Curvilinear/swath data used 2D auxiliary coordinate variables
- CRS was an afterthought

But more fundamentally, the *mental model* is different from the remote sensing world:

**Remote sensing (GeoTIFF lineage):**
- Raw image → needs geometric correction → rectified product is the goal
- Geotransform describes the *corrected* state
- You warp TO a map
- Unrectified data is an intermediate state

**Climate/ocean modeling (CF lineage):**
- Model output on a grid → the grid IS the coordinate system
- Coordinates describe where things ARE, not where they need to be transformed to
- You're already in "model space" - regridding is a deliberate analytical choice, not a correction
- The native grid is the authoritative representation

When an ocean modeler has NEMO output on a rotated pole grid with 2D `nav_lat`/`nav_lon` arrays, they're not thinking "this needs rectification." They're thinking "these are my coordinates, I might regrid to a regular lat/lon for comparison with observations, but that's an analytical operation with explicit choices about interpolation method."

The CF `coordinates` attribute pointing to 2D lat/lon arrays isn't a deficient form of georeferencing waiting to be warped—it's the *complete* description of where the data lives. The curvilinear grid is the native representation, not an intermediate state.

This explains why:
- CDO and xarray regridding are explicit operations with user-specified methods
- There's no "auto-warp on read" in the netCDF ecosystem  
- `pyresample` exists as a separate tool for when you actually want to resample
- Nobody in the CF world talks about "rectification"

**The `grid_mapping` approach:**
- A variable (often called `crs`) holds attributes describing the CRS
- Data variables reference it via `grid_mapping="crs"`
- Projection parameters spread across multiple attributes

### 4.2 The `crs_wkt` Addition

CF 1.7 (circa 2017) added support for WKT strings via `crs_wkt` attribute:

**From CF conventions:**
> "The crs_wkt attribute is intended to act as a supplement to other single-property CF grid mapping attributes; it is not intended to replace those attributes. If data producers omit the single-property grid mapping attributes in favour of the compound crs_wkt attribute, software which cannot interpret crs_wkt will be unable to use the file."

**The tension:** Some argued (GitHub issue #222) that `crs_wkt` should be authoritative since:
- Many grid_mapping implementations are botched
- Copy-paste WKT is easier than mapping parameters
- PROJ strongly discourages PROJ strings for CRS representation

**Current CF convention:** If both exist and conflict, single-property attributes take precedence.

**Assessment:** This is a reasonable compromise but creates redundancy and potential inconsistency.

### 4.3 Known Problems with CF Grid Mapping

From the CF conventions wiki (Mapping from CF Grid Mapping Attributes to CRS WKT Elements):

> "There are several parameters that do not match up and in several cases a grid mapping does not exist. This is problematic for users who wish to convert back and forth between the two."

**pyproj perspective:** The library now includes `CRS.to_cf()` and `CRS.from_cf()` methods, but these mappings remain complex and error-prone.

---

## 5. GeoZarr: Current State and Critique

### 5.1 What GeoZarr Is Trying to Do

The GeoZarr SWG Charter (OGC) states:

> "The goal of the GeoZarr specification is to establish flexible and inclusive conventions for the Zarr cloud-native format, specifically designed to meet the diverse requirements within the geospatial domain."

**Planned conformance classes include:**
- Multiple related variables with heterogeneous coordinates
- Multiple resolutions (multiscale arrays)
- Multi-dimensional optimizations
- Earth observation products (multispectral bands)

### 5.2 Current Approach

GeoZarr (draft spec) relies on CF conventions:
- `grid_mapping` variable with CF attributes + optional `crs_wkt`
- 1D coordinate variables for rectilinear grids
- 2D auxiliary coordinates for curvilinear grids
- xarray's `_ARRAY_DIMENSIONS` convention

**GDAL's Zarr driver approach:**
- Uses a `_CRS` attribute with `wkt`, `url`, or `projjson` keys
- Reads `GeoTransform` from spatial_ref variable when available
- Falls back to CF conventions

### 5.3 The Over-Spatialization Concern

**My concern (Michael's hypothesis, needs testing):**

GeoZarr may be trying to encode ALL the coordinate reference models that GeoTIFF supports (affine, GCPs, RPCs, geolocation arrays) within the format specification itself. This could be problematic because:

1. **Transformation is runtime logic, not format metadata**: The warper API decides *how* to transform, not the format
2. **Reinventing GDAL in a spec**: Format specs shouldn't duplicate transformation implementations
3. **Interoperability risk**: Different readers may implement different transformation methods
4. **Conflating two different worlds**: The CF/climate world and the remote sensing world have fundamentally different mental models (see Section 4.1)

**The two-worlds problem:**

GeoZarr sits at the intersection of two communities with different assumptions:

| Aspect | Remote Sensing World | Climate/CF World |
|--------|---------------------|------------------|
| Native state | Unrectified image | Model grid output |
| Coordinates describe | Where to transform TO | Where things ARE |
| Transformation | Correction (rectification) | Analysis choice (regridding) |
| "Success" | Producing a map-projected product | Using data on its native grid |

Trying to generalize across both by adding more transformation machinery to the spec conflates these use cases. The CF world doesn't need GCP support—they have coordinate arrays that completely describe their grids. The remote sensing world doesn't need 2D auxiliary coordinates for rectified products—they have geotransforms and RPCs.

**What GeoZarr *should* do (hypothesis):**
- Pass through CRS information cleanly (WKT2, PROJJSON)
- Pass through coordinate arrays or geotransform metadata
- Let reading libraries (GDAL, rioxarray, pyresample) handle the transformation logic
- Acknowledge two conformance classes with different assumptions rather than trying to unify them

**What GeoZarr should NOT do:**
- Define interpolation methods for non-affine transformations
- Specify what readers should do when they encounter GCPs or RPCs  
- Try to standardize warping semantics

The GeoTIFF spec got this right by explicitly putting interpolation "out of scope." GeoZarr should follow suit.

**⚠️ NEEDS VERIFICATION:** Is this a fair characterization of GeoZarr's direction? Need to review current spec drafts and SWG discussions more carefully.

---

## 6. The Spectrum of Coordinate Reference Models

### 6.1 A Taxonomy

From simplest to most complex:

1. **Bbox** - `xmin, ymin, xmax, ymax` + image dimensions → assumes axis-aligned, north-up
2. **Offset + scale** - origin point + pixel size → single tiepoint + scale equivalent
3. **6-parameter affine** - handles rotation and shear
4. **GCPs + polynomial** - sparse points, fitted transformation
5. **GCPs + TPS** - sparse points, exact interpolation
6. **RPCs** - rational polynomial sensor model
7. **Geolocation arrays** - full or subsampled lat/lon per pixel

### 6.2 What Each Format Actually Supports

| Format | Affine | Rotated Affine | GCPs | RPCs | GeoLoc Arrays |
|--------|--------|----------------|------|------|---------------|
| GeoTIFF | ✓ (baseline) | ✓ (matrix) | via multiple tiepoints (not baseline) | via metadata | ✗ |
| World file | ✓ | ✓ | ✗ | ✗ | ✗ |
| netCDF/CF | via grid_mapping | via grid_mapping | ✗ | ✗ | ✓ (2D coords) |
| HDF-EOS | ✓ | ✓ | ✓ | ✗ | ✓ |
| GDAL VRT | ✓ | ✓ | ✓ | ✓ | ✓ |

### 6.3 The Missing Middle Ground

Neither GeoTIFF nor Zarr/CF elegantly handles:
- Sparse but dense GCPs (more than a few, less than full geolocation)
- Mixed models (e.g., RPCs refined by GCPs)
- Progressive refinement

**GDAL handles all of these at runtime** but the metadata storage is fragmented.

---

## 7. Conclusions and Recommendations

### 7.1 What We Can Say With Confidence

1. **Single tiepoint + scale ≡ affine geotransform** - This is just math.

2. **Multiple tiepoints were never standardized** - The GeoTIFF spec explicitly called this "not Baseline" and recommended against interchange. The interpolation method was left undefined.

3. **GCPs and tiepoints diverged** - They solve different problems (post-hoc correction vs intrinsic georeferencing) and were handled separately in implementations.

4. **CF conventions came from a different world** - Climate modeling didn't need projected coordinates the way remote sensing did. The `grid_mapping` approach is workable but awkward.

5. **GDAL unified the practical landscape** - The transformer API handles all mechanisms through a common interface, even if the metadata storage remains fragmented.

### 7.2 Open Questions

1. **Were multiple tiepoints ever seriously used?** Need to check GDAL history and old implementations.

2. **What was the original intent for tiepoints vs GCPs?** Were they ever meant to converge?

3. **Is GeoZarr over-spatializing?** Need to review current spec and SWG discussions.

4. **Should GeoZarr just pass through CRS + coordinate metadata?** And leave transformation to runtime libraries?

### 7.3 For the Essay

**Potential framing:** "The mess we're in" is the result of:
- GeoTIFF standardizing early around the simplest case
- CF conventions emerging from a different domain with different priorities
- Format specs trying to encode runtime behavior
- No single authority for "how to georeference raster data"

**GDAL's role:** The de facto standard that papers over format differences, but this creates implicit dependencies that format specs don't capture.

**The way forward:** Format specs should be minimal (store CRS and coordinate metadata cleanly) and leave transformation semantics to libraries. GeoZarr should resist the temptation to re-specify transformation methods.

---

## Appendix A: Two Worlds, Two Mental Models

This appendix expands on the philosophical distinction between the remote sensing and climate/modeling communities that underlies much of the tension in geospatial format design.

### A.1 The Remote Sensing Mental Model

In the remote sensing world (Landsat, Sentinel, commercial satellites, aerial imagery), the fundamental workflow is:

1. **Acquire** raw imagery from a sensor
2. **Correct** geometric distortions (terrain, sensor geometry, timing)
3. **Rectify** to a standard map projection
4. **Distribute** the corrected product

The *goal* is a map-projected, geometrically corrected image. Raw or unrectified data is an intermediate state—useful for some applications, but not the end product. The geotransform describes the *corrected* relationship between pixels and coordinates.

When remote sensing practitioners encounter unrectified data with GCPs or RPCs, they think: "I need to warp this to get a usable product." The transformation is a *correction* to be applied.

### A.2 The Climate/Modeling Mental Model

In the climate and ocean modeling world (NEMO, MOM, WRF, CESM), the workflow is fundamentally different:

1. **Run** a model on a computational grid
2. **Output** results on that grid
3. **Analyze** the output, possibly regridding for comparison with other datasets

The model grid IS the coordinate system. When a climate modeler has NEMO output on a rotated pole ORCA grid with 2D `nav_lat`/`nav_lon` arrays, they're not thinking "this needs rectification." The curvilinear grid is the *authoritative* representation of where the data lives.

Regridding to a regular lat/lon grid is an *analytical choice*—necessary for some comparisons, but fundamentally a transformation that loses information about the native grid structure. It's not a correction; it's a deliberate resampling with explicit method choices (bilinear, conservative, nearest-neighbor).

### A.3 The Coordinate Array as Complete Description

This is why CF conventions treat 2D auxiliary coordinate variables as a complete georeferencing solution:

```
float temperature(time, y, x);
    temperature:coordinates = "lat lon";
    temperature:grid_mapping = "crs";

float lat(y, x);
float lon(y, x);
```

The `lat` and `lon` arrays aren't a deficient form of georeferencing waiting to be converted to a geotransform. They're the *complete* description of where each grid cell is located. The curvilinear structure is intentional—it's how the model discretized the domain.

This explains why:
- CDO and xarray regridding are explicit operations with user-specified methods
- There's no "auto-warp on read" in the netCDF ecosystem
- `pyresample` exists as a separate tool for deliberate resampling
- Nobody in the CF world talks about "rectification"

In effect, CF users are *always* in what the remote sensing world would call "coordinate space"—they work with the data on its native grid and transform only when analytically necessary.

### A.4 Implications for Format Design

These different mental models create tension when designing a format like GeoZarr that aims to serve both communities:

**Remote sensing expectations:**
- Geotransform or transformation metadata that can produce a map
- Clear path to "rectified" state
- Warp-on-read as a reasonable default

**Climate/modeling expectations:**
- Coordinate arrays that fully describe the grid
- No implicit transformation on read
- Regridding as an explicit analytical step

A format spec that tries to unify these by adding transformation machinery (polynomial methods, interpolation specifications) is likely to satisfy neither community. Better to:

1. Support both models cleanly (geotransform OR coordinate arrays)
2. Leave transformation semantics to libraries
3. Not try to specify what readers should do with non-affine georeferencing

The GeoTIFF spec's wisdom—"interpolation methods are not in scope"—applies equally well to GeoZarr.

### A.5 The "Carried" vs "Transformed" Distinction

GDAL's architecture embodies this distinction:

| Operation | Remote Sensing View | Climate/Modeling View |
|-----------|--------------------|-----------------------|
| **RasterIO** | Read pixels, carry metadata | Read data on native grid |
| **Warp** | Apply correction to get map | Deliberate resampling choice |

Both communities use RasterIO-style operations constantly. The difference is in how they view Warp:
- Remote sensing: essential correction step
- Climate/modeling: analytical transformation to be used judiciously

Format specs live in RasterIO space—they describe what's stored, not how to transform it. When GeoZarr tries to move into Warp space by specifying transformation methods, it's overstepping the natural boundary between format and library.

---

## Appendix B: Sources and References

### Primary Sources
- OGC GeoTIFF Standard 1.1: https://docs.ogc.org/is/19-008r4/19-008r4.html
- GeoTIFF Spec 1.0 (1995): http://geotiff.maptools.org/spec/geotiffhome.html
- GDAL RFC 4 (Geolocation Arrays): https://gdal.org/development/rfc/rfc4_geolocate.html
- GDAL RFC 22 (RPC): https://gdal.org/development/rfc/rfc22_rpc.html
- CF Conventions (grid mappings): https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/ch05s06.html
- GeoZarr Charter: https://github.com/zarr-developers/geozarr-spec/blob/main/CHARTER.adoc
- CF Issue #222 (CRS WKT): https://github.com/cf-convention/cf-conventions/issues/222

### Secondary Sources
- GDAL Geotransform Tutorial: https://gdal.org/tutorials/geotransforms_tut.html
- GDAL Zarr Driver: https://gdal.org/drivers/raster/zarr.html
- pyproj CF documentation: https://pyproj4.github.io/pyproj/stable/build_crs_cf.html

### To Consult
- GDAL devs and community
- GeoZarr SWG members (current direction)
- xarray/rioxarray developers (practical CF/CRS handling)

---

## Revision History

- 2026-02-04: Initial draft compiled from web research
