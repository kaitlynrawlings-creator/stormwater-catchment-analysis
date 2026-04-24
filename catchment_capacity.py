"""
catchment_capacity.py
Phase 3 — Pipe Capacity Analysis

Joins each catchment (watersheds_p_Dissolve) to its connected enclosed gravity
pipes in the storm conduits layer via gridcode → structure number → upstream/
downstream structure fields. Applies Manning's full-pipe equation to each pipe
and sums capacity per catchment.

Inputs (edit CONFIGURATION block below):
  - watersheds_p_Dissolve : catchment polygons from Phase 2
  - storm conduits layer  : conduit feature class from your asset management GDB
  - upstream structure field (CONDUIT_US_STRUCT) : upstream connected structure number
  - downstream structure field (CONDUIT_DS_STRUCT): downstream connected structure number

Output:
  - watersheds_capacity   : copy of watersheds_p_Dissolve with capacity fields added

Output fields added:
  PIPE_COUNT    - number of enclosed gravity pipes connected to this catchment's inlet
  TOTAL_CAP_AF  - total full-flow pipe capacity, acre-feet per 24-hour storm
  TOTAL_CAP_GL  - total full-flow pipe capacity, gallons
  DOM_DIA_IN    - dominant (most common) pipe diameter in inches
  AVG_SLOPE_PCT - average slope across connected pipes, percent
  CAP_PER_ACRE  - TOTAL_CAP_AF divided by catchment area in acres
                  (only computed for catchments >= MIN_CATCHMENT_AREA_M2)

Manning's equation: Q = (1/n) * A * R^(2/3) * S^(1/2)
  Units: Q in cfs, A in ft², R in ft, S dimensionless (decimal)
  Assumes full-pipe flow (theoretical maximum capacity)
  Applies to enclosed gravity circular pipes only (CONDUIT_TYPE_ENCLOSED = 1)

NOTE: Field names in the CONFIGURATION block below are placeholders. Update them
to match your asset management system's schema before running.
"""

import arcpy
import os
import math
from collections import Counter

# ── CONFIGURATION ────────────────────────────────────────────────────────────

GDB = r"C:\path\to\your\geodatabase.gdb"  # update to your GDB path

CATCHMENTS_IN  = os.path.join(GDB, "watersheds_p_dissolve")
CONDUITS       = os.path.join(GDB, "StormConduits")     # your conduit feature class
CATCHMENTS_OUT = os.path.join(GDB, "watersheds_capacity")

# ── FIELD NAME MAPPING ───────────────────────────────────────────────────────
# Update these to match your asset management system's field names

CONDUIT_TYPE_FIELD    = "CONDUIT_TYPE_CD"    # pipe type code field
CONDUIT_TYPE_ENCLOSED = 1                    # type code value for enclosed gravity pipe
CONDUIT_US_STRUCT     = "CONDUIT_US_STR"     # upstream structure number field
CONDUIT_DS_STRUCT     = "CONDUIT_DS_STR"     # downstream structure number field
CONDUIT_DIA_FIELD     = "CONDUIT_DIA"        # pipe diameter field (inches)
CONDUIT_SLOPE_FIELD   = "CONDUIT_SLOPE"      # pipe slope field (percent)
CONDUIT_MATR_FIELD    = "CONDUIT_MATR_CD"    # pipe material code field

# ── ANALYSIS PARAMETERS ──────────────────────────────────────────────────────

# Minimum catchment area in square meters to compute CAP_PER_ACRE
# Excludes raster artifact slivers from the normalized metric
MIN_CATCHMENT_AREA_M2 = 2000

# Default slope (percent) applied where slope field is null or zero
# Use your local stormwater design standard minimum as the fallback value
DEFAULT_SLOPE_PCT = 1.0  # update to match your jurisdiction's design standard

# Manning's n by material code
# Keys correspond to your asset management system's pipe material codes
# Update codes and n values to match your data
MANNINGS_N = {
    1:  0.013,   # Concrete
    2:  0.013,   # Reinforced Concrete
    3:  0.012,   # PVC
    4:  0.011,   # HDPE
    5:  0.024,   # Corrugated Metal
    6:  0.015,   # Clay
    7:  0.013,   # Cast Iron
    8:  0.014,   # Ductile Iron
    9:  0.013,   # Other
    10: 0.025,   # Brick
}
DEFAULT_N = 0.013  # fallback where material code is not in table

# ── CONSTANTS ────────────────────────────────────────────────────────────────

CFS_TO_AF_24HR = 3600 * 24 / 43560   # cfs → acre-feet per 24 hours
AF_TO_GAL      = 325851              # acre-feet → gallons

# ── FUNCTIONS ────────────────────────────────────────────────────────────────

def mannings_cfs(dia_in, slope_pct, n):
    """
    Full-pipe Manning's Q in cfs for a circular pipe.
    dia_in   : pipe diameter in inches
    slope_pct: pipe slope in percent (e.g. 1.0 = 1%)
    n        : Manning's roughness coefficient
    Returns Q in cfs, or 0.0 if inputs are invalid.
    """
    if dia_in <= 0 or slope_pct <= 0 or n <= 0:
        return 0.0
    d  = dia_in / 12.0           # inches → feet
    r  = d / 4.0                 # hydraulic radius for full circular pipe = D/4
    a  = math.pi * (d / 2.0)**2  # cross-sectional area
    s  = slope_pct / 100.0       # percent → decimal
    return (1.0 / n) * a * (r ** (2.0 / 3.0)) * (s ** 0.5)

def cfs_to_af_24hr(q_cfs):
    """Convert cfs (instantaneous full-flow) to acre-feet over 24 hours."""
    return q_cfs * CFS_TO_AF_24HR

# ── MAIN ─────────────────────────────────────────────────────────────────────

arcpy.env.overwriteOutput = True
arcpy.env.workspace = GDB

print("Phase 3 — Pipe Capacity Analysis")
print(f"  Catchments : {CATCHMENTS_IN}")
print(f"  Conduits   : {CONDUITS}")
print(f"  Output     : {CATCHMENTS_OUT}")

# Step 1: Copy catchments to output
# ────────────────────────────────
print("\nStep 1: Copying catchments to output...")
arcpy.management.CopyFeatures(CATCHMENTS_IN, CATCHMENTS_OUT)
print(f"  Copied {arcpy.management.GetCount(CATCHMENTS_OUT)[0]} catchments")

# Step 2: Add output fields
# ─────────────────────────
print("\nStep 2: Adding output fields...")
fields_to_add = [
    ("PIPE_COUNT",    "SHORT",  "Pipe count"),
    ("TOTAL_CAP_AF",  "DOUBLE", "Total capacity acre-feet"),
    ("TOTAL_CAP_GL",  "DOUBLE", "Total capacity gallons"),
    ("DOM_DIA_IN",    "SHORT",  "Dominant diameter inches"),
    ("AVG_SLOPE_PCT", "DOUBLE", "Average slope percent"),
    ("CAP_PER_ACRE",  "DOUBLE", "Capacity per acre"),
]
existing = [f.name.upper() for f in arcpy.ListFields(CATCHMENTS_OUT)]
for fname, ftype, falias in fields_to_add:
    if fname not in existing:
        arcpy.management.AddField(CATCHMENTS_OUT, fname, ftype, field_alias=falias)
        print(f"  Added field: {fname}")
    else:
        print(f"  Field already exists, skipping: {fname}")

# Step 3: Load conduit data into memory
# ─────────────────────────────────────
# Filter to enclosed gravity pipes only
# Manning's circular pipe equation applies to closed circular sections only
print("\nStep 3: Loading enclosed gravity conduit data...")

conduit_fields = [
    CONDUIT_US_STRUCT, CONDUIT_DS_STRUCT, CONDUIT_TYPE_FIELD,
    CONDUIT_DIA_FIELD, CONDUIT_SLOPE_FIELD, CONDUIT_MATR_FIELD
]

# Build lookup: structure_number → list of (dia_in, slope_pct, n)
# Each structure may have multiple pipes connected as upstream or downstream node
struct_pipes = {}

pipe_eg      = 0
pipe_no_dia  = 0
pipe_def_slp = 0

with arcpy.da.SearchCursor(
    CONDUITS,
    conduit_fields,
    where_clause=f"{CONDUIT_TYPE_FIELD} = {CONDUIT_TYPE_ENCLOSED}"
) as cur:
    for us_str, ds_str, type_cd, dia, slope, matr in cur:
        pipe_eg += 1

        # Skip pipes with no diameter — cannot compute capacity
        if dia is None or dia <= 0:
            pipe_no_dia += 1
            continue

        # Apply default slope where null or zero
        if slope is None or slope <= 0:
            slope = DEFAULT_SLOPE_PCT
            pipe_def_slp += 1

        # Manning's n from material code
        n = MANNINGS_N.get(int(matr) if matr else 0, DEFAULT_N)

        pipe_tuple = (float(dia), float(slope), n)

        # Index by both upstream and downstream structure numbers
        for struct_no in [us_str, ds_str]:
            if struct_no is not None and str(struct_no).strip() != "":
                key = str(struct_no).strip()
                if key not in struct_pipes:
                    struct_pipes[key] = []
                struct_pipes[key].append(pipe_tuple)

pipe_total = int(arcpy.management.GetCount(CONDUITS)[0])
print(f"  Total conduits in layer  : {pipe_total}")
print(f"  Enclosed gravity pipes   : {pipe_eg}")
print(f"  Skipped (null/zero dia)  : {pipe_no_dia}")
print(f"  Default slope applied    : {pipe_def_slp} ({pipe_def_slp/max(pipe_eg,1)*100:.0f}%)")
print(f"  Unique structure keys    : {len(struct_pipes)}")

# Step 4: Compute capacity per catchment
# ───────────────────────────────────────
print("\nStep 4: Computing capacity per catchment...")

update_fields = [
    "gridcode", "Shape_Area",
    "PIPE_COUNT", "TOTAL_CAP_AF", "TOTAL_CAP_GL",
    "DOM_DIA_IN", "AVG_SLOPE_PCT", "CAP_PER_ACRE"
]

stats = {"total": 0, "with_pipes": 0, "no_pipes": 0}

with arcpy.da.UpdateCursor(CATCHMENTS_OUT, update_fields) as cur:
    for row in cur:
        gridcode   = row[0]
        shape_area = row[1]  # square meters (assumes projected CRS)
        stats["total"] += 1

        key   = str(gridcode).strip() if gridcode is not None else ""
        pipes = struct_pipes.get(key, [])

        if not pipes:
            row[2] = 0      # PIPE_COUNT
            row[3] = 0.0    # TOTAL_CAP_AF
            row[4] = 0.0    # TOTAL_CAP_GL
            row[5] = 0      # DOM_DIA_IN
            row[6] = 0.0    # AVG_SLOPE_PCT
            row[7] = None   # CAP_PER_ACRE
            stats["no_pipes"] += 1
            cur.updateRow(row)
            continue

        stats["with_pipes"] += 1

        # Compute Manning's Q for each connected pipe
        total_q_cfs = 0.0
        diameters   = []
        slopes      = []

        for dia_in, slope_pct, n in pipes:
            q = mannings_cfs(dia_in, slope_pct, n)
            total_q_cfs += q
            diameters.append(int(round(dia_in)))
            slopes.append(slope_pct)

        total_af = cfs_to_af_24hr(total_q_cfs)
        total_gl = total_af * AF_TO_GAL

        # Dominant diameter = most common diameter value
        dom_dia = Counter(diameters).most_common(1)[0][0] if diameters else 0

        # Average slope
        avg_slope = sum(slopes) / len(slopes) if slopes else 0.0

        # Capacity per acre — only for catchments above minimum area threshold
        area_m2    = shape_area if shape_area else 0.0
        area_acres = area_m2 / 4046.856  # m² → acres
        if area_m2 >= MIN_CATCHMENT_AREA_M2 and area_acres > 0:
            cap_per_acre = total_af / area_acres
        else:
            cap_per_acre = None

        row[2] = len(pipes)   # PIPE_COUNT
        row[3] = total_af     # TOTAL_CAP_AF
        row[4] = total_gl     # TOTAL_CAP_GL
        row[5] = dom_dia      # DOM_DIA_IN
        row[6] = avg_slope    # AVG_SLOPE_PCT
        row[7] = cap_per_acre # CAP_PER_ACRE

        cur.updateRow(row)

# Step 5: Print summary
# ──────────────────────
print("\n── SUMMARY ──────────────────────────────────────")
print(f"  Total catchments  : {stats['total']}")
print(f"  With pipe data    : {stats['with_pipes']} ({stats['with_pipes']/stats['total']*100:.1f}%)")
print(f"  No pipe data      : {stats['no_pipes']} ({stats['no_pipes']/stats['total']*100:.1f}%)")
print(f"\nOutput written to: {CATCHMENTS_OUT}")
print("\nPhase 3 complete. Run catchment_runoff.py next.")
