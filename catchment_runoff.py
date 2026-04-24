"""
catchment_runoff.py
Phase 4 — Runoff, Risk Assessment, Open Channel Inventory, BMP Integration,
           and Depression Depth Ponding

Runs after catchment_capacity.py. Adds the following fields to watersheds_capacity:

  RUNOFF FIELDS (NRCS TR-55):
    AVG_CN         - area-weighted average curve number per catchment
    RUNOFF_IN      - runoff depth in inches (TR-55 equation)
    RUNOFF_AF      - runoff volume in acre-feet
    CAP_DEFICIT_AF - runoff volume minus pipe capacity (positive = deficit)
    RISK_FLAG      - pipe-only risk assessment
                     1  = at-risk (runoff > pipe capacity)
                     0  = adequate
                    -1  = no pipe data

  OPEN CHANNEL INVENTORY FIELDS (spatial intersect):
    DITCH_FT       - linear feet of ditch within catchment
    SWALE_FT       - linear feet of swale within catchment
    STREAM_FT      - linear feet of stream channel within catchment
    CONCCH_FT      - linear feet of concrete channel within catchment
    TRENCH_FT      - linear feet of trench drain within catchment
    OPEN_CH_FT     - total open channel footage (all types combined)

  BMP FIELDS (spatial join):
    BMP_PRESENT    - 1 if active flood control BMP intersects catchment, else 0
    BMP_TYPE       - comma-separated list of BMP type names present

  DEPRESSION DEPTH FIELDS (zonal statistics):
    MAX_POND_IN    - maximum depression depth in inches within catchment
    AVG_POND_IN    - average depression depth in inches within catchment

  ADJUSTED RISK FLAG:
    RISK_FLAG_ADJ  - risk flag accounting for supplemental infrastructure
                     1  = at-risk, no supplemental infrastructure present
                     2  = at-risk, open channel conveyance present — field verify
                     3  = at-risk, flood control BMP present — field verify
                     4  = at-risk, both open channel and BMP present — field verify
                     0  = adequate pipe capacity
                    -1  = no pipe data, no open channel or BMP present
                    -2  = no pipe data, but open channel or BMP present

NOTE: Field names in the CONFIGURATION block below are placeholders. Update them
to match your asset management system's schema before running. Design storm
parameters should reflect local NOAA Atlas 14 point frequency estimates for your
study area.
"""

import arcpy
import os
import math

# ── CONFIGURATION ────────────────────────────────────────────────────────────

GDB = r"C:\path\to\your\geodatabase.gdb"  # update to your GDB path

CATCHMENTS     = os.path.join(GDB, "watersheds_capacity")
CONDUITS       = os.path.join(GDB, "StormConduits")       # your conduit feature class
BMP_ASSETS     = os.path.join(GDB, "BMPAssets")           # your BMP asset feature class
NLCD           = r"C:\path\to\NLCD_land_cover.tif"        # NLCD annual land cover raster
SSURGO         = r"C:\path\to\SSURGO_hydrologic_soils.tif"# SSURGO hydrologic soil groups raster
DEPRESSION_DEM = os.path.join(GDB, "depression_depth")    # depression depth raster from Phase 1

# CN raster saved to main GDB — required for ZonalStatisticsAsTable
CN_RASTER = os.path.join(GDB, "cn_raster")

# ── FIELD NAME MAPPING ───────────────────────────────────────────────────────
# Update these to match your asset management system's field names

CONDUIT_TYPE_FIELD = "CONDUIT_TYPE_CD"  # pipe type code field
CONDUIT_LENGTH_FIELD = "CONDUIT_LENGTH" # pipe length field (feet)

# Open channel type codes and their corresponding output field names
# Update type code keys to match your asset management system
OPEN_CHANNEL_TYPES = {
    2: "STREAM_FT",   # stream channel
    4: "TRENCH_FT",   # trench drain
    5: "DITCH_FT",    # ditch
    6: "CONCCH_FT",   # concrete channel
    9: "SWALE_FT",    # swale
}

# BMP field names
BMP_STATUS_FIELD   = "BMP_STATUS"       # operational status field
BMP_STATUS_ACTIVE  = 1                  # value indicating operational/active status
BMP_TYPE_FIELD     = "BMP_TYPE_CD"      # BMP type code field

# BMP type codes treated as flood control detention
# Update to match your asset management system's type codes
FLOOD_CONTROL_BMP_TYPES = [1, 2, 13, 23, 24, 26]

# BMP type code → display name mapping
# Update codes and names to match your asset management system
BMP_TYPE_NAMES = {
    1: "Dry Basin", 2: "Wet Bottom", 3: "Concrete", 4: "Underground Vault",
    5: "Parking Lot", 6: "Roof", 7: "Native Vegetation", 8: "Rain Garden",
    9: "Infiltration Practice", 10: "Bioretention", 11: "Pervious Pavement",
    12: "Wetland Swale", 13: "Ext Det Wetland", 14: "Bioswale",
    15: "Media Filter", 16: "Native Veg Swale", 17: "Turf Grass Swale",
    18: "Storm Treat In Line", 19: "Forebay", 20: "Stream Corridor",
    21: "Sand Filter", 22: "Storm Treat Off Line", 23: "Dry Basin Backup",
    24: "Restrictor Plate", 25: "Native Restoration", 26: "Ext Wet Detention",
    27: "Sediment Trap", 28: "Sand Oil Separator"
}

# NLCD and SSURGO raster field names (produced by the Combine tool)
# These will reflect the filenames of your input rasters — update accordingly
NLCD_FIELD   = "NLCD_land_cover"           # NLCD value field in combined raster
SSURGO_FIELD = "SSURGO_hydrologic_soils"   # SSURGO value field in combined raster

# ── DESIGN STORM PARAMETERS ──────────────────────────────────────────────────

# Update DESIGN_STORM_P to the NOAA Atlas 14 point frequency estimate
# for your study area and desired return interval / storm duration
DESIGN_STORM_P = 5.50  # inches — update to your local NOAA Atlas 14 value
                        # (10-year 24-hour depth; see https://hdsc.nws.noaa.gov/hdsc/pfds/)

# ── CN LOOKUP TABLE ──────────────────────────────────────────────────────────
# Key: (NLCD_code, soil_value) → CN
# Soil values: 2=Group B, 3=Group C, 4=Group D, 7=Group C/D, 15=NoData fallback
# Source: USDA NRCS TR-55, Table 2-2 (USDA NRCS, 1986)

CN_TABLE = {
    (11, 2): 100, (11, 3): 100, (11, 4): 100, (11, 7): 100, (11, 15): 100,
    (21, 2): 61,  (21, 3): 74,  (21, 4): 80,  (21, 7): 77,  (21, 15): 77,
    (22, 2): 77,  (22, 3): 85,  (22, 4): 90,  (22, 7): 88,  (22, 15): 88,
    (23, 2): 85,  (23, 3): 90,  (23, 4): 92,  (23, 7): 91,  (23, 15): 91,
    (24, 2): 90,  (24, 3): 92,  (24, 4): 93,  (24, 7): 93,  (24, 15): 93,
    (31, 2): 77,  (31, 3): 86,  (31, 4): 91,  (31, 7): 89,  (31, 15): 89,
    (41, 2): 55,  (41, 3): 70,  (41, 4): 77,  (41, 7): 74,  (41, 15): 74,
    (42, 2): 55,  (42, 3): 70,  (42, 4): 77,  (42, 7): 74,  (42, 15): 74,
    (43, 2): 55,  (43, 3): 70,  (43, 4): 77,  (43, 7): 74,  (43, 15): 74,
    (52, 2): 65,  (52, 3): 76,  (52, 4): 82,  (52, 7): 79,  (52, 15): 79,
    (71, 2): 61,  (71, 3): 74,  (71, 4): 80,  (71, 7): 77,  (71, 15): 77,
    (81, 2): 69,  (81, 3): 79,  (81, 4): 84,  (81, 7): 82,  (81, 15): 82,
    (82, 2): 75,  (82, 3): 82,  (82, 4): 86,  (82, 7): 84,  (82, 15): 84,
    (90, 2): 78,  (90, 3): 83,  (90, 4): 87,  (90, 7): 85,  (90, 15): 85,
    (95, 2): 78,  (95, 3): 83,  (95, 4): 87,  (95, 7): 85,  (95, 15): 85,
}

# ── FUNCTIONS ────────────────────────────────────────────────────────────────

def tr55_runoff(p_in, cn):
    """
    NRCS TR-55 runoff equation.
    p_in : rainfall depth in inches
    cn   : curve number
    Returns runoff depth Q in inches, or 0.0 if rainfall <= initial abstraction.
    Source: USDA NRCS (1986), Urban Hydrology for Small Watersheds, TR-55.
    """
    if cn <= 0 or cn > 100:
        return 0.0
    s  = (1000.0 / cn) - 10.0   # potential maximum soil retention
    ia = 0.2 * s                 # initial abstraction
    if p_in <= ia:
        return 0.0
    return ((p_in - ia) ** 2) / (p_in + 0.8 * s)

def runoff_af(runoff_in, area_m2):
    """Convert runoff depth (inches) and catchment area (m²) to acre-feet."""
    area_acres = area_m2 / 4046.856
    return (runoff_in / 12.0) * area_acres

# ── MAIN ─────────────────────────────────────────────────────────────────────

arcpy.env.overwriteOutput = True
arcpy.env.workspace = GDB
scratch = arcpy.env.scratchGDB

print("Phase 4 — Runoff, Risk, Open Channel Inventory, BMP, Depression Depth")
print(f"  Catchments : {CATCHMENTS}")

# ════════════════════════════════════════════════════════════════════
# BLOCK 1: NRCS TR-55 RUNOFF AND RISK FLAG
# ════════════════════════════════════════════════════════════════════
print("\n── Block 1: NRCS TR-55 Runoff ──────────────────────────────────")

# Step 1a: Build CN raster via Combine
print("Step 1a: Building CN raster from NLCD and SSURGO...")
combine_raster = os.path.join(scratch, "nlcd_ssurgo_combine")
arcpy.env.extent    = CATCHMENTS
arcpy.env.snapRaster = NLCD
arcpy.sa.Combine([NLCD, SSURGO]).save(combine_raster)
print("  Combine complete")

# Step 1b: Map combined raster values to CN
print("Step 1b: Mapping combined pixel values to curve numbers...")

value_to_cn = {}
with arcpy.da.SearchCursor(combine_raster, ["Value", NLCD_FIELD, SSURGO_FIELD]) as cur:
    for val, nlcd_code, soil_code in cur:
        cn = CN_TABLE.get((int(nlcd_code), int(soil_code)), None)
        if cn is None:
            # Fallback: use soil code 15 (NoData fallback CN values)
            cn = CN_TABLE.get((int(nlcd_code), 15), 75)
        value_to_cn[int(val)] = cn

print(f"  Mapped {len(value_to_cn)} unique raster combinations to CN values")

# Step 1c: Reclassify combine raster to CN values
print("Step 1c: Reclassifying to CN raster...")
remap_list = [[k, v] for k, v in value_to_cn.items()]
remap = arcpy.sa.RemapValue(remap_list)
cn_ras = arcpy.sa.Reclassify(combine_raster, "Value", remap, "NODATA")
cn_ras.save(CN_RASTER)
print(f"  CN raster saved to: {CN_RASTER}")

# Step 1d: Zonal Statistics — average CN per catchment
print("Step 1d: Zonal Statistics — average CN per catchment...")
cn_zonal_table = os.path.join(scratch, "cn_zonal_table")
arcpy.sa.ZonalStatisticsAsTable(
    CATCHMENTS, "gridcode",
    cn_ras,
    cn_zonal_table,
    "DATA", "MEAN"
)
print("  Zonal Statistics complete")

gridcode_to_cn = {}
with arcpy.da.SearchCursor(cn_zonal_table, ["gridcode", "MEAN"]) as cur:
    for gc, mean_cn in cur:
        gridcode_to_cn[gc] = mean_cn
print(f"  CN values loaded for {len(gridcode_to_cn)} catchments")

# Step 1e: Add runoff and risk fields
print("Step 1e: Adding runoff and risk fields...")
runoff_fields = [
    ("AVG_CN",         "DOUBLE", "Average curve number"),
    ("RUNOFF_IN",      "DOUBLE", "Runoff depth inches"),
    ("RUNOFF_AF",      "DOUBLE", "Runoff volume acre-feet"),
    ("CAP_DEFICIT_AF", "DOUBLE", "Capacity deficit acre-feet"),
    ("RISK_FLAG",      "SHORT",  "Risk flag pipe only"),
]
existing = [f.name.upper() for f in arcpy.ListFields(CATCHMENTS)]
for fname, ftype, falias in runoff_fields:
    if fname.upper() not in existing:
        arcpy.management.AddField(CATCHMENTS, fname, ftype, field_alias=falias)

# Step 1f: Calculate runoff and risk per catchment
print("Step 1f: Calculating runoff and risk flag...")
update_fields = [
    "gridcode", "Shape_Area", "TOTAL_CAP_AF",
    "AVG_CN", "RUNOFF_IN", "RUNOFF_AF",
    "CAP_DEFICIT_AF", "RISK_FLAG"
]

risk_stats = {"at_risk": 0, "adequate": 0, "no_data": 0}

with arcpy.da.UpdateCursor(CATCHMENTS, update_fields) as cur:
    for row in cur:
        gridcode  = row[0]
        area_m2   = row[1] if row[1] else 0.0
        total_cap = row[2] if row[2] else 0.0

        avg_cn = gridcode_to_cn.get(gridcode, None)

        if avg_cn is None or avg_cn <= 0:
            row[3] = None
            row[4] = None
            row[5] = None
            row[6] = None
            row[7] = -1       # RISK_FLAG — no CN data
            risk_stats["no_data"] += 1
            cur.updateRow(row)
            continue

        q_in    = tr55_runoff(DESIGN_STORM_P, avg_cn)
        q_af    = runoff_af(q_in, area_m2)
        deficit = q_af - total_cap

        if total_cap == 0.0:
            risk_flag = -1
            risk_stats["no_data"] += 1
        elif q_af > total_cap:
            risk_flag = 1
            risk_stats["at_risk"] += 1
        else:
            risk_flag = 0
            risk_stats["adequate"] += 1

        row[3] = avg_cn
        row[4] = q_in
        row[5] = q_af
        row[6] = max(deficit, 0.0)
        row[7] = risk_flag
        cur.updateRow(row)

print(f"  At-risk (RISK_FLAG=1) : {risk_stats['at_risk']}")
print(f"  Adequate (RISK_FLAG=0): {risk_stats['adequate']}")
print(f"  No data (RISK_FLAG=-1): {risk_stats['no_data']}")

# ════════════════════════════════════════════════════════════════════
# BLOCK 2: OPEN CHANNEL INVENTORY
# ════════════════════════════════════════════════════════════════════
print("\n── Block 2: Open Channel Inventory ─────────────────────────────")

oc_fields = [
    ("DITCH_FT",   "DOUBLE", "Ditch linear feet"),
    ("SWALE_FT",   "DOUBLE", "Swale linear feet"),
    ("STREAM_FT",  "DOUBLE", "Stream channel linear feet"),
    ("CONCCH_FT",  "DOUBLE", "Concrete channel linear feet"),
    ("TRENCH_FT",  "DOUBLE", "Trench drain linear feet"),
    ("OPEN_CH_FT", "DOUBLE", "Total open channel linear feet"),
]
existing = [f.name.upper() for f in arcpy.ListFields(CATCHMENTS)]
for fname, ftype, falias in oc_fields:
    if fname.upper() not in existing:
        arcpy.management.AddField(CATCHMENTS, fname, ftype, field_alias=falias)

print("Step 2a: Initializing open channel fields to 0...")
oc_field_names = [f[0] for f in oc_fields]
with arcpy.da.UpdateCursor(CATCHMENTS, oc_field_names) as cur:
    for row in cur:
        cur.updateRow([0.0] * len(oc_field_names))

print("Step 2b: Intersecting open channel conduits with catchments...")
type_codes  = list(OPEN_CHANNEL_TYPES.keys())
type_filter = f"{CONDUIT_TYPE_FIELD} IN ({','.join(str(t) for t in type_codes)})"

open_ch_intersect = os.path.join(scratch, "open_ch_intersect")
arcpy.analysis.PairwiseIntersect(
    [CONDUITS, CATCHMENTS],
    open_ch_intersect,
    "ALL"
)
print("  Intersect complete")

print("Step 2c: Summing open channel footage by catchment and type...")
gc_oc = {}

with arcpy.da.SearchCursor(
    open_ch_intersect,
    ["gridcode", CONDUIT_TYPE_FIELD, CONDUIT_LENGTH_FIELD],
    where_clause=type_filter
) as cur:
    for gc, type_cd, length in cur:
        if gc is None or length is None:
            continue
        field = OPEN_CHANNEL_TYPES.get(int(type_cd))
        if field is None:
            continue
        if gc not in gc_oc:
            gc_oc[gc] = {f: 0.0 for f in OPEN_CHANNEL_TYPES.values()}
        gc_oc[gc][field] += float(length)

print(f"  Catchments with open channel footage: {len(gc_oc)}")

print("Step 2d: Writing open channel footage to catchments...")
update_fields = ["gridcode"] + [f[0] for f in oc_fields]
with arcpy.da.UpdateCursor(CATCHMENTS, update_fields) as cur:
    for row in cur:
        gc = row[0]
        if gc in gc_oc:
            data = gc_oc[gc]
            row[1] = data.get("DITCH_FT",   0.0)
            row[2] = data.get("SWALE_FT",   0.0)
            row[3] = data.get("STREAM_FT",  0.0)
            row[4] = data.get("CONCCH_FT",  0.0)
            row[5] = data.get("TRENCH_FT",  0.0)
            row[6] = sum(data.values())          # OPEN_CH_FT total
            cur.updateRow(row)

# ════════════════════════════════════════════════════════════════════
# BLOCK 3: BMP INTEGRATION
# ════════════════════════════════════════════════════════════════════
print("\n── Block 3: BMP Integration ─────────────────────────────────────")

bmp_fields = [
    ("BMP_PRESENT", "SHORT", "Flood control BMP present"),
    ("BMP_TYPE",    "TEXT",  "BMP type names present"),
]
existing = [f.name.upper() for f in arcpy.ListFields(CATCHMENTS)]
for fname, ftype, falias in bmp_fields:
    if fname.upper() not in existing:
        if ftype == "TEXT":
            arcpy.management.AddField(
                CATCHMENTS, fname, ftype,
                field_length=255, field_alias=falias
            )
        else:
            arcpy.management.AddField(CATCHMENTS, fname, ftype, field_alias=falias)

print("Step 3a: Initializing BMP fields...")
with arcpy.da.UpdateCursor(CATCHMENTS, ["BMP_PRESENT", "BMP_TYPE"]) as cur:
    for row in cur:
        cur.updateRow([0, ""])

print("Step 3b: Filtering BMPs to active flood control types...")
type_list  = ",".join(str(t) for t in FLOOD_CONTROL_BMP_TYPES)
bmp_filter = f"{BMP_STATUS_FIELD} = {BMP_STATUS_ACTIVE} AND {BMP_TYPE_FIELD} IN ({type_list})"

bmp_layer = "bmp_fc_layer"
arcpy.management.MakeFeatureLayer(BMP_ASSETS, bmp_layer, bmp_filter)
bmp_count = int(arcpy.management.GetCount(bmp_layer)[0])
print(f"  Active flood control BMPs: {bmp_count}")

print("Step 3c: Spatial join — catchments intersecting flood control BMPs...")
bmp_intersect = os.path.join(scratch, "bmp_intersect")
arcpy.analysis.PairwiseIntersect(
    [bmp_layer, CATCHMENTS],
    bmp_intersect,
    "ALL"
)

gc_bmp = {}
with arcpy.da.SearchCursor(bmp_intersect, ["gridcode", BMP_TYPE_FIELD]) as cur:
    for gc, bmp_type in cur:
        if gc is None:
            continue
        name = BMP_TYPE_NAMES.get(int(bmp_type), f"Type {bmp_type}")
        if gc not in gc_bmp:
            gc_bmp[gc] = set()
        gc_bmp[gc].add(name)

print(f"  Catchments intersecting a flood control BMP: {len(gc_bmp)}")

print("Step 3d: Writing BMP fields to catchments...")
with arcpy.da.UpdateCursor(CATCHMENTS, ["gridcode", "BMP_PRESENT", "BMP_TYPE"]) as cur:
    for row in cur:
        gc = row[0]
        if gc in gc_bmp:
            row[1] = 1
            row[2] = ", ".join(sorted(gc_bmp[gc]))
            cur.updateRow(row)

# ════════════════════════════════════════════════════════════════════
# BLOCK 4: DEPRESSION DEPTH — TERRAIN PONDING
# ════════════════════════════════════════════════════════════════════
print("\n── Block 4: Depression Depth Ponding ───────────────────────────")

pond_fields = [
    ("MAX_POND_IN", "DOUBLE", "Max depression depth inches"),
    ("AVG_POND_IN", "DOUBLE", "Avg depression depth inches"),
]
existing = [f.name.upper() for f in arcpy.ListFields(CATCHMENTS)]
for fname, ftype, falias in pond_fields:
    if fname.upper() not in existing:
        arcpy.management.AddField(CATCHMENTS, fname, ftype, field_alias=falias)

# Depression depth raster is assumed to be in meters (projected CRS)
# Converted to inches by multiplying by 39.3701
print("Step 4a: Zonal Statistics — depression depth per catchment...")
pond_zonal_max  = os.path.join(scratch, "pond_zonal_max")
pond_zonal_mean = os.path.join(scratch, "pond_zonal_mean")

arcpy.sa.ZonalStatisticsAsTable(
    CATCHMENTS, "gridcode", DEPRESSION_DEM,
    pond_zonal_max,  "DATA", "MAXIMUM"
)
arcpy.sa.ZonalStatisticsAsTable(
    CATCHMENTS, "gridcode", DEPRESSION_DEM,
    pond_zonal_mean, "DATA", "MEAN"
)
print("  Zonal Statistics complete")

gc_pond_max  = {}
gc_pond_mean = {}
with arcpy.da.SearchCursor(pond_zonal_max,  ["gridcode", "MAX"])  as cur:
    for gc, val in cur:
        gc_pond_max[gc] = val * 39.3701   # meters → inches
with arcpy.da.SearchCursor(pond_zonal_mean, ["gridcode", "MEAN"]) as cur:
    for gc, val in cur:
        gc_pond_mean[gc] = val * 39.3701

print(f"  Depression depth loaded for {len(gc_pond_max)} catchments")

print("Step 4b: Writing depression depth to catchments...")
with arcpy.da.UpdateCursor(
    CATCHMENTS, ["gridcode", "MAX_POND_IN", "AVG_POND_IN"]
) as cur:
    for row in cur:
        gc = row[0]
        row[1] = gc_pond_max.get(gc,  None)
        row[2] = gc_pond_mean.get(gc, None)
        cur.updateRow(row)

# ════════════════════════════════════════════════════════════════════
# BLOCK 5: ADJUSTED RISK FLAG
# ════════════════════════════════════════════════════════════════════
print("\n── Block 5: Adjusted Risk Flag ──────────────────────────────────")

existing = [f.name.upper() for f in arcpy.ListFields(CATCHMENTS)]
if "RISK_FLAG_ADJ" not in existing:
    arcpy.management.AddField(
        CATCHMENTS, "RISK_FLAG_ADJ", "SHORT",
        field_alias="Adjusted risk flag"
    )

print("Step 5: Computing RISK_FLAG_ADJ...")
adj_fields = ["gridcode", "RISK_FLAG", "OPEN_CH_FT", "BMP_PRESENT", "RISK_FLAG_ADJ"]

adj_stats = {1: 0, 2: 0, 3: 0, 4: 0, 0: 0, -1: 0, -2: 0}

with arcpy.da.UpdateCursor(CATCHMENTS, adj_fields) as cur:
    for row in cur:
        risk    = row[1]
        open_ch = row[2] if row[2] else 0.0
        bmp     = row[3] if row[3] else 0
        has_oc  = open_ch > 0
        has_bmp = bmp == 1

        if risk == 1:
            if has_oc and has_bmp:
                adj = 4
            elif has_bmp:
                adj = 3
            elif has_oc:
                adj = 2
            else:
                adj = 1
        elif risk == 0:
            adj = 0
        else:
            # RISK_FLAG = -1 (no pipe data)
            adj = -2 if (has_oc or has_bmp) else -1

        row[4] = adj
        adj_stats[adj] += 1
        cur.updateRow(row)

# ════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ════════════════════════════════════════════════════════════════════
print("\n── FINAL SUMMARY ────────────────────────────────────────────────")
print(f"  RISK_FLAG_ADJ = 1  (at-risk, no supplemental)          : {adj_stats[1]}")
print(f"  RISK_FLAG_ADJ = 2  (at-risk, open channel present)     : {adj_stats[2]}")
print(f"  RISK_FLAG_ADJ = 3  (at-risk, flood control BMP present): {adj_stats[3]}")
print(f"  RISK_FLAG_ADJ = 4  (at-risk, both present)             : {adj_stats[4]}")
print(f"  RISK_FLAG_ADJ = 0  (adequate)                          : {adj_stats[0]}")
print(f"  RISK_FLAG_ADJ = -1 (no pipe data, no infrastructure)   : {adj_stats[-1]}")
print(f"  RISK_FLAG_ADJ = -2 (no pipe data, infrastructure present): {adj_stats[-2]}")
print(f"\nOutput: {CATCHMENTS}")
print("\nPhase 4 complete.")
