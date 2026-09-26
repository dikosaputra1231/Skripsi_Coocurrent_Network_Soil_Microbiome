"""
Enrich EMP Cropland Metadata with NASA POWER Climate and ISRIC SoilGrids 2.0 Physicochemical Properties.
Author: Diko Duwi Saputra
"""

import os
import io
import json
import time
import requests
import pandas as pd
import numpy as np
from PIL import Image
from concurrent.futures import ThreadPoolExecutor, as_completed

MAPPING_INPUT = 'mapping_files/emp_qiime_mapping_qc_filtered.csv'
METADATA_OUTPUT = 'mapping_files/cropland_enriched_metadata.csv'
CACHE_SOILGRIDS = 'mapping_files/soilgrids_cache.json'

SOIL_PROPERTIES = ['phh2o', 'soc', 'nitrogen', 'cec', 'clay', 'sand', 'silt', 'bdod']
DEPTH_INTERVALS = ['0-5cm', '5-15cm']


def fetch_nasa_power_moisture(lat, lon, date_str):
    """
    Fetch surface soil wetness (GWETTOP) and temperature (T2M) from NASA POWER API.
    GWETTOP: 0.0 (bone dry) to 1.0 (fully saturated/flooded)
    """
    if pd.isna(lat) or pd.isna(lon):
        return None, None
    
    if pd.notna(date_str) and str(date_str).strip() != '' and str(date_str).strip() != 'nan':
        clean_date = str(date_str).strip().replace('-', '')[:8]
        start_date = clean_date
        end_date = clean_date
    else:
        start_date = '20110601'
        end_date = '20110605'

    url = (
        f"https://power.larc.nasa.gov/api/temporal/daily/point?"
        f"parameters=GWETTOP,T2M&community=AG&longitude={lon}&latitude={lat}"
        f"&start={start_date}&end={end_date}&format=JSON"
    )
    
    try:
        resp = requests.get(url, timeout=12)
        if resp.status_code == 200:
            res_json = resp.json()
            param_dict = res_json.get('properties', {}).get('parameter', {})
            gwettop_vals = [v for v in param_dict.get('GWETTOP', {}).values() if v is not None and v != -999.0]
            t2m_vals = [v for v in param_dict.get('T2M', {}).values() if v is not None and v != -999.0]
            
            gwettop = float(np.mean(gwettop_vals)) if gwettop_vals else None
            t2m = float(np.mean(t2m_vals)) if t2m_vals else None
            return gwettop, t2m
    except Exception as e:
        print(f"Error fetching NASA POWER lat={lat}, lon={lon}, date={date_str}: {e}")
    
    return None, None


def fetch_single_wcs_layer(session, lat, lon, prop, depth, max_retries=3):
    """
    Query ISRIC SoilGrids WCS 1.0.0 for a single property and depth layer.
    Includes adaptive radius expansion for coordinates near water/nodata masks (e.g. coastal paddies).
    """
    url = f"https://maps.isric.org/mapserv?map=/map/{prop}.map"
    for delta in [0.001, 0.008, 0.025]:
        params = {
            'SERVICE': 'WCS',
            'VERSION': '1.0.0',
            'REQUEST': 'GetCoverage',
            'COVERAGE': f"{prop}_{depth}_mean",
            'CRS': 'EPSG:4326',
            'BBOX': f"{lon-delta},{lat-delta},{lon+delta},{lat+delta}",
            'WIDTH': '6',
            'HEIGHT': '6',
            'FORMAT': 'GEOTIFF_INT16'
        }
        for attempt in range(max_retries):
            try:
                r = session.get(url, params=params, timeout=20)
                if r.status_code == 200:
                    img = Image.open(io.BytesIO(r.content))
                    arr = np.array(img)
                    valid_mask = (arr > 0) & (arr < 30000)
                    if np.any(valid_mask):
                        return float(np.mean(arr[valid_mask]))
                break
            except Exception:
                time.sleep(0.5 * (attempt + 1))
    return None


def fetch_all_soilgrids_for_coord(session, lat, lon):
    """
    Extract all 8 physicochemical properties for 0-5cm and 5-15cm,
    then compute weighted topsoil average (0-15 cm) with proper scientific unit scaling.
    """
    raw_depths = {p: {} for p in SOIL_PROPERTIES}
    for prop in SOIL_PROPERTIES:
        for depth in DEPTH_INTERVALS:
            val = fetch_single_wcs_layer(session, lat, lon, prop, depth)
            raw_depths[prop][depth] = val

    # Calculate weighted average 0-15 cm: (0-5cm * 5 + 5-15cm * 10) / 15
    result = {}
    for prop in SOIL_PROPERTIES:
        v0_5 = raw_depths[prop].get('0-5cm')
        v5_15 = raw_depths[prop].get('5-15cm')
        if v0_5 is not None and v5_15 is not None:
            w_mean = (v0_5 * 5.0 + v5_15 * 10.0) / 15.0
        elif v0_5 is not None:
            w_mean = v0_5
        elif v5_15 is not None:
            w_mean = v5_15
        else:
            w_mean = None

        result[prop] = w_mean

    # Convert to standard scientific units
    # phh2o: factor 10 -> pH units
    ph = round(result['phh2o'] / 10.0, 2) if result.get('phh2o') is not None else None
    # soc: factor 10 (dg/kg -> g/kg)
    soc = round(result['soc'] / 10.0, 2) if result.get('soc') is not None else None
    # nitrogen: factor 100 (cg/kg -> g/kg)
    nitrogen = round(result['nitrogen'] / 100.0, 3) if result.get('nitrogen') is not None else None
    # cec: factor 10 (mmol(c)/kg -> cmol(+)/kg meq/100g)
    cec = round(result['cec'] / 10.0, 2) if result.get('cec') is not None else None
    # clay, sand, silt: factor 10 (g/kg -> %)
    clay = round(result['clay'] / 10.0, 2) if result.get('clay') is not None else None
    sand = round(result['sand'] / 10.0, 2) if result.get('sand') is not None else None
    silt = round(result['silt'] / 10.0, 2) if result.get('silt') is not None else None
    # bdod: factor 100 (cg/cm3 -> g/cm3)
    bdod = round(result['bdod'] / 100.0, 3) if result.get('bdod') is not None else None

    return {
        'soilgrids_ph': ph,
        'soilgrids_soc_g_kg': soc,
        'soilgrids_nitrogen_g_kg': nitrogen,
        'soilgrids_cec_cmolc_kg': cec,
        'soilgrids_clay_percent': clay,
        'soilgrids_sand_percent': sand,
        'soilgrids_silt_percent': silt,
        'soilgrids_bdod_g_cm3': bdod
    }


def main():
    print(f"Loading metadata...")
    if os.path.exists(METADATA_OUTPUT):
        print(f"Found existing enriched metadata at {METADATA_OUTPUT}. Loading it...")
        df = pd.read_csv(METADATA_OUTPUT, low_memory=False)
    else:
        print(f"Loading base mapping file: {MAPPING_INPUT}...")
        base_df = pd.read_csv(MAPPING_INPUT, low_memory=False)
        df = base_df[(base_df['empo_3'] == 'Soil (non-saline)') & (base_df['env_biome'] == 'cropland biome')].copy()
        print(f"Found {len(df)} agricultural/cropland samples.")
        
        # Labeling Ecosystem Type (Wetland vs Dryland)
        df['ecosystem_type'] = df['env_feature'].apply(
            lambda x: 'Wetland_Paddy' if 'paddy' in str(x).lower() else 'Dryland_Cropland'
        )

    # Ensure NASA POWER columns exist
    if 'soil_moisture_nasa_gwettop' not in df.columns or df['soil_moisture_nasa_gwettop'].isna().all():
        print("\nExtracting unique spatial-temporal coordinates for NASA POWER queries...")
        unique_combos = df[['latitude_deg', 'longitude_deg', 'collection_timestamp']].drop_duplicates()
        print(f"Total unique NASA queries: {len(unique_combos)}")
        
        nasa_cache = {}
        for idx, row in unique_combos.iterrows():
            lat = row['latitude_deg']
            lon = row['longitude_deg']
            dt = row['collection_timestamp']
            key = (lat, lon, str(dt))
            gw, t2m = fetch_nasa_power_moisture(lat, lon, dt)
            nasa_cache[key] = (gw, t2m)
            time.sleep(0.3)
        
        df['soil_moisture_nasa_gwettop'] = [nasa_cache.get((r['latitude_deg'], r['longitude_deg'], str(r['collection_timestamp'])), (None, None))[0] for _, r in df.iterrows()]
        df['soil_temp_nasa_t2m'] = [nasa_cache.get((r['latitude_deg'], r['longitude_deg'], str(r['collection_timestamp'])), (None, None))[1] for _, r in df.iterrows()]

    # Load or initialize SoilGrids cache
    sg_cache = {}
    if os.path.exists(CACHE_SOILGRIDS):
        try:
            with open(CACHE_SOILGRIDS, 'r') as f:
                sg_cache = json.load(f)
            print(f"Loaded {len(sg_cache)} cached SoilGrids entries from {CACHE_SOILGRIDS}.")
        except Exception as e:
            print(f"Could not load cache: {e}")

    # Extract unique lat-lon coordinates
    unique_coords = df[['latitude_deg', 'longitude_deg']].drop_duplicates().dropna()
    print(f"\nTotal unique coordinates for SoilGrids query: {len(unique_coords)}")

    # Prepare HTTP session for thread pool
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(pool_connections=20, pool_maxsize=20)
    session.mount('https://', adapter)

    missing_coords = []
    for _, row in unique_coords.iterrows():
        lat = round(float(row['latitude_deg']), 5)
        lon = round(float(row['longitude_deg']), 5)
        key = f"{lat},{lon}"
        if key not in sg_cache or any(v is None for v in sg_cache[key].values()):
            missing_coords.append((lat, lon))

    print(f"Coordinates to query from ISRIC SoilGrids WCS: {len(missing_coords)}")

    if missing_coords:
        print("Starting multi-threaded SoilGrids extraction (0-15 cm topsoil)...")
        completed = 0
        total = len(missing_coords)

        def worker(coord):
            c_lat, c_lon = coord
            res = fetch_all_soilgrids_for_coord(session, c_lat, c_lon)
            return (f"{c_lat},{c_lon}", res)

        with ThreadPoolExecutor(max_workers=6) as executor:
            future_to_coord = {executor.submit(worker, c): c for c in missing_coords}
            for future in as_completed(future_to_coord):
                key, res = future.result()
                sg_cache[key] = res
                completed += 1
                print(f"[{completed}/{total}] Fetched {key} -> pH: {res['soilgrids_ph']}, SOC: {res['soilgrids_soc_g_kg']} g/kg, Clay: {res['soilgrids_clay_percent']}%")
                
                # Periodically save cache checkpoint
                if completed % 10 == 0 or completed == total:
                    with open(CACHE_SOILGRIDS, 'w') as f:
                        json.dump(sg_cache, f, indent=2)

    # Save final cache
    with open(CACHE_SOILGRIDS, 'w') as f:
        json.dump(sg_cache, f, indent=2)
    print(f"Saved complete SoilGrids cache to {CACHE_SOILGRIDS}.")

    # Map SoilGrids values back to DataFrame
    print("\nMapping SoilGrids physicochemical properties to cropland samples...")
    for prop in ['soilgrids_ph', 'soilgrids_soc_g_kg', 'soilgrids_nitrogen_g_kg', 
                 'soilgrids_cec_cmolc_kg', 'soilgrids_clay_percent', 
                 'soilgrids_sand_percent', 'soilgrids_silt_percent', 'soilgrids_bdod_g_cm3']:
        df[prop] = df.apply(
            lambda r: sg_cache.get(f"{round(float(r['latitude_deg']), 5)},{round(float(r['longitude_deg']), 5)}", {}).get(prop),
            axis=1
        )

    # Create hybrid and transparent pH imputation columns
    print("Constructing hybrid 'ph_final' and scientific provenance flagging ('ph_source', 'is_ph_imputed')...")
    # ph_final: use original measured ph if available; otherwise impute with soilgrids_ph
    df['ph_final'] = df['ph'].combine_first(df['soilgrids_ph'])
    
    # ph_source: 'measured', 'soilgrids_imputed', or 'missing'
    df['ph_source'] = np.where(
        df['ph'].notna(), 
        'measured', 
        np.where(df['soilgrids_ph'].notna(), 'soilgrids_imputed', 'missing')
    )
    df['is_ph_imputed'] = (df['ph'].isna()) & (df['ph_final'].notna())

    # Save enriched metadata
    df.to_csv(METADATA_OUTPUT, index=False)
    print(f"\n==================================================================")
    print(f"Successfully enriched and saved metadata to: {METADATA_OUTPUT}")
    print(f"Total samples: {len(df)}")
    print(f"Original measured pH present: {df['ph'].notna().sum()}")
    print(f"Missing pH filled by SoilGrids: {df['is_ph_imputed'].sum()}")
    print(f"Total valid pH_final: {df['ph_final'].notna().sum()} / {len(df)}")
    print(f"Provenance breakdown: {dict(df['ph_source'].value_counts())}")
    print(f"New SoilGrids columns added:")
    print("  - soilgrids_ph (pH units, topsoil 0-15 cm)")
    print("  - soilgrids_soc_g_kg (Soil Organic Carbon, g/kg)")
    print("  - soilgrids_nitrogen_g_kg (Total Nitrogen, g/kg)")
    print("  - soilgrids_cec_cmolc_kg (Cation Exchange Capacity, cmol(+)/kg)")
    print("  - soilgrids_clay_percent, soilgrids_sand_percent, soilgrids_silt_percent (%)")
    print("  - soilgrids_bdod_g_cm3 (Bulk Density, g/cm3)")
    print("  - ph_final, ph_source, is_ph_imputed")
    print(f"==================================================================")


if __name__ == '__main__':
    main()
