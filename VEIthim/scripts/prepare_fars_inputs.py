import pandas as pd
import numpy as np
import geopandas as gpd
import shapely.vectorized as sv
from pathlib import Path

# root dir
root = Path(r"C:\VE-ITHIM\input_development\fars")

years = ["FARS2017NationalCSV", "FARS2018NationalCSV", "FARS2019NationalCSV", "FARS2020NationalCSV", "FARS2021NationalCSV"]
crash = []
person = []
for year in years:
    # develop with one dataset
    data_dir = root.joinpath(year)
    crash_yr = pd.read_csv(data_dir.joinpath("accident.csv"), encoding='latin-1', low_memory=False)
    person_yr = pd.read_csv(data_dir.joinpath("person.csv"), encoding='latin-1', low_memory=False)

    # add year designation to data
    year = year.split("FARS")[-1].split("National")[0]
    crash_yr["year"] = year
    person_yr["year"] = year
    crash.append(crash_yr)
    person.append(person_yr)

person = pd.concat(person, ignore_index=True)
crash = pd.concat(crash, ignore_index=True)

# isolate Oregon cases
person = person[person.STATENAME.isin(["Oregon"])].copy()

# person type to occupnat categories
mode_map = {
    "Driver of a Motor Vehicle In-Transport": "auto",
    "Passenger of a Motor Vehicle In-Transport": "auto",
    "Unknown Occupant Type in a Motor Vehicle In- Transport": "auto",
    "Pedestrian": "walk",
    "Bicyclist": "bike",
    "Other Cyclist": "bike",
    "Person on Non-Motorized Personal Conveyance": "bike",  # think these are micro-mobility modes
    "Person on Personal Conveyance, Unknown if Motorized or Non-Motorized": "bike",
    "Person on Motorized Personal Conveyanced": "bike",
    "Occupant of a Motor Vehicle Not In- Transport": "other",
    "Occupant of a Non-Motor Vehicle Transport Device": "other",
    "Persons In/On Buildings": "other"
}
person["veh_mode_type"] = person.PER_TYPNAME.map(mode_map).fillna("other")

# identify fatals, other parties involved
person["uid"] = person['year'].astype(str) + "|" + person['ST_CASE'].astype(str) \
    + "|" + person['VEH_NO'].astype(str) + "|" + person['PER_NO'].astype(str)
fatals = person[person.INJ_SEV==4]
other_parties = person[person.INJ_SEV!=4]

# this is going to be slow, but let's just loop through fatal cases
all_fatal_cases = []
fatal_cases = fatals.uid.unique()
for i, uid in enumerate(fatal_cases):

    if (i % 100)==0:
        print(i, "cases complete")

    # subset
    case_fatals = fatals[fatals.uid==uid]
    others_involved = ((person.ST_CASE==case_fatals.ST_CASE.iloc[0]) & (person.VEH_NO != case_fatals.VEH_NO.iloc[0]))
    if others_involved.any():
        others_involved = person[others_involved].copy()

        # finally, get mode(s) for other involved parties
        other_modes = others_involved.groupby("ST_CASE").agg(
            strike_mode=("veh_mode_type", "unique")
        )

        strike_mode = other_modes.strike_mode.values[0]
        if "other" in strike_mode:
            strike_mode = np.setdiff1d(strike_mode, ["other"])
            if len(strike_mode) < 1:
                strike_mode = [np.nan]
        if len(strike_mode) > 1:
            if "auto" in strike_mode:
                strike_mode = ["auto"]
            if "bike" in strike_mode and "walk" in strike_mode:
                strike_mode = ["bike"]

    else:
        strike_mode = [np.nan]
    
    # keep required cols for fatal case
    keep = ["year", "ST_CASE", "AGE", "SEXNAME", "veh_mode_type"]
    fatal_case = case_fatals[keep].copy()
    fatal_case["strike_mode"] = [strike_mode]

    # append
    all_fatal_cases.append(fatal_case)

# concat cases
all_fatal_cases = pd.concat(all_fatal_cases)
all_fatal_cases["strike_mode"] = all_fatal_cases.strike_mode.apply(
    lambda x: x if isinstance(x, list) else x.tolist())
all_fatal_cases["strike_mode"] = all_fatal_cases.strike_mode.apply(
    lambda x: x[0] if isinstance(x, list) else np.nan)

# join geo
join_cols = ["year", "ST_CASE", "LATITUDE", "LONGITUD"]
all_fatal_cases = all_fatal_cases.merge(crash[join_cols], on=["year", "ST_CASE"], how="left")

# clean up schema
all_fatal_cases.columns = ['year', 'case_no', 'age', 'sex', 'cas_mode', 'strike_mode', 'lat', 'lng']

# finally, subset to crahses inside region
region = gpd.read_parquet(r"C:\VE-ITHIM\input_development\boundary.parquet")
in_region = sv.contains(region.geometry[0], all_fatal_cases.lng, all_fatal_cases.lat)

# check
fars_gdf = all_fatal_cases.copy()
fars_gdf["inside"] = in_region
fars_gdf["geometry"] = gpd.points_from_xy(fars_gdf["lng"], fars_gdf["lat"], crs=region.crs)
fars_gdf = gpd.GeoDataFrame(fars_gdf, geometry="geometry")
fars_gdf.to_parquet(r"C:\VE-ITHIM\input_development\fars.parquet")

# isolate cases in region
all_fatal_cases = all_fatal_cases[in_region].copy()

# recode instances where cas_mode is auto, strike_mode is walk/bike
recode_bike = ((all_fatal_cases.cas_mode=="auto") & (all_fatal_cases.strike_mode=="bike"))
all_fatal_cases["cas_mode"] = np.where(recode_bike, "bike", all_fatal_cases.cas_mode)
all_fatal_cases["strike_mode"] = np.where(recode_bike, "auto", all_fatal_cases.strike_mode)
recode_walk = ((all_fatal_cases.cas_mode=="auto") & (all_fatal_cases.strike_mode=="walk"))
all_fatal_cases["cas_mode"] = np.where(recode_walk, "walk", all_fatal_cases.cas_mode)
all_fatal_cases["strike_mode"] = np.where(recode_walk, "auto", all_fatal_cases.strike_mode)
all_fatal_cases.to_csv(r"C:\VE-ITHIM\input_development\fars.csv", index=False)

