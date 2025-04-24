import pandas as pd
import geopandas as gpd
from geopy.geocoders import Nominatim  # necessary installation, unnecessary import
import csv
from shapely.geometry import Point
##import geojson

"""
def camel_case_split(str):
    # adapted from https://www.geeksforgeeks.org/python-split-camelcase-string-to-individual-strings/
    words = [[str[0]]]

    for c in str[1:]:
        if words[-1][-1].islower() and c.isupper():
            words.append(list(c))
        else:
            words[-1].append(c)

    words = ["".join(word) for word in words]
    words[0] = words[0].capitalize()
    return " ".join(words)
"""

"""
def geocode_hubs(file="hubs.csv"):
    Generates a GeoJSON file from geocoded locations of hub (detailed hubs.csv)

    hubs = pd.read_csv(file)["hub"].tolist()
    hubs_cc = [camel_case_split(hub) + ", Texas" for hub in hubs]

    geohubs = gpd.tools.geocode(
        hubs_cc, provider="nominatim", user_agent="HOwDI"
    ).set_index(pd.Series(hubs, name="hub"))

    geohubs["County"] = [
        word.lstrip()
        for words in geohubs["address"].str.split(",")
        for word in words
        if ("County") in word
    ]

    return geohubs
"""
## adding a get_county function
def get_county(latitude,longitude):
    geolocator = Nominatim(user_agent="HOwDI")
    location = geolocator.reverse((latitude,longitude), exactly_one=True, timeout=5)
    if location:
        address = location.raw
        county = address["address"]["county"]
        return county
    else:
        return None
## adding a get_state function    
def get_state(latitude,longitude):
    geolocator = Nominatim(user_agent="HOwDI")
    location = geolocator.reverse((latitude,longitude), exactly_one=True, timeout=10)
    if location:
        address = location.raw
        state = address["address"]["state"]
        return state
    else:
        return None
    
def geocode_hubs(file="hub_dir/hubs.csv"):
    hubs = pd.read_csv(file)
    geometry = [Point(xy) for xy in zip(hubs['longitude'],hubs['latitude'])]

    geohubs = gpd.GeoDataFrame(hubs, geometry=geometry)
    geohubs["County"]= [get_county(latitude,longitude)
                        for latitude,longitude in zip(hubs['latitude'],hubs['longitude'])]
    geohubs["State"]= [get_state(latitude,longitude)
                        for latitude,longitude in zip(hubs['latitude'],hubs['longitude'])]
    geohubs = geohubs.drop(columns=['latitude','longitude'])
    return geohubs

def main():
    geohubs = geocode_hubs()
    geohubs.to_file("hubs.geojson", driver="GeoJSON")
    geohubs.to_csv("geohubs.csv", index=False)


if __name__ == "__main__":
    main()
