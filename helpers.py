def lat_lon_to_indices(latitude, longitude): 
    """
    Convert latitude and longitude to array indices.
    
    Args:
        latitude (float): Latitude in degrees (-90 to 90).
        longitude (float): Longitude in degrees (0 to 360).
    
    Returns:
        tuple: (latitude_index, longitude_index)
    """
    latitude_index = int((latitude + 90) / 0.25)
    longitude_index = int(longitude / 0.25)
    return latitude_index, longitude_index

def lat_lon_to_indices2(latitude, longitude): 
    """
    Convert latitude and longitude to array indices.
    
    Args:
        latitude (float): Latitude in degrees (-90 to 90).
        longitude (float): Longitude in degrees (0 to 360).
    
    Returns:
        tuple: (latitude_index, longitude_index)
    """
    longitude += 180
    latitude_index = int((latitude + 90) / 0.25)
    longitude_index = int(longitude / 0.25)
    return latitude_index, longitude_index

def indices_to_lat_lon(latitude_index, longitude_index):
    """
    Convert array indices to latitude and longitude.
    
    Args:
        latitude_index (int): Latitude index (0 to 720).
        longitude_index (int): Longitude index (0 to 1439).
    
    Returns:
        tuple: (latitude, longitude)
    """
    latitude = (latitude_index * 0.25) - 90
    longitude = longitude_index * 0.25
    return latitude, longitude