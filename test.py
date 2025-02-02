import pickle
import matplotlib.pyplot as plt

stationid='USM00091285'
with open(f"tropopauses/{stationid}.pkl", "rb") as f:
    data = pickle.load(f)

print(data['second_temp']['2020']['01'].keys())