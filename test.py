import pickle

with open("tropopauses/USM00091285.pkl", "rb") as f:
    data = pickle.load(f)

print(data)  # Check the contents
print(type(data))