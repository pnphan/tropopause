import pickle
import matplotlib.pyplot as plt

with open("monthly_means/USM00091285.pkl", "rb") as f:
    data = pickle.load(f)

print(data)  # Check the contents
print(type(data))

X = []
y = []
counter = 0
years = list(data.keys())
for year in years:
    months = list(data[year].keys())
    for month in months:
        avg = data[year][month]
        counter+=1
        X.append(counter)
        y.append(avg)

plt.figure()
plt.plot(X, y)
plt.show()