from functions import *


month = 5
day = 15
stepsize=50
differences=[]
for day in range(1,30):
    gph = get_daywise_data(2024, month)[(day,0)]['gph']
    temp = get_daywise_data(2024, month)[(day,0)]['temp']
    gph, temp = interpolate_to_points(gph, temp, kind='linear', step=stepsize)
    gph, tempgrad = temp_gradient_zangl(gph, temp)
    tropopause1 = detect_tropopause_zangl(tempgrad, gph, temp)

    gph = get_daywise_data(2024, month)[(day,0)]['gph']
    temp = get_daywise_data(2024, month)[(day,0)]['temp']
    gph, temp = interpolate_to_points(gph, temp, kind='linear', step=stepsize)
    gph, lapse, temp = lapse_rate(gph, temp)
    tropopause2 = detect_tropopause(gph, lapse)

    if tropopause1 != None and tropopause2 != None:
        differences.append(tropopause1*1000 - tropopause2)

print(np.mean(np.sqrt(np.array(differences) ** 2)))