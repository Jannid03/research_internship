import numpy as np
import matplotlib.pyplot as plt

# deathrate = 0.456713
# birthrate = 136.0967
## Werte aus R

#Normaler Gillespie Algorithmus
def run(t_final, dr):

    deathrate = dr #0.06094
    birthrate = deathrate * 297.9918


    erg = [0]
    t_erg = [0]
    speeds = []
    n = 0
    t = 0
    while t < t_final:

        rates = [birthrate, n*deathrate]

        delta = 1/np.sum(rates)*np.log(1/np.random.uniform(0,1))

        if (t+delta) > t_final:
            # print(speeds)
            return erg, t_erg, speeds
        else:
            t += delta
        

        if (np.random.uniform(0,1) < birthrate/np.sum(rates)):
            n += 1
            speeds.append((1)/delta*7/1000)
        else:
            n -= 1
            speeds.append((-1)/delta*7/1000)
        erg.append(n)
        t_erg.append(t)

        if (n <= 0): break
    
    # print(speeds)
    return erg, t_erg, speeds

sums = []
# print(range(0,1,0.1))
sequence = list(np.arange(0.0609,0.0610,0.00001)) ##Finale range
for i in [0.06094]:#sequence:
    t_final = 100
    tries = 10
    results = []
    for j in range(tries):
        xs, ts, spee = run(t_final, i)
        results.extend(spee)
        plt.plot(ts, xs)
        # if(j%10 == 0): print(j)

    results = np.sort(results)#[2500:-1500]
# print(results/)
### Vergeleich mit Originaldaten
    n,bins,bars = plt.hist(results, bins = [-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2 ,0.0 ,0.2 ,0.4,  0.6,0.8 , 1.0 , 1.2 , 1.4], density=True)
    ori = [0.01014199, 0.00000000, 0.02028398, 0.00000000, 0.10141988, 0.09127789,0.30425963, 0.51724138, 0.63894523, 0.90263692, 0.96348884, 0.63894523,0.40567951, 0.27383367, 0.09127789, 0.02028398, 0.02028398]
    sum = 0
    for j in range(len(n)):
        sum += (ori[j]-n[j])**2
    print(sum)
    sums.append(sum)

plt.plot([-1.9,-1.7,-1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1 ,0.1 ,0.3 ,0.5,  0.7,0.9 , 1.1 , 1.3 ],[0.01014199, 0.00000000, 0.02028398, 0.00000000, 0.10141988, 0.09127789,0.30425963, 0.51724138, 0.63894523, 0.90263692, 0.96348884, 0.63894523,0.40567951, 0.27383367, 0.09127789, 0.02028398, 0.02028398])
plt.show()
plt.plot(sequence, sums)
plt.show()