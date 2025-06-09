"""
Goal of the simulation:
Stochastically simulate the ODEs using a Gillespie algorithm:

dx/dt = rx(1 - (x + y)/k)
dy/dt = urx(1 - (x + y)/k)
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import sys
import math
import os
import traceback


def run_simulation(seed):
    try:
        print(f"[{seed}] Simulation started.", flush=True)

        random.seed(seed)
        np.random.seed(seed)

        output_filename = f"replicate_{seed}.txt"
        full_path = os.path.abspath(output_filename)
        print(f"[{seed}] Output will be written to: {full_path}", flush=True)

        with open(output_filename, "w") as output_file:
            output_file.write("time\twildtype_total\tmutant_total\tcurrent_size\n")

            # Initialization
            x, y, t = 1, 0, 0
            X, Y, T, S = [x], [y], [t], [x + y]
            Patches = 100
            tend = 50
            k = 10000
            r1, r2, u, d = 0.5, 0.49, 5e-5, 0.05
            alpha, epsilon = 0.5, 1.5
            count = 0
            growth_law = 'saturated'
            threshold = 50000

            #while t < tend:
            while True:
                count += 1

                growth_factor = max(0, 1 - ((x + y) / k))

                if growth_law == 'saturated':
                    g = alpha * ((epsilon * (x + y) / k) - 1)
                    g = max(0, g * Patches**(1 / 3))
                    rates = [
                        r1 * x * growth_factor * (1 - u),
                        d * x,
                        u * r1 * x * growth_factor + r2 * y * growth_factor,
                        d * y,
                        g * (k**(2 / 3))
                    ]
                elif growth_law == 'logistic':
                    rates = [
                        r1 * x * growth_factor * (1 - u),
                        d * x,
                        u * r1 * x * growth_factor + r2 * y * growth_factor,
                        d * y,
                        0
                    ]
                else:  # Exponential
                    rates = [
                        r1 * x * (1 - u),
                        d * x,
                        u * r1 * x + r2 * y,
                        d * y,
                        0
                    ]

                rate_sum = sum(rates)
                if rate_sum == 0:
                    print(f"[{seed}] Population crashed. Ending early.", flush=True)
                    print("Population crashed")
                    output_file.close()
                    os.remove(output_filename)  # Remove the bad file                  
                    return None

                tau = -math.log(random.random()) / rate_sum
                t += tau
                rand = random.uniform(0, 1) * rate_sum

                if rand <= rates[0]:
                    x += 1
                elif rand <= rates[0] + rates[1]:
                    x -= 1
                elif rand <= sum(rates[:3]):
                    y += 1
                elif rand <= sum(rates[:4]):
                    y -= 1
                else:
                    k += 1  # carrying capacity increases

                current_size = x + y
                if count == 100:
                    output_file.write(f"{t:.5f}\t{x}\t{y}\t{current_size}\n")
                    X.append(x)
                    Y.append(y)
                    T.append(t)
                    S.append(current_size)
                    count = 0

                if  current_size >= threshold: # should be 8000 equilibirum around, 7200 works
                    
                    break                 

        print(f"[{seed}] Simulation complete. Output written to {output_filename}", flush=True)
        return seed

    except Exception as e:
        print(f"[ERROR] Simulation failed for seed {seed}: {e}", flush=True)
        traceback.print_exc()
        return None

'''
if __name__ == "__main__":
    print("Script is starting...", flush=True)

    target_replicates = 10  # How many successful replicates you want
    successful_replicates = 0
    base_seed = 0
    attempts = 0

    while successful_replicates < target_replicates:
        seed = base_seed + attempts * 1000  # Or any increment to avoid collisions
        print(f"Attempting replicate {successful_replicates + 1} with seed {seed}...", flush=True)

        result = run_simulation(seed)
        if result is not None:
            successful_replicates += 1
        else:
            print(f"[{seed}] Failed. Retrying with a new seed...", flush=True)

        attempts += 1

    print("All successful replicates finished.")

'''
from multiprocessing import Pool, cpu_count
from itertools import count

if __name__ == "__main__":
    print("Script is starting...", flush=True)

    target_replicates = 100  # How many successful replicates you want - 10
    successful_replicates = 0
    base_seed = 0
    attempts = 0
    used_seeds = set()
    successful_seeds = []

    seed_gen = count(start=base_seed, step=1000)

    while successful_replicates < target_replicates:
        batch_size = min(cpu_count(), target_replicates - successful_replicates)
        batch_seeds = []

        # Generate a batch of unused seeds
        while len(batch_seeds) < batch_size:
            seed = next(seed_gen)
            if seed not in used_seeds:
                used_seeds.add(seed)
                batch_seeds.append(seed)

        print(f"Trying batch: {batch_seeds}", flush=True)

        with Pool(processes=batch_size) as pool:
            results = pool.map(run_simulation, batch_seeds)

        for seed, result in zip(batch_seeds, results):
            if result is not None:
                successful_replicates += 1
                successful_seeds.append(seed)
                print(f"[{seed}] ✅ Success ({successful_replicates}/{target_replicates})", flush=True)
            else:
                print(f"[{seed}] ❌ Failed — will try another.", flush=True)

    print("✅ All successful replicates finished.")
