# scp PCPCnttr@blp01.ccni.rpi.edu:/gpfs/u/scratch/PCPC/PCPCnttr/project/'{strats,buyers,stocks}' ./

import struct

TOTAL_BUYERS = 100
TOTAL_STOCKS = 10
NUM_T = 2

strats = []
strat_counts = {}
buyer_values = []
stock_values = []

with open("strats", "rb") as f:
    out = struct.unpack('f'*TOTAL_BUYERS*3, f.read(4*TOTAL_BUYERS*3))
    for i in range(len(out)//3):
        buy = out[i*3]
        sell = out[i*3+1]
        comm = out[i*3+2]
        strats.append((buy, sell, comm))

with open("buyers", "rb") as f:
    out = struct.unpack('f'*TOTAL_BUYERS*(NUM_T+1), f.read(4*TOTAL_BUYERS*(NUM_T+1)))
    for i in range(NUM_T+1):
        for j in range(TOTAL_BUYERS):
            buyer_values.append(out[i*TOTAL_BUYERS+j])

with open("stocks", "rb") as f:
    out = struct.unpack('f'*TOTAL_STOCKS*(NUM_T+1), f.read(4*TOTAL_STOCKS*(NUM_T+1)))
    for i in range(NUM_T+1):
        for j in range(TOTAL_STOCKS):
            stock_values.append(out[i*TOTAL_STOCKS+j])

print(len(strats), len(buyer_values), len(stock_values))

for buy, sell, comm in strats:
    if not (buy, sell) in strat_counts:
        strat_counts[(buy, sell)] = 0
    strat_counts[(buy, sell)] += 1

print(len(strat_counts), strat_counts)