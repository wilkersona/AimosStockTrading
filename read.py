# scp PCPCnttr@blp01.ccni.rpi.edu:/gpfs/u/scratch/PCPC/PCPCnttr/project/'{strats,buyers,stocks}' ./

import struct

TOTAL_BUYERS = 5000
TOTAL_STOCKS = 500
NUM_T = 100

strat_strings = {
    0: "Rise",
    1: "Fall",
    2: "Peak",
    3: "Dip",
}

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
        row = []
        for j in range(TOTAL_BUYERS):
            row.append(out[i*TOTAL_BUYERS+j])
        buyer_values.append(row)

with open("stocks", "rb") as f:
    rowmax = [-1 for _ in range(TOTAL_STOCKS)]
    rowmin = [-1 for _ in range(TOTAL_STOCKS)]
    out = struct.unpack('f'*TOTAL_STOCKS*(NUM_T+1), f.read(4*TOTAL_STOCKS*(NUM_T+1)))
    for i in range(NUM_T+1):
        row = []
        for j in range(TOTAL_STOCKS):
            if rowmax[j] < out[i*TOTAL_STOCKS+j]:
                rowmax[j] = out[i*TOTAL_STOCKS+j]
            if rowmin[j] == -1 or rowmin[j] > out[i*TOTAL_STOCKS+j]:
                rowmin[j] = out[i*TOTAL_STOCKS+j]
            row.append(out[i*TOTAL_STOCKS+j])
        stock_values.append(row)
    rowdiff = [rowmax[i]-rowmin[i] for i in range(TOTAL_STOCKS)]

if False:
    stock_order = sorted(enumerate(rowdiff), key=lambda x: x[1], reverse=True)
    top_10_stocks = [pair[0] for pair in stock_order][:10]

    with open("top10stocks.csv", "w") as outfile:
        out = "Time Step,"
        for index in top_10_stocks:
            out += f"{index},"
        out = out[:-1] + "\n"
        for t, time_step in enumerate(stock_values):
            out += f"{t},"
            for i in top_10_stocks:
                out += f"{time_step[i]},"
            out = out[:-1] + "\n"
        outfile.write(out)

if False:
    with open("top10div.csv", "w") as outfile:
        out = "Time Step,"
        for index in range(10):
            out += f"{index},"
        out = out[:-1] + "\n"
        for t, time_step in enumerate(stock_values):
            out += f"{t},"
            for i in range(10):
                out += f"{time_step[i]},"
            out = out[:-1] + "\n"
        outfile.write(out)

if True:
    strat_list = []
    strat_keys = {}
    strat_avgs = []

    for buy, sell, comm in strats:
        if not (buy, sell) in strat_keys:
            strat_keys[(buy, sell)] = {"count": 0, "index": len(strat_list)}
            strat_list.append((buy, sell))
        strat_keys[(buy, sell)]["count"] += 1

    for time_step in buyer_values:
        t_total = [0 for _ in range(len(strat_list))]
        for i, value in enumerate(time_step):
            strat = strats[i][:2]
            t_total[strat_keys[strat]["index"]] += value
        for i, total in enumerate(t_total):
            t_total[i] = total / strat_keys[strat_list[i]]["count"]
        strat_avgs.append(t_total)

    with open("strat_avgs.csv", "w") as outfile:
        out = "Time Step,"
        for buy, sell in strat_list:
            out += f"{strat_strings[buy]} {strat_strings[sell]},"
        out = out[:-1] + "\n"
        for t, time_step in enumerate(strat_avgs):
            out += f"{t},"
            for val in time_step:
                out += f"{val},"
            out = out[:-1] + "\n"
        outfile.write(out)
