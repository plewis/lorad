import sys

assert len(sys.argv) == 3, 'expecting 2 command line arguments to filter command (input and output parameter file names)'
infname  = sys.argv[1]
outfname = sys.argv[2]

filtered = {}
lines = open(infname, 'r').readlines()
headers = lines[0].strip().split('\t')
total_n = 0
for line in lines[1:]:
    parts = line.strip().split('\t')
    i = headers.index("topology")
    t = int(parts[i])
    if t in filtered.keys():
        filtered[t].append(line)
    else:
        filtered[t] = [line]
    total_n += 1

# Find most common topology
best_t = None
best_n = None
for t in filtered.keys():
    n = len(filtered[t])
    if best_n is None or n > best_n:
        best_n = n
        best_t = t
print('Most frequently sampled topology occurred in %d of %d samples' % (best_n, total_n))
best_freq = float(best_n)/total_n
print('Frequency of most frequently sampled topology is %.5f' % best_freq)

# Record samples for best topology to output file
outf = open(outfname, 'w')
outf.write(lines[0])
for line in filtered[best_t]:
    outf.write(line)
outf.close()
