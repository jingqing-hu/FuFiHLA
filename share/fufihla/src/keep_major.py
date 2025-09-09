import sys
###########################################################
# Tool Kit for Consensus Sequence Reconstruction
###########################################################

for line in sys.stdin:
    if not line: break
    if line[0] == '#' : print(line[:-1])
    else :
        llst = line.split()
        info = llst[9].split(':')
        dp = [int(x) for x in info[2].split(',')]
        if sum(dp) < 3 :
            continue
        dp_max = max(dp)
        idx = dp.index(dp_max)
        if idx == 0 : continue
        info[0] = str(idx) + '/' + str(idx)
        llst[9] = ':'.join(info)
        print('\t'.join(llst))
