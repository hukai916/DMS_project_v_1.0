import pandas as pd

def md2list(md, ref_pos_start):
    numList = []
    posRead = []
    ind = 0
    base = 0
    for index, val in enumerate(md):
        if val.isupper():
            numList.append(md[ind:index])
            ind = index + 1
    for i in range(0, len(numList)):
        posRead.append(base + int(numList[i]) + 1 + ref_pos_start)
        base = base + int(numList[i]) + 1
    return(posRead) # return a list of index of mutated site on the reference.
