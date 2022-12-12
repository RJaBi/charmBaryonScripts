"""
Takes the spreadsheet of experimental data surveyed from the PDG
and averages masses where appropriate (i.e. across charge families)
"""

import sys
import pandas as pd
import gvar as gv
import numpy as np

def aveList(row):
    """
    Averages the data in list using gvar
    """
    if isinstance(row, float):
        return None
    rowGV = []
    err = 0
    for rr in row[1:-1].split(','):
        rowGV.append(gv.mean(gv.gvar(rr)))
        err = err + gv.sdev(gv.gvar(rr))**2.0
    return gv.gvar(np.mean(rowGV), err**0.5)
        


def main(args: list):
    """
    Reads the file
    Averages
    Saves a file
    """

    df = pd.read_csv(args[0])
    df = df.rename(columns=lambda x: x.strip())
    print(df.columns)
    df['EPAve'] = df.apply(lambda row: aveList(row['EP']), axis = 1)
    df['EMAve'] = df.apply(lambda row: aveList(row['EM']), axis = 1)
    # print(df[
    print(df['EP'])
    print(df['EPAve'])
    print(df['EM'])
    print(df['EMAve'])
    #print(df)
    saveName = args[0].split('.')[0] + '_AVERAGED' + '.' + args[0].split('.')[1]
    print(saveName)
    df.to_csv(saveName)

    

    


if __name__ == '__main__':
    main(sys.argv[1:])
