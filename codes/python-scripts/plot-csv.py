#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import argparse

# try:
#     import seaborn as sns
#     sns.set(style='whitegrid',font_scale=1.2)
# except:
#     print("Run:\nsudo pip3 install seaborn")




fname='sample-data-to-plot.csv'
delim=','
SIZE=(14,10)

def plot_csv(fname, output, title=None, ylabel='', delim=',', SIZE=(8,6), ymax=None):
    print(fname)
    with open(fname)as f:
        data = f.read().split('\n')
        if len(data[-1])==0:
            data=data[:-1]
        data=[_.split(delim) for _ in data]

        header=data[0]
        xlabel=header[0]
        header=header[1:]
        data=data[1:]
        data=np.array([[float(__) for __ in _] for _ in data]).transpose()
        xs=data[0]
        data=data[1:]


        plt.figure(figsize=SIZE)
        for d,h in zip(data,header):
            plt.plot(xs,d,label=h)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel, rotation=0)
        if title is not None:
            plt.title(title)
        plt.legend()

        if ymax is not None:
            plt.ylim(0,float(ymax))

        plt.grid(True)
        # plt.axis('equal')
        plt.tight_layout()
        plt.savefig(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--title', default=None, help='plot title; default: no title')
    parser.add_argument('-d', '--delim', default=',', help='CSV delimiter, default: ,')
    parser.add_argument('--ymax', default=None, help='ymax')
    parser.add_argument('-y', '--ylabel', default=None, help='Label for Y axis; default: no label')
    parser.add_argument('-i', '--input', help='CSV, comma-separated; delimiter can be changed using --delim',required=True)
    parser.add_argument('-o', '--output', help='Output file.eps',required=True)
    args = parser.parse_args()
    plot_csv(args.input, args.output, title=args.title, ylabel=args.ylabel, delim=delim, ymax=args.ymax)

