import numpy as np
import matplotlib.pyplot as plt


inputFile = 'test.txt'
outputFile = 'fourpanel.pdf'
xAxisData = 'Bin'
#xAxisData = 'T[ms]'


# read in the data
data = np.genfromtxt(inputFile, delimiter=' ', skip_header=1)

# read in the heading, and make maps of Col->Heading and Heading->Col
ifd = open(inputFile, 'r')
colToHeading = ifd.readline()[:-1].split(' ')
headingToCol = {}
for i in range(0, len(colToHeading)):
  headingToCol[colToHeading[i]] = i
ifd.close()

# now do the plots
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 7))


axes[0][0].grid(False)
for field in [('StoM', 'red'), ('MtoM','orange'), ('MtoI', 'blue')]:
  if field[0] in headingToCol:
    axes[0][0].plot(data[:,headingToCol[xAxisData]],
                    data[:,headingToCol[field[0]]],
                    color=field[1], linestyle='solid', label=field[0])
axes[0][0].set_xlabel('$t$ [ms]')
axes[0][0].set_ylabel('Utilization')
axes[0][0].legend()


axes[0][1].grid(False)
for field in [('StoT', 'red'), ('MtoT','orange'), ('LtoT', 'blue')]:
  if field[0] in headingToCol:
    axes[0][1].plot(data[:,headingToCol[xAxisData]],
                    data[:,headingToCol[field[0]]],
                    color=field[1], linestyle='solid', label=field[0])
axes[0][1].set_xlabel('$t$ [ms]')
axes[0][1].set_ylabel('Utilization')
axes[0][1].legend()


axes[1][0].grid(False)
for field in [('MtoI', 'red'), ('ItoI','orange'), ('LtoL', 'blue')]:
  if field[0] in headingToCol:
    axes[1][0].plot(data[:,headingToCol[xAxisData]],
                    data[:,headingToCol[field[0]]],
                    color=field[1], linestyle='solid', label=field[0])
axes[1][0].set_xlabel('$t$ [ms]')
axes[1][0].set_ylabel('Utilization')
axes[1][0].legend()


axes[1][1].grid(False)
for field in [('ELCO', 'red'), ('Total','blue')]:
  if field[0] in headingToCol:
    axes[1][1].plot(data[:,headingToCol[xAxisData]],
                    data[:,headingToCol[field[0]]],
                    color=field[1], linestyle='solid', label=field[0])
axes[1][1].set_xlabel('$t$ [ms]')
axes[1][1].set_ylabel('Utilization')
axes[1][1].legend()


plt.tight_layout()
plt.savefig(outputFile)
