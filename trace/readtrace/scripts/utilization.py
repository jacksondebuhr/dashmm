#!/usr/bin/python

# so this is going to start as a fixed script, and once it works, will
# be improved into a more general use thing.


# Here are the possible parameters that will be made into arguments on the
# command line.
dbFile = 'trace.db'
nBins = 100
outFile = 'test.txt'



import sqlite3

# Create a connection to the db
conn = sqlite3.connect('trace.db')

# Now create a cursor
c = conn.cursor()

# Get the bounds of the event times
c.execute('SELECT min(TimeNs), max(TimeNs) FROM Event;')
row = c.fetchone()
assert row != None
(tMin, tMax) = row
tDelta = tMax - tMin
binDelta = tDelta / nBins
binBoundaries = [tMin + x * tDelta / nBins for x in range(0, nBins + 1)]

# the data we are collecting is a dictionary from segment type into the
# total utilization of that type.
utilization = {}

# We loop over segments, compute the bins the segment overlaps, and then
# makes contributions.
for row in conn.execute('SELECT SegmentId, StartNs, EndNs FROM Segment;'):
  # Give the segment type a set of zero bins if this is the first in that
  # segment
  if not row[0] in utilization:
    utilization[row[0]] = [0 for x in range(0, nBins)]

  # Compute which bins the segment overlaps
  sBin = (row[1] - tMin) / binDelta
  eBin = (row[2] - tMin) / binDelta
  if eBin == nBins:
    eBin = nBins - 1

  # if they are the same, all contribution goes to one bin, otherwise, in
  # each bin with overlap
  if sBin == eBin:
    utilization[row[0]][sBin] += row[2] - row[1]
  else:
    utilization[row[0]][sBin] += binBoundaries[sBin + 1] - row[1]
    utilization[row[0]][eBin] += row[2] - binBoundaries[eBin]
    for x in range(sBin + 1, eBin):
      # should I just use binDelta here?
      utilization[row[0]][x] += binBoundaries[x + 1] - binBoundaries[x]

# next get the number of workers
nWorkers = 1
for row in conn.execute('SELECT count(id) FROM Worker;'):
  nWorkers = row[0];

# Now scale the ns times into utilization fractions
scaled = {}
for k, v in utilization.iteritems():
  scaled[k] = [x / (float(binDelta) * nWorkers) for x in v]

# now we get the column names for the output
headings = {}
for row in conn.execute('SELECT * FROM Segmenttype;'):
  if row[0] != 0:
    headings[row[0]] = row[1]

# collect segment ids
segments = []
for k, v in headings.iteritems():
  if k != 0:
    segments.append(k)

# compute the total
total = [0 for x in range(0, nBins)]
for k, v in scaled.iteritems():
  if k == 7:
    continue  # We skip ELCO as that double counts some event classes
  for i in range(0, nBins):
    total[i] += v[i]

# compute the bin central time in correct units (ms)
binTimes = [(binBoundaries[i] + binBoundaries[i + 1]) / (2.0 * 1.0e6)
              for i in range(0, nBins)]

# now we put out the data
ofd = open(outFile, 'w')
ofd.write("Bin T[ms] " + " ".join([v[:-1] for (k, v) in headings.items()])
          + " Total\n")
for bidx in range(0, nBins):
  ofd.write(str(bidx))
  ofd.write(" " + str(binTimes[bidx]))
  for segtype in segments:
    ofd.write(" " + str(scaled[segtype][bidx]))
  ofd.write(" " + str(total[bidx]) + "\n")
ofd.write("\n")
ofd.close()
