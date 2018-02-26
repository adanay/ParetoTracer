Notes: These notes may be not updated.

- DTLZ1, IDTLZ1, DTLZ3 are the hardest: with n large they get very slow to solve. 
  PT also does not filter lots of weak optima. 
  The distribution of points is not that uniform as in other test functions.
- DTLZ2, IDTLZ2, CDTLZ2, DTLZ4 are easy and fast to solve even with n large.
- DTL5 is degenerated for nobj > 2. PT behaves well for nobj = 2 but not for nobj > 2.
- DTLZ6 is degenerated for nobj > 2 but PT can solve it.
- DTLZ7 is disconnected. PT behaves acceptable for nobj = 2 but runs into 
  more troubles for nobj > 2.