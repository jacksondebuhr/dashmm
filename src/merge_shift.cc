// =============================================================================
//  This file is part of:
//  Dynamic Adaptive System for Hierarchical Multipole Methods (DASHMM)
//
//  Copyright (c) 2015-2016, Trustees of Indiana University,
//  All rights reserved.
//
//  DASHMM is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DASHMM is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DASHMM. If not, see <http://www.gnu.org/licenses/>.
//
//  This software was created at the Indiana University Center for Research in
//  Extreme Scale Technologies (CREST).
// =============================================================================


/// \file src/merge_shift.cc
/// \brief Declare merge-and-shift table

#include "builtins/merge_shift.h"


namespace dashmm {


const int merge_and_shift_table[6][6][6][3] = {
  { //[0][][][]
    { // [0][0][][]
      {dall, -1, -1}, {sall, -1, -1}, {sall, -1, -1},
      {sall, -1, -1}, {sall, -1, -1}, {uall, -1, -1}
    },
    { // [0][1][][]
      {dall, -1, -1}, {wall, -1, -1}, {wall, -1, -1},
      {wall, -1, -1}, {wall, -1, -1}, {uall, -1, -1}
    },
    { // [0][2][][]
      {dall, -1, -1}, {wall, -1, -1}, {wall, -1, -1},
      {wall, -1, -1}, {wall, -1, -1}, {uall, -1, -1}
    },
    { // [0][3][][]
      {dall, -1, -1}, {wall, -1, -1}, {wall, -1, -1},
      {wall, -1, -1}, {wall, -1, -1}, {uall, -1, -1}
    },
    { // [0][4][][]
      {dall, -1, -1}, {wall, -1, -1}, {wall, -1, -1},
      {wall, -1, -1}, {wall, -1, -1}, {uall, -1, -1}
    },
    { // [0][5][][]
      {dall, -1, -1}, {nall, -1, -1}, {nall, -1, -1},
      {nall, -1, -1}, {nall, -1, -1}, {uall, -1, -1}
    }
  },
  // [1][][][]
  {
    { // [1][0][][]
      {dall, -1, -1}, {sall, -1, -1}, {sall, -1, -1},
      {sall, -1, -1}, {sall, -1, -1}, {uall, -1, -1}
    },
    { // [1][1][][]
      {dall, -1, -1}, {d5678, s34, w2}, {s3478, w2, w6},
      {s3478, w2, w6}, {u1234, s78, w6}, {uall, -1, -1}
    },
    { // [1][2][][]
      {dall, -1, -1}, {d5678, w24, -1}, {w2468, -1, -1},
      {w2468, -1, -1}, {u1234, w68, -1}, {uall, -1, -1}
    },
    { // [1][3][][]
      {dall, -1, -1}, {d5678, w24, -1}, {w2468, -1, -1},
      {w2468, -1, -1}, {u1234, w68, -1}, {uall, -1, -1}
    },
    { // [1][4][][]
      {dall, -1, -1}, {d5678, n12, w4}, {n1256, w4, w8},
      {n1256, w4, w8}, {u1234, n56, w8}, {uall, -1, -1}
    },
    { // [1][5][][]
      {dall, -1, -1}, {nall, -1, -1}, {nall, -1, -1},
      {nall, -1, -1}, {nall, -1, -1}, {uall, -1, -1}
    }
  },
  // [2][][][]
  {
    { // [2][0][][]
      {dall, -1, -1}, {sall, -1, -1}, {sall, -1, -1},
      {sall, -1, -1}, {sall, -1, -1}, {uall, -1, -1}
    },
    { // [2][1][][]
      {dall, -1, -1}, {d5678, s34, -1}, {s3478, -1, -1},
      {s3478, -1, -1}, {u1234, s78, -1}, {uall, -1, -1}
    },
    { // [2][2][][]
      {dall, -1, -1}, {d5678, -1, -1}, {-1, -1, -1},
      {-1, -1, -1}, {u1234, -1, -1}, {uall, -1, -1}
    },
    { // [2][3][][]
      {dall, -1, -1}, {d5678, -1, -1}, {-1, -1, -1},
      {-1, -1, -1}, {u1234, -1, -1}, {uall, -1, -1}
    },
    { // [2][4][][]
      {dall, -1, -1}, {d5678, n12, -1}, {n1256, -1, -1},
      {n1256, -1, -1}, {u1234, n56, -1}, {uall, -1, -1}
    },
    { // [2][5][][]
      {dall, -1, -1}, {nall, -1, -1}, {nall, -1, -1},
      {nall, -1, -1}, {nall, -1, -1}, {uall, -1, -1}
    }
  },
  // [3][][][]
  {
    { // [3][0][][]
      {dall, -1, -1}, {sall, -1, -1}, {sall, -1, -1},
      {sall, -1, -1}, {sall, -1, -1}, {uall, -1, -1}
    },
    { // [3][1][][]
      {dall, -1, -1}, {d5678, s34, -1}, {s3478, -1, -1},
      {s3478, -1, -1}, {u1234, s78, -1}, {uall, -1, -1}
    },
    { // [3][2][][]
      {dall, -1, -1}, {d5678, -1, -1}, {-1, -1, -1},
      {-1, -1, -1}, {u1234, -1, -1}, {uall, -1, -1}
    },
    { // [3][3][][]
      {dall, -1, -1}, {d5678, -1, -1}, {-1, -1, -1},
      {-1, -1, -1}, {u1234, -1, -1}, {uall, -1, -1}
    },
    { // [3][4][][]
      {dall, -1, -1}, {d5678, n12, -1}, {n1256, -1, -1},
      {n1256, -1, -1}, {u1234, n56, -1}, {uall, -1, -1}
    },
    { // [3][5][][]
      {dall, -1, -1}, {nall, -1, -1}, {nall, -1, -1},
      {nall, -1, -1}, {nall, -1, -1}, {uall, -1, -1}
    }
  },
  // [4][][][]
  {
    { // [4][0][][]
      {dall, -1, -1}, {sall, -1, -1}, {sall, -1, -1},
      {sall, -1, -1}, {sall, -1, -1}, {uall, -1, -1}
    },
    { // [4][1][][]
      {dall, -1, -1}, {d5678, s34, e1}, {s3478, e1, e5},
      {s3478, e1, e5}, {u1234, s78, e5}, {uall, -1, -1}
    },
    { // [4][2][][]
      {dall, -1, -1}, {d5678, e13, -1}, {e1357, -1, -1},
      {e1357, -1, -1}, {u1234, e57, -1}, {uall, -1, -1}
    },
    { // [4][3][][]
      {dall, -1, -1}, {d5678, e13, -1}, {e1357, -1, -1},
      {e1357, -1, -1}, {u1234, e57, -1}, {uall, -1, -1}
    },
    { // [4][4][][]
      {dall, -1, -1}, {d5678, n12, e3}, {n1256, e3, e7},
      {n1256, e3, e7}, {u1234, n56, e7}, {uall, -1, -1}
    },
    { // [4][5][][]
      {dall, -1, -1}, {nall, -1, -1}, {nall, -1, -1},
      {nall, -1, -1}, {nall, -1, -1}, {uall, -1, -1}
    }
  },
  // [5][][][]
  {
    { // [5][0][][]
      {dall, -1, -1}, {sall, -1, -1}, {sall, -1, -1},
      {sall, -1, -1}, {sall, -1, -1}, {uall, -1, -1}
    },
    { // [5][1][][]
      {dall, -1, -1}, {eall, -1, -1}, {eall, -1, -1},
      {eall, -1, -1}, {eall, -1, -1}, {uall, -1, -1}
    },
    { // [5][2][][]
      {dall, -1, -1}, {eall, -1, -1}, {eall, -1, -1},
      {eall, -1, -1}, {eall, -1, -1}, {uall, -1, -1}
    },
    { // [5][3][][]
      {dall, -1, -1}, {eall, -1, -1}, {eall, -1, -1},
      {eall, -1, -1}, {eall, -1, -1}, {uall, -1, -1}
    },
    { // [5][4][][]
      {dall, -1, -1}, {eall, -1, -1}, {eall, -1, -1},
      {eall, -1, -1}, {eall, -1, -1}, {uall, -1, -1}
    },
    { // [5][5][][]
      {dall, -1, -1}, {nall, -1, -1}, {nall, -1, -1},
      {nall, -1, -1}, {nall, -1, -1}, {uall, -1, -1}
    }
  }
};


} // dashmm
