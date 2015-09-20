#ifndef __INDEX_H__
#define __INDEX_H__


namespace dashmm {


//
class Index {
 public:
  Index(int ix, int iy, int iz, int level)
      : idx_{ix, iy, iz}, level_{level} { }

  int x() const {return idx_[0];}
  int y() const {return idx_[1];}
  int z() const {return idx_[2];}
  int level() const {return level_;}

  Index parent(int num = 1) const {
    return Index{idx_[0] >> num, idx_[1] >> num, idx_[2] >> num,
                 level_ > num ? level_ - num : 0};
  }
  Index child(int which) const {
    return Index{(idx_[0] << 1) + (which & 1 ? 1 : 0),
                 (idx_[1] << 1) + (which & 2 ? 1 : 0),
                 (idx_[2] << 1) + (which & 4 ? 1 : 0),
                 level_ + 1};
  }

 private:
  int idx_[3];
  int level_;
};


} // namespace dashmm


#endif
