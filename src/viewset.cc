#include "dashmm/viewset.h"


namespace dashmm {


void ViewSet::destroy() {
  for (size_t i = 0; i < views_.size(); ++i) {
    //delete if non-null
    if (views_[i].data != nullptr) {
      delete [] views_[i].data;
      views_[i].data == nullptr;
    }
  }
}


void ViewSet::add_view(int index) {
  add_view(index, 0, nullptr);
}


void ViewSet::add_view(int index, size_t bytes, char *data) {
  views_.push_back(View{index, bytes, data});
}


size_t ViewSet::bytes() const {
  size_t retval = 0;
  for (size_t i = 0; i < views_.size(); ++i) {
    retval += views_[i].bytes;
  }

  retval += 2 * sizeof(int);  // for count and n_digits
  retval += views_.size() * sizeof(int);
  retval += views_.size() * sizeof(size_t);

  return retval;
}


void ViewSet::serialize(WriteBuffer &buffer) {
  bool e = buffer.write(&n_digits_, sizeof(n_digits_));
  assert(e);
  e = buffer.write(count());
  assert(e);

  // then the view index and size for each
  for (int i = 0; i < count(); ++i) {
    e = buffer.write(view_index(i));
    assert(e);
    e = buffer.write(view_bytes(i));
    assert(e);
  }

  // then the data
  for (int i = 0; i < count(); ++i) {
    e = buffer.write(view_data(i), view_bytes(i));
    assert(e);
  }
}


void ViewSet::interpret(ReadBuffer &buffer) {
  // NOTE: this clears out the ViewSet
  views_.clear();

  bool e = buffer.read(&n_digits_);
  assert(e);

  int ct{0};
  e = buffer.read(&ct, sizeof(ct));
  assert(e);

  for (int i = 0; i < ct; ++i) {
    int idx{};
    size_t bts{};
    e = buffer.read(&idx, sizeof(idx));
    assert(e);
    e = buffer.read(&bts, sizeof(bts));
    assert(e);

    add_view(idx, bts, nullptr);
  }

  for (int i = 0; i < ct; ++i) {
    set_data(i, buffer.cursor());
    e = buffer.advance(view_bytes(i));
    assert(e);
  }
}


} // namespace dashmm
