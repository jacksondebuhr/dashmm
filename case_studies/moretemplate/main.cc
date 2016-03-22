#include <iostream>

// #include <hpx/hpx.h>


template <template <typename> class B>
class BigClass {
 public:
  class FirstClass {
   public:
    static void report(int a) {
      fprintf(stdout, "This is class 1 reporting %d\n", a);
    }
  };

  BigClass(int val) : member_{val} { }

  void cause_report() {
    member_.report();
  }
 private:
  B<FirstClass> member_;
};


template <typename A>
class SecondClass {
 public:
  SecondClass(int number) : data_{number} { }

  void report() {
    A::report(data_);
  }
 private:
  int data_;
};




// Now for a second approach
template <typename Source, template <typename, typename> class Method>
class dashmm {
 public:
  class SourceRef {
   public:
    SourceRef(const Source &src) {
      data_ = src;
    }

    void contribute(int val) {
      data_.value += val;
    }

    void report() {
      fprintf(stdout, "I am %d\n", data_.value);
    }

   private:
    Source data_;
  };

  dashmm() : member_{} { }

  void perform_evaluation(const Source &src, const Method<Source, SourceRef> &met) {
    member_ = met;
    SourceRef sref{src};
    sref.report();
    member_.dothework(sref);
    sref.report();
  }

 private:
  Method<Source, SourceRef> member_;
};

template <typename Source, typename SourceRef>
class somemethod {
 public:
  somemethod(int val = 0) : data_{val} { }

  void dothework(SourceRef &ref) {
    ref.contribute(data_);
  }
 private:
  int data_;
};

struct simplesource {
  int value;
};


int main(int argc, char **argv) {
  BigClass<SecondClass> test{11};
  test.cause_report();


  dashmm<simplesource, somemethod> instance{};
  simplesource bob{12};
  somemethod<simplesource,
             dashmm<simplesource, somemethod>::SourceRef> alice{14};


  instance.perform_evaluation(bob, alice);


  return 0;
}
