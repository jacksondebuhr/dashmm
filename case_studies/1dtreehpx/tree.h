#ifndef __TREE_H__
#define __TREE_H__


#include "particle.h"


//TODO turn this into a class, and then hide the details of the class 
// elsewhere. So probably just a forward declaration here. Then the moment()
// and compute_moments will return Moment * instead of the data directly.
struct moment_t {
  double mtot_;
  double xcom_;
  double Q00_;
};


class Node {
 public:
  //Constructing a node with the address directly will not transfer ownership
  Node(hpx_addr_t data = HPX_NULL);
  
  //This constructor creates an object with ownership
  Node(double lowbound, double highbound, hpx_addr_t colocate = HPX_NULL);
  
  //This does not transfer ownership
  Node(const Node &other); 
  
  //This does transfer ownership (if any)
  Node(Node &&other);
  ~Node();
  
  //Transfer ownership (if any)
  Node &operator=(Node &&other);
  
  //Copy data only, no ownership is transferred
  Node &operator=(const Node &other);
  
  //Do we want methods for construction of the backing data after the fact?
  // sync and async versions? These would likely reuse the actions of the 
  // constructor.
  
  
  hpx_addr_t data() const {return data_;}
  void drop_ownership() {owner_ = false;}
  
  void partition_sync(hpx_addr_t parts, int first, int last, int limit);
  hpx_addr_t partition(hpx_addr_t parts, int first, int last, int limit, 
                       hpx_addr_t sync = HPX_NULL);
  
  moment_t compute_moment_sync();
  hpx_addr_t compute_moment(hpx_addr_t sync = HPX_NULL);
  
  //These first two only continue a double
  double compute_potential_sync(double pos, int part, double theta_c);
  hpx_addr_t compute_potential(double pos, int part, double theta_c, 
                               hpx_addr_t sync = HPX_NULL);
  //This one will continue a double and an LCO address. It is an error for
  // sync to be HPX_NULL
  void compute_potential_and_continue(double pos, int part, double theta_c,
                                      SetApproxContinuation cont,
                                      hpx_addr_t sync);
 
 private:
  hpx_addr_t data_;
  bool owner_;          //indicates if deletion of this object implies freeing
                        // the global allocation backing this object.
};


#endif
