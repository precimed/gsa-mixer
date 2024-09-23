#include "annotate_intervals.h"

#include <vector>
#include <set>
#include <cstring>
#include <limits>

#include "intervaltree/IntervalTree.h"
#include "boost/icl/interval_map.hpp"
#include "bgmg_log.h"

int64_t annotate_intervals(int num_snp, int64_t* posvec, int num_interval, int64_t *interval_start, int64_t* interval_end, int* interval_label, std::string* lm) {
  // TBD - implement OMP parallelization here?
  LOG << ">annotate_intervals(num_snp=" << num_snp << ", num_interval=" << num_interval << ")";

  SimpleTimer timer(-1);

  typedef std::set<int> Labels;
  boost::icl::interval_map<int64_t, Labels> interval_map;
  for (int i = 0; i < num_interval; i++) {
    interval_map += std::make_pair(boost::icl::interval<int64_t>::right_open(interval_start[i], interval_end[i]), Labels({interval_label[i]}));
  }

  std::vector<Interval<int64_t, int> > unique_intervals;
  for (auto interval_iter = interval_map.begin(); interval_iter != interval_map.end(); ++interval_iter) {
    int64_t lower = interval_iter->first.lower();
    int64_t upper = interval_iter->first.upper();
    for (auto label_iter = interval_iter->second.begin(); label_iter != interval_iter->second.end(); ++label_iter) {
      unique_intervals.push_back(Interval<int64_t, int>(lower, upper, *label_iter));
    }
  }

  IntervalTree<int64_t, int> interval_tree(std::move(unique_intervals));
  LOG << " interval tree generated with " << unique_intervals.size() << " unique intervals";

  SimpleTimer log_timer(10000); // log some message each 10 seconds
  
  std::vector<int> annot_snp;
  std::vector<int> annot_label;
  for (int snp_index = 0; snp_index < num_snp; snp_index++) {
    if (log_timer.fire())
      LOG << " annotate_intervals still working, snp_index=" << annotate_intervals << " of " << num_snp;

    std::vector<Interval<int64_t, int> > results = interval_tree.findOverlapping(posvec[snp_index], posvec[snp_index]);
    for (auto interval: results) {
      if (interval.stop == posvec[snp_index]) continue; // treat intervals as open on the right, i.e. [a, b)
      annot_snp.push_back(snp_index);
      annot_label.push_back(interval.value);
    }
  }

  size_t bytesize = sizeof(int) * annot_snp.size();
  lm->resize(2 * bytesize);
  char* lm_ptr = &(*lm)[0];

  memcpy(lm_ptr, &annot_snp[0], bytesize);
  memcpy(lm_ptr + bytesize, &annot_label[0], bytesize);

  LOG << "<annotate_intervals(num_snp=" << num_snp << ", num_interval=" << num_interval << "), nnz=" << annot_snp.size()<< ", elapsed time " << timer.elapsed_ms() << "ms";
  return 2*bytesize;
}