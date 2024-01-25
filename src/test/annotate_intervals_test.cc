#include "annotate_intervals.h"

#include "gtest/gtest.h"

#include <string>
#include <vector>

namespace {

// --gtest_filter=AnnotateIntervalsTest.Basic
TEST(AnnotateIntervalsTest, Basic) {
  std::vector<int64_t> posvec({1,3,5,7,9,11,13});
  std::vector<int64_t> lower({1,1,4,7,9});
  std::vector<int64_t> upper({3,6,6,8,9});
  std::vector<int> label({1,2,2,3,4});
  std::string lm;

  int64_t size = annotate_intervals(posvec.size(), &posvec[0], label.size(), &lower[0], &upper[0], &label[0], &lm);
  ASSERT_TRUE(size == lm.size());

  char* lm_ptr = &lm[0];
  int* lm_int = reinterpret_cast<int*>(lm_ptr);

  std::vector<int> expected({0, 0, 1, 2, 3, 1, 2, 2, 2, 3});

  for (int i = 0; i < size/sizeof(int); i++) {
    ASSERT_TRUE(expected[i] == lm_int[i]);
    //std::cout << lm_int[i] << " ";
    //if ((i+1)==(size/sizeof(int)/2)) std::cout << "\n";
  }
}

}

