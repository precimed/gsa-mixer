#include "gtest/gtest.h"
#include "omp.h"

#include <set>
#include <iostream>
#include <random>
#include <algorithm>

#include "bgmg_calculator.h"

/*
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
*/

namespace {

TEST(BgmgMainTest, ShouldSucceed) {
}

/*
TEST(BgmgGzipTest, TestGzip) {
  std::ifstream file("hello.z", std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
  in.push(boost::iostreams::zlib_decompressor());
  in.push(file);
  boost::iostreams::copy(in, std::cout);
}
*/

TEST(BgmgTest, LoadData) {
  // Input data:
  // - LD structure as a set of pairwise LD r2 values
  //   i1, i2, r2 <- currently this came from plink format (e.i. lower triangular, not limited to tag SNPs)
  //   Split per chromosome
  // - LD structure in a most convenient ways
  // - Plink LD structure without differentiating tag and SNP variants
  // - Several vectors, potentially with different length (#tag, #snp)
  //   zvec, nvec, hvec, weights
  // 
}

class TestMother {
public:
  TestMother(int num_snp, int num_tag, int n) : num_snp_(num_snp), num_tag_(num_tag), rd_(), g_(12341234) {
    std::uniform_real_distribution<> mafvec_dist(0.0f, 0.5f);

    for (int i = 0; i < num_snp_; i++) tag_to_snp_complete_.push_back(i);
    for (int i = 0; i < num_snp_; i++) tag_to_snp_.push_back(i);
    std::shuffle(tag_to_snp_.begin(), tag_to_snp_.end(), g_);
    tag_to_snp_.resize(num_tag_);
    std::sort(tag_to_snp_.begin(), tag_to_snp_.end()); 

    regenerate_zvec();
    for (int i = 0; i < num_tag_; i++) n_vec_.push_back(n);
    for (int i = 0; i < num_snp_; i++) n_vec_complete_.push_back(n);

    for (int i = 0; i < num_tag_; i++) weights_.push_back(1.0f);
    for (int i = 0; i < num_snp_; i++) weights_complete_.push_back(0);
    for (int i = 0; i < num_tag_; i++) weights_complete_[tag_to_snp_[i]] = weights_[i];

    for (int i = 0; i < num_snp_; i++) mafvec_.push_back((num_snp == 1) ? 0.5 : mafvec_dist(g_));
    for (int i = 0; i < num_snp_; i++) chrnumvec_.push_back(1);
  }

  std::vector<int>* tag_to_snp() { return &tag_to_snp_; }
  std::vector<float>* zvec() { return &z_vec_; }
  std::vector<float>* nvec() { return &n_vec_; }
  std::vector<float>* mafvec() { return &mafvec_; }
  std::vector<int>* chrnumvec() { return &chrnumvec_; }
  std::vector<float>* weights() { return &weights_; }

  // functions with complete=true return logically the same data, but suitable for use_complete_tag_indices representation.
  std::vector<int>* tag_to_snp(bool complete) { return complete ? &tag_to_snp_complete_ : &tag_to_snp_; }
  std::vector<float>* zvec(bool complete) { return complete ? &z_vec_complete_ : &z_vec_; }
  std::vector<float>* nvec(bool complete) { return complete ? &n_vec_complete_ : &n_vec_; }
  std::vector<float>* weights(bool complete) { return complete ? &weights_complete_ : &weights_; }

  void regenerate_zvec() {
    z_vec_.clear();
    z_vec_complete_.clear();
    std::normal_distribution<float> norm_dist(0.0, 1.5);
    for (int i = 0; i < num_tag_; i++) z_vec_.push_back((num_snp_ == 1) ? 2.0 : norm_dist(g_));
    for (int i = 0; i < num_snp_; i++) z_vec_complete_.push_back(NAN);
    for (int i = 0; i < num_tag_; i++) z_vec_complete_[tag_to_snp_[i]] = z_vec_[i];
  }

  // note that in this file tag_index usage remain inconsistent - tag_index variable may refer to values that go from 0 to num_snp-1 (i.e. runs across all SNPs, just just tag SNPs)
  void make_r2(int num_r2, std::vector<int>* snp_index, std::vector<int>* tag_index, std::vector<float>* r2) {
    if (num_snp_ == 1) {
      return;
    }

    std::uniform_real_distribution<> dis2(0.0f, 1.0f);
    std::vector<float> rand_list;
    for (int i = 0; i < 1000; i++) rand_list.push_back(((i%3==0) ? -1.0f : 1.0f) * dis2(g_));
    
    std::set<std::tuple<int, int>> pairs;

    std::uniform_int_distribution<> dis(0, num_snp_ - 1);
    int num_try = 0;
    while ((tag_index->size() < num_r2) && (num_try < 10*num_r2)) {
      num_try++;
      int tag = dis(g_);
      int snp = dis(g_);
      if (tag <= snp) continue;
      if (pairs.find(std::make_tuple(tag, snp)) != pairs.end()) continue;
      tag_index->push_back(tag);
      snp_index->push_back(snp);
      pairs.insert(std::make_tuple(tag, snp));
      r2->push_back(rand_list[r2->size() % rand_list.size()]);
    }
  }

  std::mt19937& random_engine() { return g_; }
private:
  int num_snp_;
  int num_tag_;
  std::vector<int> tag_to_snp_;
  std::vector<int> tag_to_snp_complete_;
  std::vector<float> z_vec_;
  std::vector<float> z_vec_complete_;
  std::vector<float> mafvec_;
  std::vector<int> chrnumvec_;
  std::vector<float> n_vec_;
  std::vector<float> n_vec_complete_;
  std::vector<float> weights_;
  std::vector<float> weights_complete_;
  std::random_device rd_;
  std::mt19937 g_;
};

// --gtest_filter=LdTest.ValidateMultipleChromosomes
TEST(LdTest, ValidateMultipleChromosomes) {
  int num_snp = 60;
  int num_tag = 40;
  int kmax = 20; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  int chr_label = 1;
  TestMother tm(60, 40, N);
  
  std::vector<int> chrnumvec; 
  for (int i = 0; i < 40; i++) chrnumvec.push_back(1);
  for (int i = 0; i < 20; i++) chrnumvec.push_back(2);
  
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 1);
  calc.set_option("r2min", 0.15f);

  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &chrnumvec[0]);

  std::vector<int> snp_index_all, tag_index_all;
  std::vector<float> r2_all;
  tm.make_r2(200, &snp_index_all, &tag_index_all, &r2_all);

  for (int chr_label = 1; chr_label <= 2; chr_label++) {
    std::vector<int> snp_index, tag_index;
    std::vector<float> r2;
    int offset = (chr_label == 2) ? 40 : 0;
    for (int i = 0; i < r2_all.size(); i++) {
      if ((chr_label == 1) != (snp_index_all[i] < 40)) continue;
      snp_index.push_back(snp_index_all[i] - offset);
      tag_index.push_back(tag_index_all[i] - offset);
      r2.push_back(r2_all[i]);
    }

    calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  }

  calc.set_ld_r2_csr();  // finalize csr structure
  calc.set_option("diag", 1);

  // retrieve cached and non-cached LDr2 sum, and compare the result
  calc.set_option("cache_tag_r2sum", 0);
  std::vector<float> tag_r2_sum(num_tag*kmax, 0.0f);
  calc.retrieve_tag_r2_sum(0, 1, num_tag*kmax, &tag_r2_sum[0]);
  calc.set_option("cache_tag_r2sum", 1);
  std::vector<float> tag_r2_sum_cached(num_tag*kmax, 0.0f);
  calc.retrieve_tag_r2_sum(0, 1, num_tag*kmax, &tag_r2_sum_cached[0]);
  for (int i = 0; i < num_tag*kmax; i++) ASSERT_FLOAT_EQ(tag_r2_sum[i], tag_r2_sum_cached[i]);
}

void UgmgTest_CalcLikelihood(float r2min, int trait_index) {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 10;
  int num_tag = 5;
  int kmax = 20; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  int chr_label = 1;
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 1);
  calc.set_option("cache_tag_r2sum", 1);
  calc.set_option("r2min", r2min);

  calc.set_zvec(trait_index, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait_index, num_tag, &tm.nvec()->at(0));
  calc.set_weights(num_tag, &tm.weights()->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(20, &snp_index, &tag_index, &r2);
  
  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure

  if (true) {
    //calc.find_snp_order();

    // Now, let's calculate log likelihood.
    calc.find_tag_r2sum(0, 2.1);
    calc.find_tag_r2sum(0, 1.1);
    calc.find_tag_r2sum(0, 3.1);
    calc.find_tag_r2sum(0, 3);
    calc.find_tag_r2sum(0, 4);
    calc.find_tag_r2sum(0, 3.9);
  }

  double cost = calc.calc_univariate_cost(trait_index, 0.2, 1.2, 0.1);
  double cost_nocache = calc.calc_univariate_cost_nocache(trait_index, 0.2, 1.2, 0.1);
  ASSERT_TRUE(std::isfinite(cost));
  ASSERT_FLOAT_EQ(cost, cost_nocache);

  std::vector<float> zvec_grid, zvec_pdf, zvec_pdf_nocache;
  for (float z = 0; z < 15; z += 0.1) {
    zvec_grid.push_back(z);
    zvec_pdf.push_back(0.0f);
    zvec_pdf_nocache.push_back(0.0f);
  }

  const float sig2_beta = 0.1f;
  const float pi_val = 0.2f;
  const float sig2_zeroL = pi_val * sig2_beta;
  const float sig2_zeroA = 1.2f;
  const float sig2_zeroC = 1.0f;
  std::vector<float> pi_vec(num_snp, pi_val);
  std::vector<float> sig2_vec(num_snp, sig2_beta);

  calc.calc_univariate_pdf(trait_index, pi_val, sig2_zeroA, sig2_beta, zvec_grid.size(), &zvec_grid[0], &zvec_pdf[0]);
  calc.set_option("diag", 0.0);

  calc.set_option("cache_tag_r2sum", 0);
  calc.calc_univariate_pdf(trait_index, pi_val, sig2_zeroA, sig2_beta, zvec_grid.size(), &zvec_grid[0], &zvec_pdf_nocache[0]);

  for (int i = 0; i < zvec_pdf_nocache.size(); i++) {
    ASSERT_NEAR(zvec_pdf[i], zvec_pdf_nocache[i], 2e-7);  // 4.93722e-05 vs 4.9372218e-05 due to approximation of float as uint16_t
  }

  //int64_t BgmgCalculator::calc_univariate_power(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, float zthresh, int length, float* nvec, float* svec) {
  std::vector<float> nvec;
  for (int n = 10; n < 1000; n += 10) nvec.push_back(n);
  std::vector<float> svec(nvec.size(), 0.0f);
  float zthresh = 5.45f;

  calc.calc_univariate_power(trait_index, 0.2, 1.2, 0.1, zthresh, nvec.size(), &nvec[0], &svec[0]);
  for (int i = 1; i < svec.size(); i++) ASSERT_TRUE(svec[i] > svec[i-1]);
  if (r2min != 0) { ASSERT_NEAR(svec.front(), 5.36914813e-05, 1e-6); ASSERT_NEAR(svec.back(), 0.685004711, 1e-6); }
  else {            ASSERT_NEAR(svec.front(), 5.23356393e-05, 1e-6); ASSERT_NEAR(svec.back(), 0.682822824, 1e-6); }

  std::vector<float> c0(num_tag, 0.0), c1(num_tag, 0.0), c2(num_tag, 0.0);
  calc.calc_univariate_delta_posterior(trait_index, 0.2, 1.2, 0.1, num_tag, &c0[0], &c1[0], &c2[0]);
  for (int i = 0; i < num_tag; i++) {
    ASSERT_TRUE(c0[i] > 0);
    ASSERT_TRUE(c1[i] != 0);
    ASSERT_TRUE(c2[i] > 0);
    break;
  }

  calc.set_option("fast_cost", 1);
  cost = calc.calc_univariate_cost(trait_index, 0.2, 1.2, 0.1);
  ASSERT_TRUE(std::isfinite(cost));
}

// --gtest_filter=UgmgTest.CalcLikelihood
TEST(UgmgTest, CalcLikelihood) {
  const float r2min = 0.0; 
  const int trait_index = 2; // use second trait for calculations; should work...
  UgmgTest_CalcLikelihood(r2min, trait_index);
}

// --gtest_filter=UgmgTest.CalcLikelihood_with_r2min
TEST(UgmgTest, CalcLikelihood_with_r2min) {
  const float r2min = 0.2;
  const int trait_index = 1;
  UgmgTest_CalcLikelihood(r2min, trait_index);
}

void UgmgTest_CalcLikelihood_OneSnpRef() {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 1;
  int num_tag = 1;
  int num_tag_internal = 1;
  int kmax = 1; // #permutations
  int N = 10;  // gwas sample size, constant across all variants
  int chr_label = 1;
  const float r2min = 0;
  const bool use_complete_tag_indices = true;
  const int trait_index = 1;
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag_internal, &tm.tag_to_snp(use_complete_tag_indices)->at(0));
  calc.set_option("seed", 0);
  calc.set_option("kmax", kmax);
  calc.set_option("use_complete_tag_indices", use_complete_tag_indices);
  calc.set_option("threads", 1);

  calc.set_zvec(trait_index, num_tag_internal, &tm.zvec(use_complete_tag_indices)->at(0));
  calc.set_nvec(trait_index, num_tag_internal, &tm.nvec(use_complete_tag_indices)->at(0));
  calc.set_weights(num_tag_internal, &tm.weights(use_complete_tag_indices)->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  // tm.make_r2(0, &snp_index, &tag_index, &r2);
  
  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure
  
  const float pi_val = 1.0f;
  const float sig2_beta = 0.1f;
  const float sig2_zeroL = 0.0f; //pi_val * sig2_beta;
  const float sig2_zeroA = 1.2f;
  const float sig2_zeroC = 1.0f;

  std::vector<float> pi_vec(num_snp, pi_val);
  std::vector<float> sig2_vec(num_snp, sig2_beta);
  std::vector<float> gradients(num_snp + 2, 0.0);

  calc.set_option("cost_calculator", 0);
  calc.set_option("aux_option", 4);
  double cost_unified_sampling = calc.calc_unified_univariate_cost_sampling(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA, sig2_zeroC, sig2_zeroL, &gradients[0], nullptr);
  const float gradient_analytic = gradients[num_snp];
  calc.set_option("aux_option", 0);
  
  const double sig2_zeroA_eps = 0.001;
  double cost_unified_sampling_plus = calc.calc_unified_univariate_cost_sampling(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA+sig2_zeroA_eps, sig2_zeroC, sig2_zeroL, nullptr, nullptr);
  double cost_unified_sampling_minus = calc.calc_unified_univariate_cost_sampling(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA-sig2_zeroA_eps, sig2_zeroC, sig2_zeroL, nullptr, nullptr);
  const float gradient_numeric = (-cost_unified_sampling_plus + cost_unified_sampling_minus)/(2*sig2_zeroA_eps);
  ASSERT_NEAR(gradient_analytic, gradient_numeric, 1e-3); 
}

// --gtest_filter=UgmgTest.CalcConvolveLikelihoodOneSnpRef
TEST(UgmgTest, CalcConvolveLikelihoodOneSnpRef) {
  UgmgTest_CalcLikelihood_OneSnpRef();
}

void UgmgTest_CalcLikelihood_testConvolution(float r2min, float z1max, int trait_index, float pi_val, double costvec[5], bool use_complete_tag_indices) {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 10;
  int num_tag = 5;
  int num_tag_internal = (use_complete_tag_indices ? num_snp : num_tag);
  int kmax = 20000; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  int chr_label = 1;
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag_internal, &tm.tag_to_snp(use_complete_tag_indices)->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 1);
  calc.set_option("cache_tag_r2sum", 1);
  calc.set_option("r2min", r2min);
  calc.set_option("use_complete_tag_indices", use_complete_tag_indices);
  calc.set_option("threads", 1);
  calc.set_option("z1max", z1max);
  calc.set_option("z2max", z1max); // set this for both traits

  calc.set_zvec(trait_index, num_tag_internal, &tm.zvec(use_complete_tag_indices)->at(0));
  calc.set_nvec(trait_index, num_tag_internal, &tm.nvec(use_complete_tag_indices)->at(0));
  calc.set_weights(num_tag_internal, &tm.weights(use_complete_tag_indices)->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(20, &snp_index, &tag_index, &r2);
  
  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure
  calc.set_option("diag", 0);
  
  const float sig2_beta = 0.1f;
  const float sig2_zeroL = pi_val * sig2_beta;
  const float sig2_zeroA = 1.2f;
  const float sig2_zeroC = 1.0f;

  calc.set_option("cost_calculator", 0);
  double cost_sampling = calc.calc_univariate_cost(trait_index, pi_val, sig2_zeroA, sig2_beta);
  calc.set_option("cost_calculator", 1);
  double cost_gaussian = calc.calc_univariate_cost(trait_index, pi_val, sig2_zeroA, sig2_beta);
  calc.set_option("cost_calculator", 2);
  double cost_convolve = calc.calc_univariate_cost(trait_index, pi_val, sig2_zeroA, sig2_beta);

  std::vector<float> pi_vec(num_snp, pi_val);
  std::vector<float> sig2_vec(num_snp, sig2_beta);
  std::vector<float> gradients(num_snp + 2, 0.0);
  double cost_unified_gaussian = calc.calc_unified_univariate_cost_gaussian(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA, sig2_zeroC, sig2_zeroL, nullptr);
  calc.set_option("aux_option", 4);
  double cost_unified_sampling = calc.calc_unified_univariate_cost_sampling(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA, sig2_zeroC, sig2_zeroL, &gradients[0], nullptr);
  calc.set_option("aux_option", 0);
  // std::cout << "gradients: ";  for (int i = 0; i < gradients.size(); i++) std::cout << gradients[i] << " "; std::cout << "\n";
  double cost_unified_smplfast = calc.calc_unified_univariate_cost_smplfast(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA, sig2_zeroC, sig2_zeroL, nullptr, nullptr);

  std::cout << std::setprecision(9) << cost_sampling << "(s), " << cost_gaussian << "(g), " << cost_convolve << "(c), " << cost_unified_gaussian << "(ug), " << cost_unified_sampling << "(us), " << cost_unified_smplfast << "(usf), " << std::endl;
  return;

  ASSERT_TRUE(std::isfinite(cost_sampling));
  ASSERT_TRUE(std::isfinite(cost_gaussian));
  ASSERT_TRUE(std::isfinite(cost_convolve));
  ASSERT_TRUE(std::isfinite(cost_unified_gaussian));
  ASSERT_TRUE(std::isfinite(cost_unified_sampling));  
  ASSERT_TRUE(std::isfinite(cost_unified_smplfast));  

  // compare that unified gaussian approximation gives the same answer as fast cost function
  ASSERT_NEAR(cost_gaussian, cost_unified_gaussian, 1e-5); 
  
  if (pi_val==1.0f) {
    ASSERT_FLOAT_EQ(cost_unified_sampling, cost_unified_gaussian); 
  }

  ASSERT_FLOAT_EQ(costvec[0], cost_sampling);
  ASSERT_FLOAT_EQ(costvec[1], cost_gaussian);
  ASSERT_FLOAT_EQ(costvec[2], cost_convolve);
  ASSERT_FLOAT_EQ(costvec[3], cost_unified_gaussian);
  ASSERT_FLOAT_EQ(costvec[4], cost_unified_sampling);
  ASSERT_FLOAT_EQ(costvec[5], cost_unified_smplfast);

  std::vector<float> zvec_grid, zvec_pdf_unified;
  for (float z = 0; z < 15; z += 0.1) {
    zvec_grid.push_back(z);
    zvec_pdf_unified.push_back(0.0f);
  }
  calc.calc_unified_univariate_pdf(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA, sig2_zeroC, sig2_zeroL, zvec_grid.size(), &zvec_grid[0], &zvec_pdf_unified[0]);
  for (int i = 0; i < zvec_pdf_unified.size(); i++) {
    ASSERT_TRUE(std::isfinite(zvec_pdf_unified[i]));
  }

  std::vector<float> nvec;
  for (int n = 10; n < 1000; n += 10) nvec.push_back(n);
  std::vector<float> svec_num_unified(nvec.size(), 0.0f);
  std::vector<float> svec_denom_unified(nvec.size(), 0.0f);
  float zthresh = 5.45f;

  if (pi_val==1.0f) calc.set_option("kmax", 1);

  calc.calc_unified_univariate_power(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA, sig2_zeroC, sig2_zeroL, zthresh, nvec.size(), &nvec[0], &svec_num_unified[0], &svec_denom_unified[0]);
  for (int i = 1; i < nvec.size(); i++) ASSERT_TRUE(svec_num_unified[i]/svec_denom_unified[i] > svec_num_unified[i-1]/svec_denom_unified[i-1]);

  if (pi_val != 1.0f) {
    std::vector<float> svec(nvec.size(), 0.0f);
    calc.calc_univariate_power(trait_index, pi_val, sig2_zeroA, sig2_beta, zthresh, nvec.size(), &nvec[0], &svec[0]);
    for (int i = 1; i < svec.size(); i++) ASSERT_TRUE(svec[i] > svec[i-1]);
    //for (int i = 0; i < svec.size(); i++) std::cout << svec[i] << "\t" << svec_unified[i] << std::endl;
  }

  std::vector<float> c0_unified(num_tag_internal, 0.0), c1_unified(num_tag_internal, 0.0), c2_unified(num_tag_internal, 0.0);
  calc.calc_unified_univariate_delta_posterior(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA, sig2_zeroC, sig2_zeroL, num_tag_internal, &c0_unified[0], &c1_unified[0], &c2_unified[0]);
  for (int i = 0; i < num_tag_internal; i++) {
    if (tm.weights(use_complete_tag_indices)->at(i) == 0) continue;
    ASSERT_TRUE(c0_unified[i] > 0); ASSERT_TRUE(c1_unified[i] != 0); ASSERT_TRUE(c2_unified[i] > 0);
    //std::cout << c0_unified[i] << "\t" << c1_unified[i] << "\t" << c2_unified[i] << std::endl;
  }
  
  if (pi_val != 1.0f) {
    std::vector<float> c0(num_tag_internal, 0.0), c1(num_tag_internal, 0.0), c2(num_tag_internal, 0.0);
    calc.calc_univariate_delta_posterior(trait_index, pi_val, sig2_zeroA, sig2_beta, num_tag_internal, &c0[0], &c1[0], &c2[0]);
    for (int i = 0; i < num_tag_internal; i++) {
      if (tm.weights(use_complete_tag_indices)->at(i) == 0) continue;
      ASSERT_TRUE(c0[i] > 0); ASSERT_TRUE(c1[i] != 0); ASSERT_TRUE(c2[i] > 0);
      //std::cout << c0_unified[i] << "\t" << c1_unified[i] << "\t" << c2_unified[i] << "\t\t" << c0[i] << "\t" << c1[i] << "\t" << c2[i] << std::endl;
    }
  }
}

double calcLikelihoodUnifiedGaussian(float r2min, int trait_index, bool use_complete_tag_indices, float pi_val) {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 10;
  int num_tag = 5;
  int num_tag_internal = (use_complete_tag_indices ? num_snp : num_tag);
  int N = 100;  // gwas sample size, constant across all variants
  int chr_label = 1;
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag_internal, &tm.tag_to_snp(use_complete_tag_indices)->at(0));
  calc.set_option("seed", 0);
  calc.set_option("r2min", r2min);
  calc.set_option("use_complete_tag_indices", use_complete_tag_indices);
  calc.set_option("threads", 1);

  calc.set_zvec(trait_index, num_tag_internal, &tm.zvec(use_complete_tag_indices)->at(0));
  calc.set_nvec(trait_index, num_tag_internal, &tm.nvec(use_complete_tag_indices)->at(0));
  calc.set_weights(num_tag_internal, &tm.weights(use_complete_tag_indices)->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(5, &snp_index, &tag_index, &r2);
  
  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure

  const float sig2_beta = 0.1f;
  const float sig2_zeroL = pi_val * sig2_beta;
  const float sig2_zeroA = 1.2f;
  const float sig2_zeroC = 1.0f;

  std::vector<float> pi_vec(num_snp, pi_val);
  std::vector<float> sig2_vec(num_snp, sig2_beta);
  return calc.calc_unified_univariate_cost_gaussian(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA, sig2_zeroC, sig2_zeroL, nullptr);
}

// --gtest_filter=UgmgTest.CalcConvolveLikelihood
TEST(UgmgTest, CalcConvolveLikelihood) {
  const float r2min = 0.0; const float z1max = 100000;
  const int trait_index = 2; // use second trait for calculations; should work...
  //double costvec[5] = {16.0114786, 15.8589964, 15.9299297, 15.8589949, 15.9333976};
  double costvec[6] = {7.73306899, 7.68907727, 7.69623015, 7.68907783, 7.69566212, 7.69951471};
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 0.2f, costvec, true);
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 0.2f, costvec, false);
}

// --gtest_filter=UgmgTest.CalcConvolveLikelihood_z1max
TEST(UgmgTest, CalcConvolveLikelihood_z1max) {
  const float r2min = 0.0; const float z1max = 1.2;
  const int trait_index = 2; // use second trait for calculations; should work...
  //double costvec[5] = {13.0109558, 12.789616, 12.9043533, 12.789611, 12.9075708};
  double costvec[6] = {6.25358767, 6.16207771, 6.20770748, 6.16207692, 6.20474518, 6.20653857 };
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 0.2f, costvec, true);
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 0.2f, costvec, false);
}

// --gtest_filter=UgmgTest.CalcConvolveLikelihoodInft
TEST(UgmgTest, CalcConvolveLikelihoodInft) {
  const float r2min = 0.0; const float z1max = 100000;
  const int trait_index = 2; // use second trait for calculations; should work...
  //double costvec[5] = {1e+100, 20.5189491, 20.5189537, 20.518953, 20.5189529};
  double costvec[6] = { 1e+100, 9.93359483, 9.93360452, 9.93360435, 9.93360429, 9.93360429} ;
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 1.0f, costvec, true);
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 1.0f, costvec, false);
}

// --gtest_filter=UgmgTest.CalcConvolveLikelihoodInft_z1max
TEST(UgmgTest, CalcConvolveLikelihoodInft_z1max) {
  const float r2min = 0.0; const float z1max = 1.2;
  const int trait_index = 2; // use second trait for calculations; should work...
  //double costvec[5] = {1e+100, 16.9341954, 16.9341894, 16.9341892, 16.934189};
  double costvec[6] = {1e+100, 8.05745833, 8.05746337, 8.05746326, 8.05746319, 8.0574632};
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 1.0f, costvec, true);
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 1.0f, costvec, false);
}

// --gtest_filter=UgmgTest.CalcConvolveLikelihood_with_r2min
TEST(UgmgTest, CalcConvolveLikelihood_with_r2min) {
  const float r2min = 0.2; const float z1max = 100000;
  const int trait_index = 1;
  //double costvec[5] = {16.00285, 15.840247, 15.9186396, 15.8402427, 15.925907};  
  double costvec[6] = {7.73092555, 7.72741705, 7.71021994, 7.72741465, 7.71344399, 7.71267116};
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 0.2f, costvec, true);
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 0.2f, costvec, false);
}

// --gtest_filter=UgmgTest.CalcConvolveLikelihood_with_r2min_inft
TEST(UgmgTest, CalcConvolveLikelihood_with_r2min_inft) {
  const float r2min = 0.2; const float z1max = 100000;
  const int trait_index = 1;
  //double costvec[5] = {1e+100, 20.5189491, 20.5189472, 20.5189468, 20.5189467};  
  double costvec[6] = {1e+100, 9.93359483, 9.93359888, 9.9335989, 9.93359867, 9.93359858};
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 1.0f, costvec, true);
  UgmgTest_CalcLikelihood_testConvolution(r2min, z1max, trait_index, 1.0f, costvec, false);
}

// --gtest_filter=UgmgTest.CalcUnifiedGaussianLikelihood
TEST(UgmgTest, CalcUnifiedGaussianLikelihood) {
  float v1 = calcLikelihoodUnifiedGaussian(0.0, 1, true, 0.2f);
  float v2 = calcLikelihoodUnifiedGaussian(0.0, 1, false, 0.2f);
  ASSERT_FLOAT_EQ(v1, v2);

  v1 = calcLikelihoodUnifiedGaussian(0.0, 1, true, 1.0f);
  v2 = calcLikelihoodUnifiedGaussian(0.0, 1, false, 1.0f);
  ASSERT_FLOAT_EQ(v1, v2);

  v1 = calcLikelihoodUnifiedGaussian(0.2, 2, true, 0.2f);
  v2 = calcLikelihoodUnifiedGaussian(0.2, 2, false, 0.2f);
  ASSERT_FLOAT_EQ(v1, v2);

  v1 = calcLikelihoodUnifiedGaussian(0.2, 2, true, 1.0f);
  v2 = calcLikelihoodUnifiedGaussian(0.2, 2, false, 1.0f);
  ASSERT_FLOAT_EQ(v1, v2);
}

void BgmgTest_CalcLikelihood_testConvolution(float r2min, float z1max, float z2max, float* pi_vec, double costvec[5], bool use_complete_tag_indices) {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 10;
  int num_tag = 5;
  int num_tag_internal = (use_complete_tag_indices ? num_snp : num_tag);
  int kmax = 20000; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag_internal, &tm.tag_to_snp(use_complete_tag_indices)->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 3);
  calc.set_option("cache_tag_r2sum", 1);
  calc.set_option("r2min", r2min);
  calc.set_option("use_complete_tag_indices", use_complete_tag_indices);
  calc.set_option("threads", 1);
  calc.set_option("z1max", z1max);
  calc.set_option("z2max", z1max);

  int trait = 1;
  int chr_label = 1;
  calc.set_zvec(trait, num_tag_internal, &tm.zvec(use_complete_tag_indices)->at(0));
  calc.set_nvec(trait, num_tag_internal, &tm.nvec(use_complete_tag_indices)->at(0));

  trait = 2; tm.regenerate_zvec();
  calc.set_zvec(trait, num_tag_internal, &tm.zvec(use_complete_tag_indices)->at(0));
  calc.set_nvec(trait, num_tag_internal, &tm.nvec(use_complete_tag_indices)->at(0));

  calc.set_weights(num_tag_internal, &tm.weights(use_complete_tag_indices)->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(20, &snp_index, &tag_index, &r2);

  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure
  //calc.set_weights_randprune(20, 0.25, 0.0, false);
  calc.set_option("diag", 0);

  float sig2_beta[] = { 0.5, 0.3 };
  float rho_beta = 0.8;
  float sig2_zero[] = { 1.1, 1.2 };
  float rho_zero = 0.1;

  calc.set_option("cost_calculator", 0);
  double cost_sampling = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  calc.set_option("cost_calculator", 1);
  double cost_gaussian = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  calc.set_option("cost_calculator", 2);
  double cost_convolve = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);

  // num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, aux
  std::vector<float> pi_unified(3*num_snp, 0);
  for (int i = 0; i < num_snp; i++) {pi_unified[i] = pi_vec[0]; pi_unified[num_snp+i] = pi_vec[1];pi_unified[2*num_snp+i] = pi_vec[2]; }
  std::vector<float> sig2_unified(2*num_snp, 0);
  for (int i = 0; i < num_snp; i++) {sig2_unified[i] = sig2_beta[0]; sig2_unified[num_snp+i] = sig2_beta[1]; }  
  std::vector<float> rho_unified(num_snp, rho_beta);
  float sig2_zeroC[] = { 1.0, 1.0 };
  float sig2_zeroL[] = { (pi_vec[0] + pi_vec[2]) * sig2_beta[0], (pi_vec[1] + pi_vec[2]) * sig2_beta[1] };
  float rho_zeroL = rho_beta * pi_vec[2] / sqrt((pi_vec[0]+pi_vec[2]) * (pi_vec[1]+pi_vec[2]));
  double cost_gaussian_unified = calc.calc_unified_bivariate_cost_gaussian(num_snp, &pi_unified[0], &sig2_unified[0], &rho_unified[0], sig2_zero, sig2_zeroC, sig2_zeroL, rho_zero, rho_zeroL, nullptr);
  double cost_sampling_unified = calc.calc_unified_bivariate_cost_sampling(num_snp, &pi_unified[0], &sig2_unified[0], &rho_unified[0], sig2_zero, sig2_zeroC, sig2_zeroL, rho_zero, rho_zeroL, nullptr, nullptr);
  double cost_smplfast_unified = calc.calc_unified_bivariate_cost_smplfast(num_snp, &pi_unified[0], &sig2_unified[0], &rho_unified[0], sig2_zero, sig2_zeroC, sig2_zeroL, rho_zero, rho_zeroL, nullptr, nullptr);

  ASSERT_TRUE(std::isfinite(cost_sampling));
  ASSERT_TRUE(std::isfinite(cost_gaussian));
  ASSERT_TRUE(std::isfinite(cost_convolve));
  ASSERT_TRUE(std::isfinite(cost_gaussian_unified));  
  ASSERT_TRUE(std::isfinite(cost_sampling_unified));  
  ASSERT_TRUE(std::isfinite(cost_smplfast_unified));  

  // There are subtle details in the ways "pi" coefficients are treated in old sampling, new sampling and convolve
  // In the old sampling we sample exactly (N*pi) causal SNPs -- no variation here, while if we strictly follow the model this needs to follow binomial distributino
  // New sampling is asymptotically correct, and gives the same answer as convolve.
  // There is always a difference between cost_sampling and cost_sampling_unified due to the reasons highlighted above.
  std::cout << std::setprecision(9) << cost_sampling << "(s), " << cost_gaussian << "(g), " << cost_convolve << "(c), " << cost_gaussian_unified << "(ug), " << cost_sampling_unified << "(us), " << cost_smplfast_unified << "(usf), " << std::endl;

  ASSERT_FLOAT_EQ(costvec[0], cost_sampling);
  ASSERT_FLOAT_EQ(costvec[1], cost_gaussian);
  ASSERT_FLOAT_EQ(costvec[2], cost_convolve);
  ASSERT_FLOAT_EQ(costvec[3], cost_gaussian_unified);
  ASSERT_FLOAT_EQ(costvec[4], cost_sampling_unified);
  ASSERT_FLOAT_EQ(costvec[5], cost_smplfast_unified);

  if (pi_vec[2] == 1 && r2min == 0) {
    // can't validate this for r2min != 0, see "TBD: apply ld_tag_sum_r2_below_r2min as an infinitesimal model"  in BgmgCalculator::calc_bivariate_cost_fast.
    ASSERT_NEAR(cost_gaussian, cost_gaussian_unified, 1e-4);
  }

  std::vector<float> c00(num_tag_internal, 0.0), c10(num_tag_internal, 0.0), c01(num_tag_internal, 0.0), c20(num_tag_internal, 0.0), c11(num_tag_internal, 0.0), c02(num_tag_internal, 0.0);
  calc.calc_unified_bivariate_delta_posterior(num_snp, &pi_unified[0], &sig2_unified[0], &rho_unified[0], sig2_zero, sig2_zeroC, sig2_zeroL, rho_zero, rho_zeroL, num_tag_internal, &c00[0], &c10[0], &c01[0], &c20[0], &c11[0], &c02[0]);
  for (int i = 0; i < num_tag_internal; i++) {
    if (tm.weights(use_complete_tag_indices)->at(i) == 0) continue;
    ASSERT_TRUE(c00[i] > 0);
    ASSERT_TRUE(c10[i] != 0);
    ASSERT_TRUE(c01[i] != 0);
    ASSERT_TRUE(c20[i] > 0);
    ASSERT_TRUE(c11[i] != 0);
    ASSERT_TRUE(c02[i] > 0);
  }

  std::vector<float> zvec1_grid, zvec2_grid, zvec_pdf;
  for (float z1 = -10; z1 < 10; z1 += 0.2) {
    for (float z2 = -10; z2 < 10; z2 += 0.2) {
      zvec1_grid.push_back(z1);
      zvec2_grid.push_back(z2);
      zvec_pdf.push_back(0.0f);
    }
  }
 
  calc.set_option("kmax", 200);
  calc.calc_unified_bivariate_pdf(num_snp, &pi_unified[0], &sig2_unified[0], &rho_unified[0], sig2_zero, sig2_zeroC, sig2_zeroL, rho_zero, rho_zeroL, zvec_pdf.size(), &zvec1_grid[0], &zvec2_grid[0], &zvec_pdf[0]);
  for (int i = 0; i < num_tag_internal; i++) {
    if (tm.weights(use_complete_tag_indices)->at(i) == 0) continue;
    ASSERT_TRUE(c00[i] > 0);
    ASSERT_TRUE(c10[i] != 0);
    ASSERT_TRUE(c01[i] != 0);
    ASSERT_TRUE(c20[i] > 0);
    ASSERT_TRUE(c11[i] != 0);
    ASSERT_TRUE(c02[i] > 0);
  }

}

void BgmgTest_CalcLikelihood(float r2min) {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 10;
  int num_tag = 5;
  int kmax = 20; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 3);
  calc.set_option("cache_tag_r2sum", 1);
  calc.set_option("r2min", r2min);

  int trait = 1;
  int chr_label = 1;
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

  trait = 2; tm.regenerate_zvec();
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

  calc.set_weights(num_tag, &tm.weights()->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(20, &snp_index, &tag_index, &r2);

  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure
  calc.set_weights_randprune(20, 0.25, 0.0, false);

  float pi_vec[] = { 0.1, 0.2, 0.15 };
  float sig2_beta[] = { 0.5, 0.3 };
  float rho_beta = 0.8;
  float sig2_zero[] = { 1.1, 1.2 };
  float rho_zero = 0.1;

  double cost = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  double cost_nocache = calc.calc_bivariate_cost_nocache(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  ASSERT_TRUE(std::isfinite(cost));
  ASSERT_FLOAT_EQ(cost, cost_nocache);

  calc.set_option("diag", 0.0);
  calc.set_option("fast_cost", 1);
  cost = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  ASSERT_TRUE(std::isfinite(cost));

  std::vector<float> zvec1_grid, zvec2_grid, zvec_pdf, zvec_pdf_nocache;
  for (float z1 = -10; z1 < 10; z1 += 0.2) {
    for (float z2 = -10; z2 < 10; z2 += 0.2) {
      zvec1_grid.push_back(z1);
      zvec2_grid.push_back(z2);
      zvec_pdf.push_back(0.0f);
      zvec_pdf_nocache.push_back(0.0f);
    }
  }

  calc.calc_bivariate_pdf(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero, zvec_pdf.size(), &zvec1_grid[0], &zvec2_grid[0], &zvec_pdf[0]);

  std::vector<float> c00(num_tag, 0.0), c10(num_tag, 0.0), c01(num_tag, 0.0), c20(num_tag, 0.0), c11(num_tag, 0.0), c02(num_tag, 0.0);
  calc.calc_bivariate_delta_posterior(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero, num_tag, &c00[0], &c10[0], &c01[0], &c20[0], &c11[0], &c02[0]);
  for (int i = 0; i < num_tag; i++) {
    ASSERT_TRUE(c00[i] > 0);
    ASSERT_TRUE(c10[i] != 0);
    ASSERT_TRUE(c01[i] != 0);
    ASSERT_TRUE(c20[i] > 0);
    ASSERT_TRUE(c11[i] != 0);
    ASSERT_TRUE(c02[i] > 0);
    break;
  }

  calc.set_option("cache_tag_r2sum", 0);
  calc.calc_bivariate_pdf(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero, zvec_pdf.size(), &zvec1_grid[0], &zvec2_grid[0], &zvec_pdf_nocache[0]);

  for (int i = 0; i < zvec_pdf_nocache.size(); i++)
    ASSERT_FLOAT_EQ(zvec_pdf[i], zvec_pdf_nocache[i]);
}

// --gtest_filter=BgmgTest.CalcConvolveLikelihood_inft
TEST(BgmgTest, CalcConvolveLikelihood_inft) {
  const float r2min = 0.0; const float z1max = 10000; const float z2max = 10000;
  float pi_vec[] = { 0.0, 0.0, 1.0 }; 
  //double costvec[5] = {1e+100, 23.7376571, 23.7376779, 23.7376838, 23.7376843};
  double costvec[6] = {1e+100, 24.5043969, 24.5040033, 24.5043765, 24.5043765, 24.5043765};
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, true);
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, false);
}

// --gtest_filter=BgmgTest.CalcConvolveLikelihood_inft_z1max
TEST(BgmgTest, CalcConvolveLikelihood_inft_z1max) {
  const float r2min = 0.0; const float z1max = 1.2; const float z2max = 0.5;
  float pi_vec[] = { 0.0, 0.0, 1.0 }; 
  //double costvec[5] = {1e+100, 15.4232014, 15.4232295, 15.4232234, 15.4232239};
  double costvec[6] = {1e+100, 9.92858297, 9.92856951, 9.92855016, 9.92855022, 9.92855012 };
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, true);
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, false);
}

// --gtest_filter=BgmgTest.CalcConvolveLikelihood
TEST(BgmgTest, CalcConvolveLikelihood) {
  const float r2min = 0.0; const float z1max = 10000; const float z2max = 10000;
  float pi_vec[] = { 0.1, 0.2, 0.15 };
  //double costvec[5] = {18.9650083, 18.0064621, 18.7262241, 20.3703903, 18.7258614};
  double costvec[6] = { 20.1751637, 19.5427144, 19.9417243, 21.1827832, 19.9471259, 19.9283836 };
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, true);
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, false);
}

// --gtest_filter=BgmgTest.CalcConvolveLikelihood_z1max
TEST(BgmgTest, CalcConvolveLikelihood_z1max) {
  const float r2min = 0.0; const float z1max = 1.2; const float z2max = 0.5;
  float pi_vec[] = { 0.1, 0.2, 0.15 };
  //double costvec[5] = {11.7194747, 10.7965121, 11.4819549, 13.1884063, 11.489059};
  double costvec[6] = {7.79574923, 7.28866735, 7.64816694, 8.63425681, 7.65108184, 7.63906013};
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, true);
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, false);
}

// --gtest_filter=BgmgTest.CalcConvolveLikelihood_with_r2min_inft
TEST(BgmgTest, CalcConvolveLikelihood_with_r2min_inft) {
  const float r2min = 0.2; const float z1max = 10000; const float z2max = 10000;
  float pi_vec[] = { 0.0, 0.0, 1.0 };
  //double costvec[5] = {1e+100, 23.5455545, 23.7376709, 23.7376766, 23.7376769};
  double costvec[6] = {1e+100, 24.2959385, 24.5039971, 24.5043688, 24.5043702, 24.5043707};
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, true);
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, false);
}

// --gtest_filter=BgmgTest.CalcConvolveLikelihood_with_r2min
TEST(BgmgTest, CalcConvolveLikelihood_with_r2min) {
  const float r2min = 0.2; const float z1max = 10000; const float z2max = 10000;
  float pi_vec[] = { 0.1, 0.2, 0.15 };
  //double costvec[5] = {18.9611203, 17.834824, 18.8095303, 20.3703838, 18.8040619};
  double costvec[6] = {20.1445922, 19.4030883, 20.0373171, 21.1827779, 20.0454545, 20.0298678};
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, true);
  BgmgTest_CalcLikelihood_testConvolution(r2min, z1max, z2max, pi_vec, costvec, false);
}

// bgmg-test.exe --gtest_filter=BgmgTest.CalcLikelihood
TEST(BgmgTest, CalcLikelihood) {
  const float r2min = 0.0;
  BgmgTest_CalcLikelihood(r2min);
}

// bgmg-test.exe --gtest_filter=BgmgTest.CalcLikelihood_with_r2min
TEST(BgmgTest, CalcLikelihood_with_r2min) {
  const float r2min = 0.20;
  BgmgTest_CalcLikelihood(r2min);
}


// bgmg-test.exe --gtest_filter=Test.RandomSeedAndThreading
TEST(Test, RandomSeedAndThreading) {
  int num_snp = 100;
  int num_tag = 50;
  int kmax = 200; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  int num_r2 = 20;
  TestMother tm(num_snp, num_tag, N);
  TestMother tm2(num_snp, num_tag, N);
  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(num_r2, &snp_index, &tag_index, &r2);

  int64_t seed_list[3] = { 123123123, 456456456, 123123123 };

  double ugmg_costs[3];
  double bgmg_costs[3];
  for (int num_threads = 1; num_threads <= 16; num_threads *= 2) {
    for (int i = 0; i < 3; i++) {
      BgmgCalculator calc;
      calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
      calc.set_option("max_causals", num_snp);
      calc.set_option("kmax", kmax);
      calc.set_option("num_components", 3);
      calc.set_option("cache_tag_r2sum", 1);
      calc.set_option("threads", num_threads);
      calc.set_seed(seed_list[i]);

      int num_threads_real;
#pragma omp parallel
      {
        num_threads_real = omp_get_num_threads();
      }
      ASSERT_EQ(num_threads_real, num_threads);

      int trait = 1;
      int chr_label = 1;
      calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
      calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));
      trait = 2;
      calc.set_zvec(trait, num_tag, &tm2.zvec()->at(0));
      calc.set_nvec(trait, num_tag, &tm2.nvec()->at(0));

      calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
      calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
      calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
      calc.set_ld_r2_csr();  // finalize csr structure
      calc.set_weights_randprune(20, 0.25, 0.0, false);

      float pi_vec[] = { 0.1, 0.2, 0.15 };
      float sig2_beta[] = { 0.5, 0.3 };
      float rho_beta = 0.8;
      float sig2_zero[] = { 1.1, 1.2 };
      float rho_zero = 0.1;
      int trait_index = 1;
      double ugmg_cost = calc.calc_univariate_cost(trait_index, pi_vec[0], sig2_zero[0], sig2_beta[0]);
      if (num_threads == 1) ugmg_costs[i] = ugmg_cost;
      else ASSERT_FLOAT_EQ(ugmg_costs[i], ugmg_cost);
      ASSERT_FLOAT_EQ(ugmg_cost, calc.calc_univariate_cost(trait_index, pi_vec[0], sig2_zero[0], sig2_beta[0]));  // check calc twice => the same cost
      ASSERT_FLOAT_EQ(ugmg_cost, calc.calc_univariate_cost_nocache(trait_index, pi_vec[0], sig2_zero[0], sig2_beta[0]));
      
      double bgmg_cost = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
      if (num_threads == 1) bgmg_costs[i] = bgmg_cost;
      else ASSERT_FLOAT_EQ(bgmg_costs[i], bgmg_cost);
      ASSERT_FLOAT_EQ(bgmg_costs[i], calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero));  // check calc twice => the same cost
      ASSERT_FLOAT_EQ(bgmg_costs[i], calc.calc_bivariate_cost_nocache(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero));

      if ((num_threads == 1) && (i==0)) {
        // test that pdf calculation has the same logic as log likelihood calculation
        std::vector<float> weights(num_tag, 0.0);
        calc.retrieve_weights(num_tag, &weights[0]);

        for (int test_caching = 0; test_caching < 2; test_caching++) {
          calc.set_option("cache_tag_r2sum", test_caching == 0);
          std::vector<float> ugmg_pdf(num_tag, 0.0);
          std::vector<float> bgmg_pdf(num_tag, 0.0);

          for (int tag_pdf_index = 0; tag_pdf_index < num_tag; tag_pdf_index++) {
            std::vector<float> weights2(num_tag, 0.0);
            weights2[tag_pdf_index] = 1.0;
            calc.set_weights(num_tag, &weights2[0]);
            calc.calc_univariate_pdf(trait_index, pi_vec[0], sig2_zero[0], sig2_beta[0], 1, &tm.zvec()->at(tag_pdf_index), &ugmg_pdf[tag_pdf_index]);
            calc.calc_bivariate_pdf(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero, 1, &tm.zvec()->at(tag_pdf_index), &tm2.zvec()->at(tag_pdf_index), &bgmg_pdf[tag_pdf_index]);
          }
          calc.set_weights(num_tag, &weights[0]); // restore weights back to the original statse

          double ugmg_cost_from_pdf = 0.0, bgmg_cost_from_pdf = 0.0;
          for (int tag_index = 0; tag_index < num_tag; tag_index++) {
            if (weights[tag_index] == 0) continue;
            ugmg_cost_from_pdf += -std::log(static_cast<double>(ugmg_pdf[tag_index])) * weights[tag_index];
            bgmg_cost_from_pdf += -std::log(static_cast<double>(bgmg_pdf[tag_index])) * weights[tag_index];
          }

          ASSERT_FLOAT_EQ(ugmg_cost, ugmg_cost_from_pdf);
          ASSERT_FLOAT_EQ(bgmg_cost, bgmg_cost_from_pdf);
        }
      }
    }

    ASSERT_FLOAT_EQ(ugmg_costs[0], ugmg_costs[2]);
    ASSERT_TRUE(abs(ugmg_costs[0] - ugmg_costs[1]) > 1e-4);
    ASSERT_FLOAT_EQ(bgmg_costs[0], bgmg_costs[2]);
    ASSERT_TRUE(abs(bgmg_costs[0] - bgmg_costs[1]) > 1e-4);
  }
}

void test_tag_r2_caching() {
  int trait_index = 1;
  int num_snp = 100;
  int num_tag = 50;
  int kmax = 200; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  int num_r2 = 20;
  TestMother tm(num_snp, num_tag, N);
  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(num_r2, &snp_index, &tag_index, &r2);
  
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 1);
  calc.set_option("cache_tag_r2sum", 1);
  calc.set_seed(123123123);
  int trait = 1;
  int chr_label = 1;
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

  int sequence_length = 10;
  std::vector<float> num_causal_sequence;
  std::uniform_real_distribution<float> rng(51.0, 99.5);
  for (int i = 0; i < sequence_length; i++) num_causal_sequence.push_back(rng(tm.random_engine()));
  
  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure
  calc.set_weights_randprune(20, 0.25, 0.0, false);

  int repeats_count = 100;
  std::vector<double> costs(sequence_length, 0.0);
  for (int j = 0; j < repeats_count; j++) {  // repeat the sequence 100 times and validate that we got the same cost.
    for (int i = 0; i < sequence_length; i++) {
      float pi = static_cast<float>(num_causal_sequence[i]) / static_cast<float>(num_snp);
      double cost = calc.calc_univariate_cost(trait_index, pi, 0.2, 0.15);
      if (j == 0) costs[i] = cost;
      ASSERT_FLOAT_EQ(cost, costs[i]);
    }
  }
}

// bgmg-test.exe --gtest_filter=Test.tag_r2_caching
TEST(Test, tag_r2_caching) {
  test_tag_r2_caching();
}

// --gtest_filter=Test.perform_ld_clump
TEST(Test, perform_ld_clump) {
  int num_snp = 100;
  int num_tag = 50;
  int N = 100; 
  int num_r2 = 2500;
  int chr_label = 1;
  TestMother tm(num_snp, num_tag, N);
  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(num_r2, &snp_index, &tag_index, &r2);
  
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure
  calc.set_option("diag", 0);
  
  std::vector<float> buffer;
  std::uniform_real_distribution<float> rng(0.0f, 1.0f);
  for (int i = 0; i < num_tag; i++) buffer.push_back(rng(tm.random_engine()));
  // for (int i = 0; i < num_tag; i++) std::cout << buffer[i] << " "; std::cout << std::endl;
  calc.perform_ld_clump(0.6, num_tag, &buffer[0]);
  // for (int i = 0; i < num_tag; i++) std::cout << buffer[i] << " "; std::cout << std::endl;

  int pass_clump = 0;
  for (int i = 0; i < num_tag; i++) if (std::isfinite(buffer[i])) pass_clump ++;
  ASSERT_EQ(pass_clump, 18); // hardcode how many pass clumping on this specific dataset
}

// --gtest_filter=Test.performance
TEST(Test, performance) {
  return;
  int trait_index = 1;
  // ideally test with scale = 100. To speedup, use 10 or 1.
  for (int scale = 1; scale <= 100; scale *= 10) {
    SimpleTimer timer_prep(-1);
    std::cout << "scale      : " << scale << "\n";
    int num_snp = 100000 * scale; 
    int num_tag = 10000 * scale;
    int kmax = 10 * scale; // #permutations
    int N = 100000;  // gwas sample size, constant across all variants
    int num_r2 = 10000000 * scale;
    TestMother tm(num_snp, num_tag, N);
    std::vector<int> snp_index, tag_index;
    std::vector<float> r2;
    tm.make_r2(num_r2, &snp_index, &tag_index, &r2);

    BgmgCalculator calc;
    calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
    calc.set_option("max_causals", 0.02 * static_cast<float>(num_snp));
    calc.set_option("kmax", kmax);
    calc.set_option("num_components", 1);
    calc.set_option("cache_tag_r2sum", 1);
    calc.set_option("threads", 32);
    calc.set_seed(123123123);
    int trait = 1;
    int chr_label = 1;
    calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
    calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

    calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
    calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
    calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
    calc.set_ld_r2_csr();  // finalize csr structure
    calc.set_weights_randprune(20, 0.25, 0.0, false);

    std::cout << "preparation: " << timer_prep.elapsed_ms() << " ms\n";
    float pivec[5] = { 0.0001, 0.0003, 0.001, 0.003, 0.01 };
    for (int repeat = 0; repeat < 5; repeat++) {
      SimpleTimer timer(-1);
      double cost_float = calc.calc_univariate_cost_nocache_float(trait_index, pivec[repeat], 0.2, 0.15);
      int time_with_float = timer.elapsed_ms();
      std::cout << "float  cost: " << cost_float << ", time: " << time_with_float << " ms\n";
      ASSERT_TRUE(std::isfinite(cost_float));

      SimpleTimer timer2(-1);
      double cost_double = calc.calc_univariate_cost_nocache_double(trait_index, pivec[repeat], 0.2, 0.15);
      int time_with_double = timer2.elapsed_ms();
      std::cout << "double cost: " << cost_double << ", time: " << time_with_double << " ms\n";
      ASSERT_TRUE(std::isfinite(cost_double));

      std::cout << "cost diff  : " << cost_double - cost_float << "\n";
    }
  }
}

void BgmgTest_CalcSampling_testPerformance(float r2min, float z1max, float z2max, float* pi_vec, int type) {
  int num_snp = 1000;
  int num_tag = 100;
  bool use_complete_tag_indices = false;
  int num_tag_internal = (use_complete_tag_indices ? num_snp : num_tag);
  int kmax = 20000; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag_internal, &tm.tag_to_snp(use_complete_tag_indices)->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 3);
  calc.set_option("cache_tag_r2sum", 0);
  calc.set_option("r2min", r2min);
  calc.set_option("use_complete_tag_indices", use_complete_tag_indices);
  calc.set_option("threads", 6);
  calc.set_option("z1max", z1max);
  calc.set_option("z2max", z1max);

  int trait = 1;
  int chr_label = 1;
  calc.set_zvec(trait, num_tag_internal, &tm.zvec(use_complete_tag_indices)->at(0));
  calc.set_nvec(trait, num_tag_internal, &tm.nvec(use_complete_tag_indices)->at(0));

  trait = 2; tm.regenerate_zvec();
  calc.set_zvec(trait, num_tag_internal, &tm.zvec(use_complete_tag_indices)->at(0));
  calc.set_nvec(trait, num_tag_internal, &tm.nvec(use_complete_tag_indices)->at(0));

  calc.set_weights(num_tag_internal, &tm.weights(use_complete_tag_indices)->at(0));

  // dense LD r2 matrix, all in perfect LD
  size_t num_r2 = num_snp * (num_snp-1) / 2;
  std::vector<int> snp_index(num_r2, 0), tag_index(num_r2, 0);
  std::vector<float> r2(num_r2, 1.0f);
  int index_r2 = 0;
  for (int i = 0; i < num_snp; i++) {
    for (int j = (i+1); j < num_snp; j++) {
      snp_index[index_r2] = i;
      tag_index[index_r2] = j;
      index_r2++;
    }
  }
    
  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure
  //calc.set_weights_randprune(20, 0.25, 0.0, false);
  calc.set_option("diag", 0);

  float sig2_beta[] = { 0.5, 0.3 };
  float rho_beta = 0.8;
  float sig2_zero[] = { 1.1, 1.2 };
  float rho_zero = 0.1;

  // num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, aux
  std::vector<float> pi_unified(3*num_snp, 0);
  for (int i = 0; i < num_snp; i++) {pi_unified[i] = pi_vec[0]; pi_unified[num_snp+i] = pi_vec[1];pi_unified[2*num_snp+i] = pi_vec[2]; }
  std::vector<float> sig2_unified(2*num_snp, 0);
  for (int i = 0; i < num_snp; i++) {sig2_unified[i] = sig2_beta[0]; sig2_unified[num_snp+i] = sig2_beta[1]; }  
  std::vector<float> rho_unified(num_snp, rho_beta);
  float sig2_zeroC[] = { 1.0, 1.0 };
  float sig2_zeroL[] = { (pi_vec[0] + pi_vec[2]) * sig2_beta[0], (pi_vec[1] + pi_vec[2]) * sig2_beta[1] };
  float rho_zeroL = rho_beta * pi_vec[2] / sqrt((pi_vec[0]+pi_vec[2]) * (pi_vec[1]+pi_vec[2]));
  
  const int nrep = 1;
  for (int repi = 0; repi < nrep; repi++) {
    double cost;
    if (type == 0) cost = calc.calc_unified_bivariate_cost_sampling(num_snp, &pi_unified[0], &sig2_unified[0], &rho_unified[0], sig2_zero, sig2_zeroC, sig2_zeroL, rho_zero, rho_zeroL, nullptr, nullptr);
    if (type == 1) cost = calc.calc_unified_bivariate_cost_smplfast(num_snp, &pi_unified[0], &sig2_unified[0], &rho_unified[0], sig2_zero, sig2_zeroC, sig2_zeroL, rho_zero, rho_zeroL, nullptr, nullptr);
    if (type == 2) { calc.set_option("cost_calculator", 0); cost = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero); }
    ASSERT_TRUE(std::isfinite(cost));  
  }
}

// --gtest_filter=BgmgTest.CalcSamplingPerformanceSampling
TEST(BgmgTest, CalcSamplingPerformanceSampling) {
  const float r2min = 0.0; const float z1max = 10000; const float z2max = 10000;
  float pi_vec[] = { 0.2, 0.1, 0.15 }; 
  BgmgTest_CalcSampling_testPerformance(r2min, z1max, z2max, pi_vec, 0);
}

// --gtest_filter=BgmgTest.CalcSamplingPerformanceSmplfast
TEST(BgmgTest, CalcSamplingPerformanceSmplfast) {
  const float r2min = 0.0; const float z1max = 10000; const float z2max = 10000;
  float pi_vec[] = { 0.2, 0.1, 0.15 }; 
  BgmgTest_CalcSampling_testPerformance(r2min, z1max, z2max, pi_vec, 1);
}

// --gtest_filter=BgmgTest.CalcSamplingPerformanceLegacy
TEST(BgmgTest, CalcSamplingPerformanceLegacy) {
  const float r2min = 0.0; const float z1max = 10000; const float z2max = 10000;
  float pi_vec[] = { 0.2, 0.1, 0.15 }; 
  BgmgTest_CalcSampling_testPerformance(r2min, z1max, z2max, pi_vec, 2);
}


}  // namespace
