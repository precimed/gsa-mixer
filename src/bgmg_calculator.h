/*
  bgmg - tool to calculate log likelihood of BGMG and UGMG mixture models
  Copyright (C) 2018 Oleksandr Frei 

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <stdint.h>

#include <unordered_map>
#include <memory>

#include <iostream>
#include <sstream>
#include <vector>
#include <valarray>

#include "boost/utility.hpp"

#include "bgmg_log.h"
#include "bgmg_parse.h"
#include "bgmg_serialize.h"
#include "ld_matrix_csr.h"

enum CostCalculator {
  CostCalculator_Sampling = 0,
  CostCalculator_Gaussian = 1,
  CostCalculator_Convolve = 2,
  CostCalculator_Smplfast = 3,
  CostCalculator_MAX = 4,
}; 

enum AuxOption {
  AuxOption_None = 0,
  AuxOption_Ezvec2 = 1,
  AuxOption_TagPdf = 2,
  AuxOption_TagPdfErr = 3,
  AuxOption_Gradients = 4,
  AuxOption_MAX = 5,
};

#define DEFAULT_ZMAX 1e10

class MultinomialSampler;

// Singleton class to manage a collection of objects, identifiable with some integer ID.
template<class Type>
class TemplateManager : boost::noncopyable {
 public:
  static TemplateManager<Type>& singleton() {
    static TemplateManager<Type> manager;
    return manager;
  }

  std::shared_ptr<Type> Get(int id) {
    if (map_.find(id) == map_.end()) {
      LOG << " Create new context (id=" << id << ")";
      map_[id] = std::make_shared<Type>();
    }
    return map_[id];
  }

  void Erase(int id) {
    LOG << " Dispose context (id=" << id << ")";
    map_.erase(id);
  }

  void Clear() {
    LOG << " Dispose all context ids";
    map_.clear();
  }

 private:
  TemplateManager() { }  // Singleton (make constructor private)
  std::unordered_map<int, std::shared_ptr<Type>> map_;
};

template<typename T>
class DenseMatrix {
 public:
  explicit DenseMatrix(size_t no_rows = 0, size_t no_columns = 0, bool store_by_rows = true)
    : no_rows_(no_rows),
    no_columns_(no_columns),
    store_by_rows_(store_by_rows),
    data_(nullptr) {
    if (no_rows > 0 && no_columns > 0) {
      data_ = new T[no_rows_ * no_columns_];
    }
  }

  DenseMatrix(const DenseMatrix<T>& src_matrix) {
    no_rows_ = src_matrix.no_rows();
    no_columns_ = src_matrix.no_columns();
    store_by_rows_ = src_matrix.store_by_rows_;
    if (no_columns_ >0 && no_rows_ > 0) {
      data_ = new T[no_rows_ * no_columns_];

      for (size_t i = 0; i < no_rows_ * no_columns_; ++i) {
        data_[i] = src_matrix.get_data()[i];
      }
    } else {
      data_ = nullptr;
    }
  }

  virtual ~DenseMatrix() {
    delete[] data_;
  }

  void InitializeZeros() {
    memset(data_, 0, sizeof(T)* no_rows_ * no_columns_);
  }

  T& operator() (size_t index_row, size_t index_col) {
    if (store_by_rows_) {
      return data_[index_row * no_columns_ + index_col];
    }
    return data_[index_col * no_rows_ + index_row];
  }

  const T& operator() (size_t index_row, size_t index_col) const {
    assert(index_row < no_rows_);
    assert(index_col < no_columns_);
    if (store_by_rows_) {
      return data_[index_row * no_columns_ + index_col];
    }
    return data_[index_col * no_rows_ + index_row];
  }

  DenseMatrix<T>& operator= (const DenseMatrix<T>& src_matrix) {
    no_rows_ = src_matrix.no_rows();
    no_columns_ = src_matrix.no_columns();
    store_by_rows_ = src_matrix.store_by_rows_;
    if (data_ != nullptr) {
      delete[] data_;
    }
    if (no_columns_ >0 && no_rows_ > 0) {
      data_ = new  T[no_rows_ * no_columns_];

      for (size_t i = 0; i < no_rows_ * no_columns_; ++i) {
        data_[i] = src_matrix.get_data()[i];
      }
    } else {
      data_ = nullptr;
    }

    return *this;
  }

  size_t no_rows() const { return no_rows_; }
  size_t no_columns() const { return no_columns_; }
  size_t size() const { return no_rows_ * no_columns_; }
  bool is_equal_size(const DenseMatrix<T>& rhs) const {
    return no_rows_ == rhs.no_rows_ && no_columns_ == rhs.no_columns_;
  }

  std::string to_str() {
    int rows_to_str = std::min<int>(5, no_rows_ - 1);
    int cols_to_str = std::min<int>(5, no_columns_ - 1);
    std::stringstream ss;
    ss << "[";
    for (int i = 0; i < rows_to_str; i++) {
      bool last_row = (i == (rows_to_str - 1));
      for (int j = 0; j < cols_to_str; j++) {
        bool last_col = (j == (cols_to_str - 1));
        ss << (*this)(i, j);
        if (last_col) ss << ", ...";
        else ss << ", ";
      }
      if (last_row) ss << "; ...";
      else ss << "; ";
    }
    ss << "]";

    size_t nnz = 0;
    for (int i = 0; i < no_rows_; i++)
      for (int j = 0; j < no_columns_; j++)
        if ((*this)(i, j) != 0) nnz++;
    ss << ", nnz=" << nnz;

    return ss.str();
  }


  T* get_data() {
    return data_;
  }

  const T* get_data() const {
    return data_;
  }

 private:
  size_t no_rows_;
  size_t no_columns_;
  bool store_by_rows_;
  T* data_;
};

class BgmgCalculator : public TagToSnpMapping {
 public:
  BgmgCalculator();
  virtual ~BgmgCalculator() {}

  int64_t init(std::string bim_file, std::string frq_file, std::string chr_labels, std::string trait1_file, std::string trait2_file, std::string exclude, std::string extract, std::string exclude_ranges);
  int64_t read_trait_file(int trait_index, std::string trait1_file, std::string exclude, std::string extract, std::string exclude_ranges);
  int64_t convert_plink_ld(std::string plink_ld_gz, std::string plink_ld_bin);  // require init() to be called first, e.i. doesn't work after set_tag_indices.

  // num_snp = total size of the reference (e.i. the total number of genotyped variants)
  // num_tag = number of tag variants to include in the inference (must be a subset of the reference)
  // indices = array of size num_tag, containing indices from 0 to num_snp-1
  // NB: all tag variants must have defined zvec, nvec, mafvec and weights.
  int64_t set_tag_indices(int num_snp, int num_tag, int* tag_indices);
  int64_t retrieve_tag_indices(int num_tag, int* tag_indices);

  int64_t set_chrnumvec(int num_snp, const int* chrlabel);
  int64_t retrieve_chrnumvec(int length, int* buffer);

  // consume input in plink format, e.i. lower triangular LD r2 matrix
  // - snp_index is less than tag_index;
  // - does not contain r2=1.0 of the variant with itself
  // must be called after set_tag_indices
  // must be called one for each chromosome, sequentially, starting from lowest chromosome number
  // non-tag variants will be ignored
  int64_t set_ld_r2_coo(int chr_label, int64_t length, int* snp_index, int* tag_index, float* r);
  int64_t set_ld_r2_coo(int chr_label, const std::string& filename);
  int64_t set_ld_r2_csr(int chr_label = -1);  // finalize

  int64_t num_ld_r_tag(int tag_index);
  int64_t retrieve_ld_r_tag(int tag_index, int length, int* snp_index, float* r);
  int64_t num_ld_r_chr(int chr_label);
  int64_t retrieve_ld_r_chr(int chr_label, int64_t length, int* tag_index, int* snp_index, float* r);
  int64_t num_ld_r_tag_range(int tag_index_from, int tag_index_to);
  int64_t retrieve_ld_r_tag_range(int tag_index_from, int tag_index_to, int length, int* tag_index, int* snp_index, float* r);

  // must be called after set_ld_r2, as it adjusts r2 matrix
  // one value for each snp (tag and non-tag)
  int64_t set_mafvec(int length, float* values);
  
  // zvec, nvec, weights for tag variants
  // all values must be defined
  int64_t set_zvec(int trait, int length, float* values);
  int64_t set_nvec(int trait, int length, float* values);
  int64_t set_causalbetavec(int trait, int length, float* values);
  int64_t set_weights(int length, float* values);
  int64_t set_weights_hardprune(float hardprune_subset, float hardprune_r2, float hardprune_maf, bool use_w_ld);
  int64_t set_weights_randprune(int n, float r2, float maf, bool use_w_ld);   // alternative to set_weights; calculates weights based on random pruning from LD matrix
  int64_t set_weights_randprune(int n, float r2, float maf, bool use_w_ld, std::string exclude, std::string extract);   // alternative to set_weights; calculates weights based on random pruning from LD matrix
  int64_t perform_ld_clump(float r2, int length, float* buffer);

  int64_t retrieve_zvec(int trait, int length, float* buffer);
  int64_t retrieve_nvec(int trait, int length, float* buffer);
  int64_t retrieve_causalbetavec(int trait, int length, float* buffer);
  int64_t retrieve_mafvec(int length, float* buffer);
  int64_t retrieve_weights(int length, float* buffer);

  int64_t set_option(char* option, double value);
  
  int64_t retrieve_ld_sum_r2(int length, float* buffer);
  int64_t retrieve_ld_sum_r4(int length, float* buffer);
  int64_t retrieve_fixed_effect_delta(int trait_index, int length, float* delta);
  
  void log_diagnostics();

  // a note about 'smplfast': this is a special implementation of sampling, valid when 'pi' is constant across SNPs, and it use chunks_forward_ (while all other functions in unified implementation use chunks_reverse_).
  double calc_unified_univariate_cost(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux);
  double calc_unified_univariate_cost_gaussian(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux);
  double calc_unified_univariate_cost_convolve(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux);
  double calc_unified_univariate_cost_sampling(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux, const float* weights);
  double calc_unified_univariate_cost_smplfast(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux, const float* weights);
  int64_t calc_unified_univariate_pdf(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, int length, float* zvec, float* pdf);
  int64_t calc_unified_univariate_power(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float zthresh, int length, float* nvec, float* svec_num, float* svec_denom);
  int64_t calc_unified_univariate_delta_posterior(int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, int length, float* c0, float* c1, float* c2);

  // pi_vec     : num_components X num_snp,  - weights of the three mixture components (num_components = 3)
  // sig2_vec   : num_traits X num_snp,      - variance of cauasal effects for the two traits (num_traits = 2)
  // rho_vec    : 1 x num_snps,              - correlation of genetic effects
  // sig2_zeroX : 1 x num_traits,            - inflation parameters (additive, multiplicative, truncation of the LD structure)
  // rho_zeroA, sig2_zeroL                   - correlation of sig2_zeroA and sig2_zeroL across the two traits
  // aux        : 3 x num_tag                - expore auxilary information, as defined by aux_option (for AuxOption_Ezvec2 aux will store three vectors: E[z20], E[z11], E[z02], in this order)
  double calc_unified_bivariate_cost(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux);
  double calc_unified_bivariate_cost_gaussian(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux);
  double calc_unified_bivariate_cost_convolve(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux);
  double calc_unified_bivariate_cost_sampling(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux, const float* weights);
  double calc_unified_bivariate_cost_smplfast(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux, const float* weights);
  int64_t calc_unified_bivariate_pdf(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, int length, float* zvec1, float* zvec2, float* pdf);
  int64_t calc_unified_bivariate_delta_posterior(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL,
                                                 int length, float* c00, float* c10, float* c01, float* c20, float* c11, float* c02);

  int64_t seed() { return seed_; }
  void set_seed(int64_t seed) { seed_ = seed; }

  virtual int num_snp() { return num_snp_; }
  virtual int num_tag() { return num_tag_; }
  virtual bool has_complete_tag_indices() { return num_snp_ == num_tag_; }
  virtual bool disable_snp_to_tag_map() { return disable_snp_to_tag_map_ || (num_snp_ == num_tag_); }
  virtual const std::vector<int>& tag_to_snp() { return  tag_to_snp_; }
  virtual const std::vector<int>& snp_to_tag() { return snp_to_tag_; }
  virtual const std::vector<char>& is_tag() { return is_tag_; }
  virtual const std::vector<int>& chrnumvec() { return chrnumvec_; }
  virtual const std::vector<float>& mafvec() { return mafvec_; }
  virtual std::vector<float>* mutable_mafvec() { return &mafvec_; }

  int64_t save(const char* file);
  int64_t load(const char* file);

 private:
  int num_snp_;
  int num_tag_;
  std::vector<int> tag_to_snp_; // 0..num_snp_-1, size=num_tag_
  std::vector<int> snp_to_tag_; // 0..num_tag-1,  size=num_snp_, -1 indicate non-tag snp
  std::vector<char> is_tag_;    // true or false, size=num_snp_, is_tag_[i] == (snp_to_tag_[i] != -1)
  std::vector<int> chrnumvec_;  // vector of chromosome labels, one per snp

  LdMatrixCsr ld_matrix_csr_;

  // all stored for for tag variants (only)
  std::vector<float> zvec1_;
  std::vector<float> nvec1_;
  std::vector<float> zvec2_;
  std::vector<float> nvec2_;
  std::vector<float> weights_;
  std::vector<float> mafvec_;

  std::vector<float> causalbetavec1_;  // assumed causal betas, added as fixed effects component in the model
  std::vector<float> causalbetavec2_;  // e.g. delta_j = \sqrt N_j \sum_i \sqrt H_i r_ij \beta_i   <- here "beta_i" is causalbetavec
  
  std::vector<float>* get_zvec(int trait_index);
  std::vector<float>* get_nvec(int trait_index);
  std::vector<float>* get_causalbetavec(int trait_index);

  // options, and what do they affect
  int k_max_;                      // number of samples used in: log-likelihood, posterior delta; also everywhere in legacy implementation
  int k_max_pdf_;                  // number of samples used in: qq-plot, power-curves
  int64_t seed_;
  bool use_complete_tag_indices_;  // an option that indicates that all SNPs are TAG (i.e. num_snp_ == num_tag_).
  bool disable_snp_to_tag_map_;    // save memory by disabling ld_matrix_csr_.chunks_forward_ (not needed in unified implementation)
  bool allow_ambiguous_snps_;      // NB, not serialized

  float r2_min_;
  float z1max_;
  float z2max_;
  CostCalculator cost_calculator_;
  double cubature_abs_error_;
  double cubature_rel_error_;
  int cubature_max_evals_;
  AuxOption aux_option_;       // controls auxilary data stored in "aux" array of calc_unified_univariate_cost_*** calls
  int ld_format_version_;      // overwrite format version for LD matrix files. Default -1. Set this to 0 to read from MiXeR v1.0 LD files.
  int retrieve_ld_sum_type_;   // control behaviour of retrieve_ld_sum_r2() and retrieve_ld_sum_r4()
                               // 0 = above r2min; 1 = below r2min; 2 = above r2min adjusted for hvec; 3 = below r2min adjusted for hvec; 
                               // For retrieve_ld_sum_r4() the "below r2min" option is not available, therefore 1 will work the same as 0, and 3 will work same as 2.

  void check_num_snp(int length);
  void check_num_tag(int length);

  void find_unified_univariate_tag_delta_sampling(int num_components, int k_max, float* pi_vec, float* sig2_vec, float sig2_zeroC, int tag_index, const float* nvec, const float* hvec, std::vector<float>* tag_delta2, MultinomialSampler* subset_sampler, LdMatrixRow* ld_matrix_row);
  void find_unified_bivariate_tag_delta_sampling(int num_snp, int k_max, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, int tag_index, const float* nvec1, const float* nvec2, const float* hvec, std::vector<float>* tag_delta20, std::vector<float>* tag_delta02, std::vector<float>* tag_delta11, MultinomialSampler* subset_sampler, LdMatrixRow* ld_matrix_row);

  void find_unified_univariate_tag_delta_smplfast(int num_components, int k_max, float* pi_vec, float* sig2_vec, float sig2_zeroC, int k_index, const float* nvec, const float* hvec, std::vector<float>* tag_delta2, MultinomialSampler* subset_sampler, LdMatrixRow* ld_matrix_row);
  void find_unified_bivariate_tag_delta_smplfast(int num_snp, int k_max, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, int k_index, const float* nvec1, const float* nvec2, const float* hvec, std::vector<float>* tag_delta20, std::vector<float>* tag_delta02, std::vector<float>* tag_delta11, MultinomialSampler* subset_sampler, LdMatrixRow* ld_matrix_row);

  void find_unified_univariate_gradients_sampling(int num_components, float* pi_vec, float* sig2_vec, float sig2_zeroC, int tag_index, const float* nvec, const float* hvec, std::valarray<float>* gradients,  MultinomialSampler* subset_sampler, LdMatrixRow* ld_matrix_row, const std::vector<float>& gradient_inc);

  void calc_fixed_effect_delta_from_causalbetavec(int trait_index, std::valarray<float>* delta);
  void find_z_minus_fixed_effect_delta(int trait_index, std::vector<float>* z_minus_fixed_effect_delta);

  int find_deftag_indices(const float* weights, std::vector<int>* deftag_indices);

  BimFile bim_file_;

  void dump(BgmgSerializer* serializer);  // save or load the object

  static void calc_bivariate_delta_posterior_integrals(float a, float b, float c, float i, float j, float k, float z1, float z2,
                                                       float* c00, float* c10, float* c01, float* c20, float* c11, float* c02);

  // ============= legacy stuff starts ============= 
 public:
  int64_t set_snp_order(int component_id, int64_t length, const int* buffer);
  int64_t retrieve_snp_order(int component_id, int64_t length, int* buffer);
  int64_t retrieve_k_pdf(int length, double* buffer);

  int64_t find_snp_order();  // private - only for testing
  int64_t find_tag_r2sum(int component_id, float num_causals);  // private - only for testing
  void find_tag_r2sum_no_cache(int component_id, float num_causal, int k_index, std::vector<float>* buffer); // private - only for testing
  void clear_tag_r2sum(int component_id);
  int64_t retrieve_tag_r2_sum(int component_id, float num_causal, int length, float* buffer);

  void clear_state();

  double calc_univariate_cost(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);
  double calc_univariate_cost_cache(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);
  double calc_univariate_cost_nocache(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);        // default precision (see FLOAT_TYPE in bgmg_calculator.cc)
  double calc_univariate_cost_nocache_float(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);  // for testing single vs double precision
  double calc_univariate_cost_nocache_double(int trait_index, float pi_vec, float sig2_zero, float sig2_beta); // for testing single vs double precision
  int64_t calc_univariate_pdf(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int length, float* zvec, float* pdf);
  int64_t calc_univariate_power(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, float zthresh, int length, float* nvec, float* svec);
  int64_t calc_univariate_delta_posterior(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, int length, float* c0, float* c1, float* c2);
  int64_t calc_bivariate_delta_posterior(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero,
                                         int length, float* c00, float* c10, float* c01, float* c20, float* c11, float* c02);

  double calc_bivariate_cost(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero);
  double calc_bivariate_cost_nocache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero);
  double calc_bivariate_cost_cache(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero);
  int64_t calc_bivariate_pdf(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero, int length, float* zvec1, float* zvec2, float* pdf);

 private:
  // vectors with one value for each component in the mixture
  // snp_order_ gives the order of how SNPs are considered to be causal 
  // tag_r2_sum_ gives cumulated r2 across causal SNPs, according to snp_order, where last_num_causals_ define the actual number of causal variants.
  int max_causals_;
  int num_components_;
  bool cache_tag_r2sum_;  
  bool calc_k_pdf_;            // a flag indicating whether we should calculate k_pdf_
  std::vector<std::shared_ptr<DenseMatrix<int>>> snp_order_;  // permutation matrix; #rows = pimax*num_snp; #cols=k_max_
  std::vector<std::shared_ptr<DenseMatrix<float>>> tag_r2sum_;
  std::vector<float>                               last_num_causals_;
  std::vector<double> k_pdf_;  // the log-likelihood cost calculated independently for each of 0...k_max-1 selections of causal variants.            

  int find_deftag_indices_znw(int trait_index, std::vector<int>* deftag_indices);
  int find_deftag_indices_znw(std::vector<int>* deftag_indices);
  int find_deftag_indices_nw(int trait_index, std::vector<int>* deftag_indices);
  int find_deftag_indices_nw(std::vector<int>* deftag_indices);
  int find_deftag_indices_w(std::vector<int>* deftag_indices);

  double calc_univariate_cost_fast(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);
  double calc_bivariate_cost_fast(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero);
  double calc_univariate_cost_convolve(int trait_index, float pi_vec, float sig2_zero, float sig2_beta);
  double calc_bivariate_cost_convolve(int pi_vec_len, float* pi_vec, int sig2_beta_len, float* sig2_beta, float rho_beta, int sig2_zero_len, float* sig2_zero, float rho_zero);

  template<typename T>
  friend double calc_univariate_cost_nocache_template(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, BgmgCalculator& rhs);
  // ============= legacy stuff ends ============= 
};

typedef TemplateManager<BgmgCalculator> BgmgCalculatorManager;
