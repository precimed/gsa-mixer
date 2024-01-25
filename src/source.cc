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

#include <stdexcept>

#include <bgmg.h>
#include "bgmg_calculator.h"
#include "annotate_intervals.h"
#include "ld_matrix.h"
#include "bgmg_log.h"


static std::string last_error_;
static const char* last_error() { return last_error_.c_str(); }
static void set_last_error(const std::string& error) { last_error_.assign(error); }
const char* bgmg_get_last_error() { return last_error(); }

static std::string last_buffer_;  // this could be boost::thread_specific_ptr<std::string>, but we don't need that
static std::string* last_buffer() {
  return &last_buffer_;
}

static char* StringAsArray(std::string* str) {
  return str->empty() ? NULL : &*str->begin();
}

#define CATCH_EXCEPTIONS                                                       \
catch (const std::runtime_error& e) {                                          \
  LOG << " runtime_error:  " << std::string(e.what());                         \
  set_last_error("runtime_error:  " + std::string(e.what()));                  \
  return -1;                                                                   \
} catch (const std::invalid_argument& e) {                                     \
  LOG << " invalid_argument:  " << std::string(e.what());                      \
  set_last_error("invalid_argument:  " + std::string(e.what()));               \
  return -1;                                                                   \
} catch (const std::exception& e) {                                            \
  LOG << " exception:  " << std::string(e.what());                             \
  set_last_error("exception:  " + std::string(e.what()));                      \
  return -1;                                                                   \
} catch (...) {                                                                \
  LOG << " unknown critical error";                                            \
  set_last_error("unknown critical error");                                    \
  return -1;                                                                   \
}

// validation logic
template<typename T> void fix_pi_vec(T *pi_vec) { if (*pi_vec < 0) { LOG << " FIX: pi_vec < 0"; *pi_vec = 0; } }
template<typename T> void fix_pi_vec(int numel, T *pi_vec) { for (int i = 0; i < numel; i++) if (pi_vec[i] < 0) { LOG << " FIX: pi_vec < 0"; pi_vec[i] = 0; } }
template<typename T> void fix_num_causal(T *num_causal) { if (*num_causal < 0) { LOG << " FIX: num_causal < 0"; *num_causal = 0; } }
template<typename T> void fix_rho(T *rho) { 
  if (*rho < -1) { LOG << " FIX: rho < -1"; *rho = -1; }; 
  if (*rho > 1) { LOG << " FIX: rho > 1"; *rho = 1; }
}
template<typename T> void fix_rho(int numel, T *rho) {for ( int i = 0; i < numel; i++) fix_rho(&rho[i]); }

void check_trait_index(int trait_index) { if ((trait_index != 1) && (trait_index != 2)) { BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2")); } }
template<typename T> void check_is_positive(T arg) { if (arg <= 0) { LOG << " ERROR: arg <= 0 (value: " << arg << ")"; BGMG_THROW_EXCEPTION(::std::runtime_error("arg <= 0")); } }
template<typename T> void check_is_nonnegative(T arg) { if (arg < 0) { LOG << " ERROR: arg < 0 (value: " << arg << ")"; BGMG_THROW_EXCEPTION(::std::runtime_error("arg < 0")); } }
template<typename T> void check_is_nonnegative(int numel, T* arg) { for (int i = 0; i < numel; i++)  if (arg[i] < 0) { LOG << " ERROR: arg[i] < 0 (i: " << i << ", value: " << arg[i] << ")"; BGMG_THROW_EXCEPTION(::std::runtime_error("arg[i] < 0")); } }
template<typename T> void check_is_not_null(T* ptr) { if (ptr == nullptr) { BGMG_THROW_EXCEPTION(::std::runtime_error("ptr == nullptr")); } }
template<typename T> void check_r2(T arg) { if (arg < 0 | arg > 1) { LOG << " ERROR: arg < 0 | arg > 1 (value: " << arg << ")"; BGMG_THROW_EXCEPTION(::std::runtime_error("arg < 0 | arg > 1")); } }

int64_t bgmg_set_zvec(int context_id, int trait, int length, float* values) {
  try {
    set_last_error(std::string());
    check_trait_index(trait); check_is_positive(length); check_is_not_null(values);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_zvec(trait, length, values);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_nvec(int context_id, int trait, int length, float* values) {
  try {
    set_last_error(std::string());
    check_trait_index(trait); check_is_positive(length); check_is_not_null(values);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_nvec(trait, length, values);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_causalbetavec(int context_id, int trait, int length, float* values) {
  try {
    set_last_error(std::string());
    check_trait_index(trait); check_is_positive(length); check_is_not_null(values);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_causalbetavec(trait, length, values);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_mafvec(int context_id, int length, float* values) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(values);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_mafvec(length, values);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_weights(int context_id, int length, float* values) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(values);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_weights(length, values);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_tag_indices(int context_id, int num_snp, int num_tag, int* tag_indices) {
  try {
    if (!LoggerImpl::singleton().is_initialized()) LoggerImpl::singleton().init("bgmg.log");
    set_last_error(std::string());
    check_is_positive(num_snp); check_is_positive(num_tag); check_is_not_null(tag_indices);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_tag_indices(num_snp, num_tag, tag_indices);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_get_num_tag(int context_id) {
  try {
    set_last_error(std::string());
    return BgmgCalculatorManager::singleton().Get(context_id)->num_tag();
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_get_num_snp(int context_id) {
  try {
    set_last_error(std::string());
    return BgmgCalculatorManager::singleton().Get(context_id)->num_snp();
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_tag_indices(int context_id, int num_tag, int* tag_indices) {
  try {
    set_last_error(std::string());
    check_is_positive(num_tag); check_is_not_null(tag_indices);
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_tag_indices(num_tag, tag_indices);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_ld_r2_coo(int context_id, int chr_label, int64_t length, int* snp_index, int* snp_other_index, float* r) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(snp_index); check_is_not_null(snp_other_index); check_is_not_null(r);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_ld_r2_coo(length, chr_label, snp_index, snp_other_index, r);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_ld_r2_coo_from_file(int context_id, int chr_label, const char* filename) {
  try {
    set_last_error(std::string());
    check_is_not_null(filename);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_ld_r2_coo(chr_label, filename);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_ld_r2_csr(int context_id, int chr_label) {
  try {
    set_last_error(std::string());
    return BgmgCalculatorManager::singleton().Get(context_id)->set_ld_r2_csr(chr_label);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_num_ld_r_tag(int context_id, int tag_index) {
  try {
    set_last_error(std::string());
    return BgmgCalculatorManager::singleton().Get(context_id)->num_ld_r_tag(tag_index);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_ld_r_tag(int context_id, int tag_index, int length, int* snp_index, float* r) {
  try {
    set_last_error(std::string());
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_ld_r_tag(tag_index, length, snp_index, r);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_num_ld_r_chr(int context_id, int chr_label) {
  try {
    set_last_error(std::string());
    return BgmgCalculatorManager::singleton().Get(context_id)->num_ld_r_chr(chr_label);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_ld_r_chr(int context_id, int chr_label, int64_t length, int* tag_index, int* snp_index, float* r) {
  try {
    set_last_error(std::string());
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_ld_r_chr(chr_label, length, tag_index, snp_index, r);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_num_ld_r_tag_range(int context_id, int tag_index_from, int tag_index_to) {
  try {
    set_last_error(std::string());
    return BgmgCalculatorManager::singleton().Get(context_id)->num_ld_r_tag_range(tag_index_from, tag_index_to);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_ld_r_tag_range(int context_id, int tag_index_from, int tag_index_to, int64_t length, int* tag_index, int* snp_index, float* r) {
  try {
    set_last_error(std::string());
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_ld_r_tag_range(tag_index_from, tag_index_to, length, tag_index, snp_index, r);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_ld_sum_r2(int context_id, int length, float* buffer) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(buffer);
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_ld_sum_r2(length, buffer);
  } CATCH_EXCEPTIONS;
}
int64_t bgmg_retrieve_ld_sum_r4(int context_id, int length, float* buffer) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(buffer);
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_ld_sum_r4(length, buffer);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_zvec(int context_id, int trait, int length, float* buffer) {
  try {
    set_last_error(std::string());
    check_trait_index(trait); check_is_positive(length); check_is_not_null(buffer);
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_zvec(trait, length, buffer);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_nvec(int context_id, int trait, int length, float* buffer) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(buffer);
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_nvec(trait, length, buffer);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_causalbetavec(int context_id, int trait, int length, float* buffer) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(buffer);
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_causalbetavec(trait, length, buffer);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_fixed_effect_delta(int context_id, int trait, int length, float* buffer) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(buffer);
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_fixed_effect_delta(trait, length, buffer);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_mafvec(int context_id, int length, float* buffer) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(buffer);
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_mafvec(length, buffer);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_weights(int context_id, int length, float* buffer) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(buffer);
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_weights(length, buffer);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_perform_ld_clump(int context_id, float r2, int length, float* buffer) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(buffer);
    return BgmgCalculatorManager::singleton().Get(context_id)->perform_ld_clump(r2, length, buffer);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_option(int context_id, char* option, double value) {
  try {
    set_last_error(std::string());
    check_is_not_null(option);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_option(option, value);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_weights_hardprune(int context_id, float hardprune_subset, float hardprune_r2, float hardprune_maf, int use_w_ld) {
  try {
    set_last_error(std::string());
    check_is_nonnegative(hardprune_subset);
    check_r2(hardprune_r2);
    check_is_nonnegative(hardprune_maf);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_weights_hardprune(hardprune_subset, hardprune_r2, hardprune_maf, use_w_ld);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_weights_randprune(int context_id, int n, float r2, float maf, int use_w_ld, const char* exclude, const char* extract) {
  try {
    set_last_error(std::string());
    check_is_positive(n); check_r2(r2); check_is_nonnegative(maf);
    check_is_not_null(exclude);
    check_is_not_null(extract);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_weights_randprune(n, r2, maf, use_w_ld, exclude, extract);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_dispose(int context_id) {
  try {
    set_last_error(std::string());
    if (context_id == -1) BgmgCalculatorManager::singleton().Clear();
    else BgmgCalculatorManager::singleton().Erase(context_id);
    return 0;
  } CATCH_EXCEPTIONS;
}

const char* bgmg_status(int context_id) {
  return "";
}

void bgmg_init_log(const char* file) {
  check_is_not_null(file);
  LoggerImpl::singleton().init(file);
  Logger::singleton() << "=bgmg plugin version: " << VERSION;
}

void bgmg_log_message(const char* message) {
  if (!LoggerImpl::singleton().is_initialized()) {
    LoggerImpl::singleton().init("bgmg.log");
  }

  std::vector<std::string> tokens = Logger::tokenize_message(message);
  for (auto token: tokens)
    Logger::singleton() << "=" << token;
}

int64_t bgmg_set_chrnumvec(int context_id, int length, int* values) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(values);
    return BgmgCalculatorManager::singleton().Get(context_id)->set_chrnumvec(length, values);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_retrieve_chrnumvec(int context_id, int length, int* buffer) {
  try {
    set_last_error(std::string());
    check_is_positive(length); check_is_not_null(buffer);
    return BgmgCalculatorManager::singleton().Get(context_id)->retrieve_chrnumvec(length, buffer);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_init(int context_id, const char* bim_file, const char* frq_file, const char* chr_labels, const char* trait1_file, const char* trait2_file, const char* exclude, const char* extract, const char* exclude_ranges) {
  try {
    set_last_error(std::string());
    check_is_not_null(bim_file);
    check_is_not_null(frq_file);
    check_is_not_null(chr_labels);
    check_is_not_null(trait1_file);
    check_is_not_null(trait2_file);
    check_is_not_null(exclude);
    check_is_not_null(extract);
    check_is_not_null(exclude_ranges);    
    return BgmgCalculatorManager::singleton().Get(context_id)->init(bim_file, frq_file, chr_labels, trait1_file, trait2_file, exclude, extract, exclude_ranges);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_read_trait_file(int context_id, int trait_index, const char* trait_file, const char* exclude, const char* extract, const char* exclude_ranges) {
  try {
    set_last_error(std::string());
    check_is_not_null(trait_file);
    check_is_not_null(exclude);
    check_is_not_null(extract);
    check_is_not_null(exclude_ranges);    
    return BgmgCalculatorManager::singleton().Get(context_id)->read_trait_file(trait_index, trait_file, exclude, extract, exclude_ranges);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_convert_plink_ld(int context_id, const char* plink_ld_gz, const char* plink_ld_bin) {
  try {
    set_last_error(std::string());
    check_is_not_null(plink_ld_gz);
    check_is_not_null(plink_ld_bin);
    return BgmgCalculatorManager::singleton().Get(context_id)->convert_plink_ld(plink_ld_gz, plink_ld_bin);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_calc_ld_matrix(const char* bfile, const char* outfile, double r2min, double ldscore_r2min, int ld_window, float ld_window_kb) {
  try {
    if (!LoggerImpl::singleton().is_initialized()) LoggerImpl::singleton().init("bgmg.log");
    set_last_error(std::string());
    check_is_not_null(bfile); check_is_not_null(outfile);
    generate_ld_matrix_from_bed_file(bfile, r2min, ldscore_r2min, ld_window, ld_window_kb, outfile);
    return 0;
  } CATCH_EXCEPTIONS;
}

void check_and_fix_unified_univariate(int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL) {
  check_is_positive(num_components); check_is_positive(num_snp); 
  fix_pi_vec(num_snp*num_components, pi_vec); check_is_nonnegative(num_snp*num_components, sig2_vec);
  check_is_nonnegative(sig2_zeroA); check_is_nonnegative(sig2_zeroC); check_is_nonnegative(sig2_zeroL); 
}

void check_and_fix_unified_bivariate(int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, 
                                     float* sig2_zeroL, float *rho_zeroA, float *rho_zeroL) {
  check_is_positive(num_snp);
  fix_pi_vec(num_snp*3, pi_vec);
  check_is_nonnegative(num_snp*2, sig2_vec);
  fix_rho(num_snp, rho_vec);
  check_is_nonnegative(2, sig2_zeroA); check_is_nonnegative(2, sig2_zeroC); check_is_nonnegative(2, sig2_zeroL); 
  fix_rho(rho_zeroA);
  fix_rho(rho_zeroL);
}

double bgmg_calc_unified_univariate_cost(int context_id, int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float* aux) {
  try {
    set_last_error(std::string());
    check_trait_index(trait_index); check_and_fix_unified_univariate(num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL);
    return BgmgCalculatorManager::singleton().Get(context_id)->calc_unified_univariate_cost(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, aux);
  } CATCH_EXCEPTIONS;
}
  
int64_t bgmg_calc_unified_univariate_pdf(int context_id, int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, int length, float* zvec, float* pdf) {
  try {
    set_last_error(std::string());
    check_trait_index(trait_index); check_and_fix_unified_univariate(num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL); check_is_positive(length); check_is_not_null(zvec); check_is_not_null(pdf);
    return BgmgCalculatorManager::singleton().Get(context_id)->calc_unified_univariate_pdf(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, length, zvec, pdf);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_calc_unified_univariate_power(int context_id, int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, float zthresh, int length, float* nvec, float* svec_num, float* svec_denom) {
  try {
    set_last_error(std::string());
    check_trait_index(trait_index); check_and_fix_unified_univariate(num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL); check_is_positive(length); check_is_not_null(nvec); check_is_not_null(svec_num); check_is_not_null(svec_denom);
    return BgmgCalculatorManager::singleton().Get(context_id)->calc_unified_univariate_power(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, zthresh, length, nvec, svec_num, svec_denom);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_calc_unified_univariate_delta_posterior(int context_id, int trait_index, int num_components, int num_snp, float* pi_vec, float* sig2_vec, float sig2_zeroA, float sig2_zeroC, float sig2_zeroL, int length, float* c0, float* c1, float* c2) {
  try {
    set_last_error(std::string());
    check_trait_index(trait_index); check_and_fix_unified_univariate(num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL); check_is_positive(length); check_is_not_null(c0); check_is_not_null(c1); check_is_not_null(c2);
    return BgmgCalculatorManager::singleton().Get(context_id)->calc_unified_univariate_delta_posterior(trait_index, num_components, num_snp, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, length, c0, c1, c2);
  } CATCH_EXCEPTIONS;
}

double bgmg_calc_unified_bivariate_cost(int context_id, int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, float* aux) {
  try {
    set_last_error(std::string());
    check_and_fix_unified_bivariate(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC,  sig2_zeroL, &rho_zeroA, &rho_zeroL);
    return BgmgCalculatorManager::singleton().Get(context_id)->calc_unified_bivariate_cost(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC,  sig2_zeroL, rho_zeroA, rho_zeroL, aux);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_calc_unified_bivariate_pdf(int context_id, int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL, int length, float* zvec1, float* zvec2, float* pdf) {
  try {
    set_last_error(std::string());
    check_and_fix_unified_bivariate(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC,  sig2_zeroL, &rho_zeroA, &rho_zeroL);
    check_is_positive(length); check_is_not_null(zvec1); check_is_not_null(zvec2); check_is_not_null(pdf);
    return BgmgCalculatorManager::singleton().Get(context_id)->calc_unified_bivariate_pdf(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC,  sig2_zeroL, rho_zeroA, rho_zeroL, length, zvec1, zvec2, pdf);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_calc_unified_bivariate_delta_posterior(int context_id, int num_snp, float* pi_vec, float* sig2_vec, float* rho_vec, float* sig2_zeroA, float* sig2_zeroC, float* sig2_zeroL, float rho_zeroA, float rho_zeroL,
                                                    int length, float* c00, float* c10, float* c01, float* c20, float* c11, float* c02) {
  try {
    set_last_error(std::string());
    check_and_fix_unified_bivariate(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC,  sig2_zeroL, &rho_zeroA, &rho_zeroL);
    check_is_positive(length); check_is_not_null(c00); check_is_not_null(c10); check_is_not_null(c01); check_is_not_null(c20); check_is_not_null(c11); check_is_not_null(c02);
    return BgmgCalculatorManager::singleton().Get(context_id)->calc_unified_bivariate_delta_posterior(num_snp, pi_vec, sig2_vec, rho_vec, sig2_zeroA, sig2_zeroC,  sig2_zeroL, rho_zeroA, rho_zeroL, length, c00, c10, c01, c20, c11, c02);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_save_context(int context_id, const char* file) {
  try {
    set_last_error(std::string());
    return BgmgCalculatorManager::singleton().Get(context_id)->save(file);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_load_context(int context_id, const char* file) {
  try {
    set_last_error(std::string());
    return BgmgCalculatorManager::singleton().Get(context_id)->load(file);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_copy_request_impl(int64_t length, char* address, std::string* source) {
  try {
    if (source == nullptr) {
      std::stringstream ss;
      ss << "There is no data to copy; check if bgmg_request_xxx method is executed before copying the result";
      set_last_error(ss.str().c_str());
      return -1;
    }

    if (length != static_cast<int64_t>(source->size())) {
      std::stringstream ss;
      ss << "Invalid 'length' parameter ";
      ss << "(" << source->size() << " expected, found " << length << ").";
      set_last_error(ss.str());
      return -1;
    }

    memcpy(address, StringAsArray(source), length);

    return 0;
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_copy_requested_object(int64_t length, char* address) {
  LOG << " bgmg_copy_requested_object is copying " << length << " bytes...";
  return bgmg_copy_request_impl(length, address, last_buffer());
}

int64_t bgmg_request_annotate_intervals(int num_snp, int64_t* posvec, int num_interval, int64_t *interval_start, int64_t* interval_end, int* interval_label) {
  try {
    set_last_error(std::string());
    return annotate_intervals(num_snp, posvec, num_interval, interval_start, interval_end, interval_label, last_buffer());
  } CATCH_EXCEPTIONS;
}