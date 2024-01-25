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

#include "bgmg_calculator_impl.h"

#include <immintrin.h>  // _mm_setcsr, _mm_getcsr

std::vector<float>* BgmgCalculator::get_zvec(int trait_index) {
  if ((trait_index != 1) && (trait_index != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  return (trait_index == 1) ? &zvec1_ : &zvec2_;
}

std::vector<float>* BgmgCalculator::get_nvec(int trait_index) {
  if ((trait_index != 1) && (trait_index != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  return (trait_index == 1) ? &nvec1_ : &nvec2_;
}

std::vector<float>* BgmgCalculator::get_causalbetavec(int trait_index) {
  if ((trait_index != 1) && (trait_index != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  return (trait_index == 1) ? &causalbetavec1_ : &causalbetavec2_;
}

BgmgCalculator::BgmgCalculator() : num_snp_(-1), num_tag_(-1), k_max_(100), k_max_pdf_(100), seed_(0), aux_option_(AuxOption_Ezvec2),
    use_complete_tag_indices_(false), disable_snp_to_tag_map_(false), allow_ambiguous_snps_(false), r2_min_(0.0), z1max_(DEFAULT_ZMAX), z2max_(DEFAULT_ZMAX), ld_format_version_(-1), retrieve_ld_sum_type_(0), num_components_(1), 
    max_causals_(100000), cost_calculator_(CostCalculator_Sampling), cache_tag_r2sum_(false), ld_matrix_csr_(*this),
    cubature_abs_error_(0), cubature_rel_error_(1e-4), cubature_max_evals_(0), calc_k_pdf_(false) {
  boost::posix_time::ptime const time_epoch(boost::gregorian::date(1970, 1, 1));
  seed_ = (boost::posix_time::microsec_clock::local_time() - time_epoch).ticks();

  // flush denormals to zero --- implemented only for GCC. Need the same for clang and MS VS.
  // https://stackoverflow.com/questions/9314534/why-does-changing-0-1f-to-0-slow-down-performance-by-10x
  // https://carlh.net/plugins/denormals.php
  // #if defined(__clang__)
  // #include <fenv.h>
  // fesetenv(FE_DFL_DISABLE_SSE_DENORMS_ENV);
  // #elif  defined(_MSC_VER)
  // #include <immintrin.h>
  // _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  // _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

  _mm_setcsr( _mm_getcsr() | (1<<15) | (1<<6));
}

void BgmgCalculator::check_num_snp(int length) {
  if (num_snp_ == -1) BGMG_THROW_EXCEPTION(::std::runtime_error("call set_tag_indices first"));
  if (num_snp_ != length) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_snps_"));
}

void BgmgCalculator::check_num_tag(int length) {
  if (num_tag_ == -1) BGMG_THROW_EXCEPTION(::std::runtime_error("call set_tag_indices first"));
  if (num_tag_ != length) BGMG_THROW_EXCEPTION(::std::runtime_error("length != num_snps_"));
}

int64_t BgmgCalculator::set_zvec(int trait, int length, float* values) {
  if ((trait != 1) && (trait != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));

  int num_undef = 0;
  for (int i = 0; i < length; i++) if (!std::isfinite(values[i])) num_undef++;
  LOG << " set_zvec(trait=" << trait << "); num_undef=" << num_undef;
  check_num_tag(length);
  get_zvec(trait)->assign(values, values + length);
  return 0;
}

int64_t BgmgCalculator::set_nvec(int trait, int length, float* values) {
  if ((trait != 1) && (trait != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  
  int num_undef = 0;
  for (int i = 0; i < length; i++) if (!std::isfinite(values[i])) num_undef++;
  LOG << " set_nvec(trait=" << trait << "); num_undef=" << num_undef;
  check_num_tag(length);
  get_nvec(trait)->assign(values, values + length);
  return 0;
}

int64_t BgmgCalculator::set_causalbetavec(int trait, int length, float* values) {
  if ((trait != 1) && (trait != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  
  int num_undef = 0;
  for (int i = 0; i < length; i++) if (!std::isfinite(values[i])) num_undef++;
  if (num_undef > 0) BGMG_THROW_EXCEPTION(::std::runtime_error("undefined values not allowed in causalbetavec, use zero instead"));
  LOG << " set_causalbetavec(trait=" << trait << "); num_undef=" << num_undef;
  check_num_snp(length);
  get_causalbetavec(trait)->assign(values, values + length);
  return 0;
}

int64_t BgmgCalculator::set_weights(int length, float* values) {
  check_num_tag(length);

  int nnz = 0;
  float sum = 0.0f;
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
      if (values[i] != 0) {
        if ((!zvec1_.empty()) && !std::isfinite(zvec1_[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("can not set a non-zero weight for missing value in zvec1"));
        if ((!nvec1_.empty()) && !std::isfinite(nvec1_[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("can not set a non-zero weight for missing value in nvec1"));
        if ((!zvec2_.empty()) && !std::isfinite(zvec2_[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("can not set a non-zero weight for missing value in zvec2"));
        if ((!nvec2_.empty()) && !std::isfinite(nvec2_[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("can not set a non-zero weight for missing value in nvec2"));
      nnz++;
      sum += values[i];
    }
  }

  LOG << " set_weights(length=" << length << "), nnz=" << nnz << ", sum=" << sum << "; ";
  weights_.assign(values, values+length);
  return 0;
}

int64_t BgmgCalculator::set_option(char* option, double value) {
  LOG << " set_option(" << option << "=" << value << "); ";

  if (!strcmp(option, "diag")) {
    log_diagnostics(); return 0;
  } else if (!strcmp(option, "kmax")) {
    clear_state(); k_max_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "kmax_pdf")) {
    clear_state(); k_max_pdf_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "r2min")) {
    r2_min_ = value; return 0;
  } else if (!strcmp(option, "max_causals")) {
    if (!snp_order_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't change max_causals after find_snp_order"));
    clear_state(); max_causals_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "num_components")) {
    if (!snp_order_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can't change num_components after find_snp_order"));
    clear_state(); num_components_ = static_cast<int>(value); return 0;
  } else if (!strcmp(option, "seed")) {
    seed_ = static_cast<int64_t>(value); return 0;
  } else if (!strcmp(option, "calc_k_pdf")) {
    calc_k_pdf_ = (value != 0); return 0;
  } else if (!strcmp(option, "cubature_max_evals")) {
    cubature_max_evals_ = static_cast<int64_t>(value); return 0;
  } else if (!strcmp(option, "cubature_abs_error")) {
    cubature_abs_error_ = value; return 0;
  } else if (!strcmp(option, "cubature_rel_error")) {
    cubature_rel_error_ = value; return 0;
  } else if (!strcmp(option, "fast_cost")) {
    cost_calculator_ = (value != 0) ? CostCalculator_Sampling : CostCalculator_Gaussian; return 0;
  } else if (!strcmp(option, "cost_calculator")) {
    int int_value = (int)value;
    if (int_value < 0 || int_value >= CostCalculator_MAX) BGMG_THROW_EXCEPTION(::std::runtime_error("cost_calculator value must be 0 (Sampling), 1 (Gaussian), 2 (Convolve) or 3 (Smplfast)"));
    cost_calculator_ = (CostCalculator)int_value; return 0;
  } else if (!strcmp(option, "aux_option")) {
    int int_value = (int)value;
    if (int_value < 0 || int_value >= AuxOption_MAX) BGMG_THROW_EXCEPTION(::std::runtime_error("aux_option value must be 0 (None), 1 (Ezvec2), 2 (TagPdf), or 3 (TagPdfErr)"));
    aux_option_ = (AuxOption)int_value; return 0;
  } else if (!strcmp(option, "z1max")) {
    if (value <= 0) BGMG_THROW_EXCEPTION(::std::runtime_error("zmax must be positive"));
    z1max_ = value; return 0;
  } else if (!strcmp(option, "z2max")) {
    if (value <= 0) BGMG_THROW_EXCEPTION(::std::runtime_error("zmax must be positive"));
    z2max_ = value; return 0;
  } else if (!strcmp(option, "ld_format_version")) {
    ld_format_version_ = int(value); return 0;
  } else if (!strcmp(option, "retrieve_ld_sum_type")) {
    retrieve_ld_sum_type_ = int(value); return 0;
  } else if (!strcmp(option, "use_complete_tag_indices")) {
    use_complete_tag_indices_ = (value != 0); return 0;
  } else if (!strcmp(option, "disable_snp_to_tag_map")) {
    disable_snp_to_tag_map_ = (value != 0); return 0;
  } else if (!strcmp(option, "allow_ambiguous_snps")) {
    allow_ambiguous_snps_ = (value != 0); return 0;
  } else if (!strcmp(option, "threads")) {
    if (value > 0) {
      LOG << " omp_set_num_threads(" << static_cast<int>(value) << ")";
      omp_set_num_threads(static_cast<int>(value));
    }
    return 0;
  } else if (!strcmp(option, "cache_tag_r2sum")) {
    cache_tag_r2sum_ = (value != 0);
    for (int component_id = 0; component_id < num_components_; component_id++) clear_tag_r2sum(component_id);
    return 0;
  }

  BGMG_THROW_EXCEPTION(::std::runtime_error("unknown option"));
  return 0;
}

int64_t BgmgCalculator::set_tag_indices(int num_snp, int num_tag, int* tag_indices) {
  if (num_snp_ != -1 || num_tag_ != -1) BGMG_THROW_EXCEPTION(::std::runtime_error("can not call set_tag_indices twice"));

  LOG << " set_tag_indices(num_snp=" << num_snp << ", num_tag=" << num_tag << "); ";
  num_snp_ = num_snp;
  num_tag_ = num_tag;

  is_tag_.resize(num_snp, 0);
  snp_to_tag_.resize(num_snp, -1);
  tag_to_snp_.assign(tag_indices, tag_indices + num_tag);
  for (int i = 0; i < tag_to_snp_.size(); i++) {
    CHECK_SNP_INDEX((*this), tag_to_snp_[i]);
    is_tag_[tag_to_snp_[i]] = 1;
    snp_to_tag_[tag_to_snp_[i]] = i;
  }

  return 0;
}

int64_t BgmgCalculator::retrieve_tag_indices(int num_tag, int* tag_indices) {
  if (num_tag != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (num_tag != tag_to_snp_.size()) BGMG_THROW_EXCEPTION(::std::runtime_error("num_tag != tag_to_snp_.size()"));
  LOG << " retrieve_tag_indices()";
  for (int i = 0; i < num_tag_; i++) tag_indices[i] = tag_to_snp_[i];
  return 0;
}

int64_t BgmgCalculator::set_ld_r2_coo(int chr_label, int64_t length, int* snp_index, int* tag_index, float* r) {
  return ld_matrix_csr_.set_ld_r2_coo(chr_label, length, snp_index, tag_index, r, r2_min_);
}

int64_t BgmgCalculator::set_ld_r2_coo(int chr_label, const std::string& filename) {
  if (ld_format_version_ == 0)
    return ld_matrix_csr_.set_ld_r2_coo_version0(chr_label, filename, r2_min_);

  if (mafvec_.empty()) {
    LOG << " initialize mafvec";
    mafvec_.assign(num_snp_, NAN);
  }

  return ld_matrix_csr_.set_ld_r2_coo_version1plus(chr_label, filename, r2_min_);
}

int64_t BgmgCalculator::set_ld_r2_csr(int chr_label) {
  int64_t retval = ld_matrix_csr_.set_ld_r2_csr(r2_min_, chr_label);
  return retval;
}

int64_t BgmgCalculator::set_mafvec(int length, float* values) {
  for (int i = 0; i < length; i++) {
    if (!std::isfinite(values[i])) BGMG_THROW_EXCEPTION(::std::runtime_error("encounter undefined values"));
  }

  if (!mafvec_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can not set mafvec twice"));

  LOG << ">set_mafvec(" << length << "); ";
  check_num_snp(length);
  mafvec_.assign(values, values + length);
  LOG << "<set_mafvec(" << length << "); ";
  return 0;
}

int64_t BgmgCalculator::retrieve_ld_sum_r2(int length, float* buffer) {
  check_num_snp(length);
  LOG << " retrieve_ld_sum_r2()";
  for (int snp_index = 0; snp_index < num_snp_; snp_index++) {
    if (retrieve_ld_sum_type_ == 0) buffer[snp_index] = ld_matrix_csr_.ld_sum()->ld_sum_r2_above_r2min()[snp_index];
    else if (retrieve_ld_sum_type_ == 1) buffer[snp_index] = ld_matrix_csr_.ld_sum()->ld_sum_r2_below_r2min()[snp_index];
    else if (retrieve_ld_sum_type_ == 2) buffer[snp_index] = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_sum_r2_above_r2min()[snp_index];
    else if (retrieve_ld_sum_type_ == 3) buffer[snp_index] = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_sum_r2_below_r2min()[snp_index];
    else BGMG_THROW_EXCEPTION(::std::runtime_error("invalid retrieve_ld_sum_type_"));
  }
  return 0;
}

int64_t BgmgCalculator::retrieve_ld_sum_r4(int length, float* buffer) {
  check_num_snp(length);
  LOG << " retrieve_ld_sum_r4()";
  for (int snp_index = 0; snp_index < num_snp_; snp_index++) {
    if (retrieve_ld_sum_type_ == 0 || retrieve_ld_sum_type_ == 1) buffer[snp_index] = ld_matrix_csr_.ld_sum()->ld_sum_r4_above_r2min()[snp_index];
    else if (retrieve_ld_sum_type_ == 2 || retrieve_ld_sum_type_ == 3) buffer[snp_index] = ld_matrix_csr_.ld_sum_adjust_for_hvec()->ld_sum_r4_above_r2min()[snp_index];
    else BGMG_THROW_EXCEPTION(::std::runtime_error("invalid retrieve_ld_sum_type_"));
  }
  return 0;
}

template<typename T>
std::string std_vector_to_str(const std::vector<T>& vec) {
  std::stringstream ss;
  int max_el = std::min<int>(5, vec.size() - 1);
  ss << "[";
  for (int i = 0; i < max_el; i++) {
    bool last = (i == (max_el - 1));
    ss << vec[i];
    if (last) ss << ", ...";
    else ss << ", ";
  }
  ss << "]";

  size_t nnz = 0;
  T sum = (T)0;
  for (size_t i = 0; i < vec.size(); i++) if (std::isfinite(vec[i]) && (vec[i] != 0)) { nnz++; sum += vec[i]; }
  ss << ", nnz=" << nnz << ", sum=" << sum;
  return ss.str();
}

void BgmgCalculator::log_diagnostics() {
  size_t mem_bytes = 0, mem_bytes_total = 0;
  LOG << " diag: num_snp_=" << num_snp_;
  LOG << " diag: num_tag_=" << num_tag_;
  mem_bytes_total += ld_matrix_csr_.log_diagnostics();
  LOG << " diag: zvec1_.size()=" << zvec1_.size();
  LOG << " diag: zvec1_=" << std_vector_to_str(zvec1_);
  LOG << " diag: nvec1_.size()=" << nvec1_.size();
  LOG << " diag: nvec1_=" << std_vector_to_str(nvec1_);
  LOG << " diag: causalbetavec1_.size()=" << causalbetavec1_.size();
  LOG << " diag: causalbetavec1_=" << std_vector_to_str(causalbetavec1_);
  LOG << " diag: zvec2_.size()=" << zvec2_.size();
  LOG << " diag: zvec2_=" << std_vector_to_str(zvec2_);
  LOG << " diag: nvec2_.size()=" << nvec2_.size();
  LOG << " diag: nvec2_=" << std_vector_to_str(nvec2_);
  LOG << " diag: causalbetavec2_.size()=" << causalbetavec2_.size();
  LOG << " diag: causalbetavec2_=" << std_vector_to_str(causalbetavec2_);
  LOG << " diag: weights_.size()=" << weights_.size();
  LOG << " diag: weights_=" << std_vector_to_str(weights_);
  LOG << " diag: mafvec_.size()=" << mafvec_.size();
  LOG << " diag: mafvec_=" << std_vector_to_str(mafvec_);
  for (int i = 0; i < snp_order_.size(); i++) {
    mem_bytes = snp_order_[i]->size() * sizeof(int); mem_bytes_total += mem_bytes;
    LOG << " diag: snp_order_[" << i << "].shape=[" << snp_order_[i]->no_rows() << ", " << snp_order_[i]->no_columns() << "]" << " (mem usage = " << mem_bytes << " bytes)";
    LOG << " diag: snp_order_[" << i << "]=" << snp_order_[i]->to_str();
  }
  for (int i = 0; i < tag_r2sum_.size(); i++) {
    mem_bytes = tag_r2sum_[i]->size() * sizeof(float); mem_bytes_total += mem_bytes;
    LOG << " diag: tag_r2sum_[" << i << "].shape=[" << tag_r2sum_[i]->no_rows() << ", " << tag_r2sum_[i]->no_columns() << "]" << " (mem usage = " << mem_bytes << " bytes)";
    LOG << " diag: tag_r2sum_[" << i << "]=" << tag_r2sum_[i]->to_str();
  }
  for (int i = 0; i < last_num_causals_.size(); i++) 
    LOG << " diag: last_num_causals_[" << i << "]=" << last_num_causals_[i];
  LOG << " diag: options.k_max_=" << k_max_;
  LOG << " diag: options.k_max_pdf_=" << k_max_pdf_;  
  LOG << " diag: options.use_complete_tag_indices_=" << use_complete_tag_indices_;
  LOG << " diag: options.disable_snp_to_tag_map_=" << disable_snp_to_tag_map_;
  LOG << " diag: options.allow_ambiguous_snps_=" << allow_ambiguous_snps_;
  LOG << " diag: options.max_causals_=" << max_causals_;
  LOG << " diag: options.num_components_=" << num_components_;
  LOG << " diag: options.r2_min_=" << r2_min_;
  LOG << " diag: options.z1max_=" << z1max_;
  LOG << " diag: options.z2max_=" << z2max_;
  LOG << " diag: options.cost_calculator_=" << ((int)cost_calculator_) <<
    ((cost_calculator_==CostCalculator_Sampling) ? " (Sampling)" :
     (cost_calculator_==CostCalculator_Gaussian) ? " (Gaussian)" :
     (cost_calculator_==CostCalculator_Convolve) ? " (Convolve)" : " (Unknown)");
  LOG << " diag: options.aux_option_=" << ((int)aux_option_) <<
    ((aux_option_==AuxOption_None) ? " (None)" :
     (aux_option_==AuxOption_Ezvec2) ? " (Ezvec2)" :
     (aux_option_==AuxOption_TagPdf) ? " (TagPdf)" : 
     (aux_option_==AuxOption_TagPdfErr) ? " (TagPdfErr)" : " (Unknown)");
  LOG << " diag: options.cache_tag_r2sum_=" << (cache_tag_r2sum_ ? "yes" : "no");
  LOG << " diag: options.seed_=" << (seed_);
  LOG << " diag: options.cubature_abs_error_=" << (cubature_abs_error_);
  LOG << " diag: options.cubature_rel_error_=" << (cubature_rel_error_);
  LOG << " diag: options.cubature_max_evals_=" << (cubature_max_evals_);
  LOG << " diag: options.calc_k_pdf_=" << (calc_k_pdf_);
  LOG << " diag: options.ld_format_version_=" << (ld_format_version_);
  LOG << " diag: options.retrieve_ld_sum_type_=" << (retrieve_ld_sum_type_);
  LOG << " diag: Estimated memory usage (total): " << mem_bytes_total << " bytes";
}

void apply_extract(std::string extract, const BimFile& bim_file, std::vector<int> *defvec) {
  SnpList extract_object;
  extract_object.read(extract);
  for (int i = 0; i < bim_file.size(); i++) {
    if (!extract_object.contains(bim_file.snp()[i]))
      defvec->at(i) = 0;
  }

  LOG << " constrain analysis to " << std::accumulate(defvec->begin(), defvec->end(), 0) << " tag variants (due to extract='" << extract << "')";
}

void apply_exclude(std::string exclude, const BimFile& bim_file, std::vector<int>* defvec) {
  SnpList exclude_object;
  exclude_object.read(exclude);
  for (int i = 0; i < bim_file.size(); i++) {
    if (exclude_object.contains(bim_file.snp()[i]))
      defvec->at(i) = 0;
  }

  LOG << " constrain analysis to " << std::accumulate(defvec->begin(), defvec->end(), 0) << " tag variants (due to exclude='" << exclude << "')";
}

void apply_exclude_ranges(std::string exclude_ranges, const BimFile& bim_file, std::vector<int>* defvec) {
  const std::string separators = ", \t\n\r";
  boost::trim_if(exclude_ranges, boost::is_any_of(separators));
  boost::to_lower(exclude_ranges);
  std::vector<std::string> str_vec;
  boost::split(str_vec, exclude_ranges, boost::is_any_of(separators), boost::token_compress_on);

  for (auto str_range : str_vec) {
    int removed_count = 0;

    std::vector<std::string> tokens;
    boost::split(tokens, str_range, boost::is_any_of(":-"), boost::token_compress_on);
    int chrnum, from_bp, to_bp;
    try {
      chrnum = stoi(tokens[0]);
      from_bp = stoi(tokens[1]);
      to_bp = stoi(tokens[2]);
    }
    catch (...) {
      std::stringstream error_str;
      error_str << "Error parsing: " << str_range;
      throw std::invalid_argument(error_str.str());
    }

    ChrBpRange chr_bp_range(chrnum, from_bp, to_bp);
    for (int i = 0; i < bim_file.size(); i++) {
      if (chr_bp_range.contains(bim_file.chr_label()[i], bim_file.bp()[i])) {
        defvec->at(i) = 0;
        removed_count++;
      }
    }

    LOG << " constrain analysis to " << std::accumulate(defvec->begin(), defvec->end(), 0) << " tag variants (" << removed_count << " removed due to exclude_ranges='" << str_range << "')";
  }
}

void apply_maf_threshold(float maf_threshold, const std::vector<float>& mafvec, std::vector<int> *defvec) {
  for (int i = 0; i < mafvec.size(); i++) {
    if (mafvec[i] < maf_threshold || mafvec[i] > (1.0f-maf_threshold))
      defvec->at(i) = 0;
  }

  LOG << " constrain analysis to " << std::accumulate(defvec->begin(), defvec->end(), 0) << " tag variants (due to maf='" << maf_threshold << "')";
}

int64_t BgmgCalculator::perform_ld_clump(float r2_threshold, int length, float* buffer) {
  check_num_tag(length);
  if (r2_threshold < r2_min_) BGMG_THROW_EXCEPTION(::std::runtime_error("perform_ld_clump: r2 < r2_min_"));
  LOG << ">perform_ld_clump(length=" << length << ", r2=" << r2_threshold << ")";
  SimpleTimer timer(-1);

  LdMatrixRow ld_matrix_row;

  std::vector<std::tuple<float, int>> sorted_buffer;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (std::isfinite(buffer[tag_index])) {
      sorted_buffer.push_back(std::tuple<float, int>(buffer[tag_index], tag_index));
    }
  }
  std::sort(sorted_buffer.rbegin(), sorted_buffer.rend());  // use reverse iterators to sort in descending order

  std::vector<char> processed_tag_indices(num_tag_, 0);
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (!std::isfinite(buffer[tag_index]))
      processed_tag_indices[tag_index] = 1;
  }

  int count = 0;
  for (int sorted_index = 0; sorted_index < sorted_buffer.size(); sorted_index++) {
    const int tag_index = std::get<1>(sorted_buffer[sorted_index]);
    if (processed_tag_indices[tag_index]) continue;
    
    count++;
    processed_tag_indices[tag_index] = 1;

    // Step 3. eliminate all unprocessed tag SNPs in LD with max_index
    ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), &ld_matrix_row);
    auto iter_end = ld_matrix_row.end();
    for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
      const int snp_index_in_ld = iter.index();
      const int tag_index_in_ld = snp_to_tag_[snp_index_in_ld];
      if (tag_index_in_ld==-1) continue;
      const float r2_value = iter.r2();  // here we are interested in r2 (hvec is irrelevant)        
      if (r2_value <= r2_threshold) continue;
      if (processed_tag_indices[tag_index_in_ld]) continue;
      processed_tag_indices[tag_index_in_ld] = 1;
      buffer[tag_index_in_ld] = NAN;
    }   
  }

  LOG << "<perform_ld_clump(length=" << length << ", r2=" << r2_threshold << "), elapsed time " << timer.elapsed_ms() << "ms, keep " << count << " tag SNPs";
  return 0;
}

int64_t BgmgCalculator::set_weights_randprune(int n, float r2_threshold, float maf_threshold, bool use_w_ld) {
  return set_weights_randprune(n, r2_threshold, maf_threshold, use_w_ld, std::string(), std::string());
}

int64_t BgmgCalculator::set_weights_randprune(int n, float r2_threshold, float maf_threshold, bool use_w_ld, std::string exclude, std::string extract) {
  std::stringstream ss;
  ss << "set_weights_randprune(n=" << n << ", r2=" << r2_threshold << ", maf=" << maf_threshold << ", use_w_ld=" << (use_w_ld ? "yes" : "no") << ", exclude=" << exclude << ", extract=" << extract << ")";
  LOG << ">" << ss.str();

  if (r2_threshold < r2_min_) BGMG_THROW_EXCEPTION(::std::runtime_error("set_weights_randprune: r2 < r2_min_"));
  if (n <= 0) BGMG_THROW_EXCEPTION(::std::runtime_error("set_weights_randprune: n <= 0"));
  if (!weights_.empty() && (weights_.size() != num_tag_)) BGMG_THROW_EXCEPTION(::std::runtime_error("weights_ size does not match num_tag_"));
  SimpleTimer timer(-1);

  std::valarray<float> weight_random_pruning(0.0f, num_tag_);

  // random pruning should operate only on SNPs with defined zvec and nvec.
  // if weights are already defined, it should retain zero elements
  std::vector<int> defvec(num_snp_, 0);
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if ((!zvec1_.empty()) && !std::isfinite(zvec1_[tag_index])) continue;
    if ((!nvec1_.empty()) && !std::isfinite(nvec1_[tag_index])) continue;
    if ((!zvec2_.empty()) && !std::isfinite(zvec2_[tag_index])) continue;
    if ((!nvec2_.empty()) && !std::isfinite(nvec2_[tag_index])) continue;
    if ((!weights_.empty()) && (weights_[tag_index] == 0)) continue;
    defvec[tag_to_snp()[tag_index]] = 1;
  }
  LOG << " constrain analysis to " << std::accumulate(defvec.begin(), defvec.end(), 0) << " tag variants (due to undefined zvec or nvec)";

  if (bim_file_.size() > 0) {
    if (bim_file_.size() != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("bim_file_ is incompatible with num_snp_"));
    if (!extract.empty()) apply_extract(extract, bim_file_, &defvec);
    if (!exclude.empty()) apply_exclude(exclude, bim_file_, &defvec);
  }

  if (maf_threshold > 0) {
    if (mafvec_.size() != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("mafvec_ is incompatible with num_snp_"));
    apply_maf_threshold(maf_threshold, mafvec_, &defvec);
  }

  if (std::accumulate(defvec.begin(), defvec.end(), 0) == 0)
    BGMG_THROW_EXCEPTION(::std::runtime_error("no tag SNPs remain after applying specified filters"));

#pragma omp parallel
  {
    std::valarray<float> weight_random_pruning_local(0.0f, num_tag_);  // count the 1/w_ld weight of tag SNPs that has passed random pruning
    LdMatrixRow ld_matrix_row;

#pragma omp for schedule(dynamic)
    for (int prune_i = 0; prune_i < n; prune_i++) {
      std::mt19937_64 random_engine;
      random_engine.seed(seed_ + prune_i);

      std::vector<int> defvec_local(num_snp_, 0);
      std::vector<int> candidate_tag_indices(num_tag_, 0);
      std::vector<char> processed_tag_indices(num_tag_, 0);
      for (int i = 0; i < num_tag_; i++) candidate_tag_indices[i] = i;
      std::set<int> non_processed_tag_indices(candidate_tag_indices.begin(), candidate_tag_indices.end());

      // find tag SNPs that should be excluded
      for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
        if (!defvec[tag_to_snp()[tag_index]]) {
          processed_tag_indices[tag_index] = 1;         // mark as processed, and
          non_processed_tag_indices.erase(tag_index);   // remove from the set
        }
      }

      while (candidate_tag_indices.size() > 0) {
        // Here is the logic:
        // 1. select a random element X from the candidate_tag_indices
        // 2. if X is present in processed_tag_indices (collision case):
        //    - re-generate candidate_tag_indices from the set of non_processed_tag_indices
        //    - continue while loop.
        // 3. add X to passed_random_pruning
        // 4. query LD matrix for everything in LD with X (we asume that X will be part of that list). Then, for each Y in LD with X:
        //    - add Y to processed_tag_indices
        //    - remove Y from non_processed_tag_indices

        const int random_candidate_index = std::uniform_int_distribution<int>(0, candidate_tag_indices.size() - 1)(random_engine);
        const int random_tag_index = candidate_tag_indices[random_candidate_index];
        if (processed_tag_indices[random_tag_index]) {
          candidate_tag_indices.assign(non_processed_tag_indices.begin(), non_processed_tag_indices.end());
          // Validate that non_processed_tag_indices is consistent with processed_tag_indices.
          // for (int i = 0; i < num_tag_; i++) {
          //  const bool is_processed = (non_processed_tag_indices.find(i) == non_processed_tag_indices.end());
          //  if (processed_tag_indices[i] != is_processed) {
          //    LOG << " set_weights_randprune is stuck, processed_tag_indices inconsistent with non_processed_tag_indices. Cancel random pruning iteration " << prune_i;
          //    candidate_tag_indices.clear();
          //    break;
          //  }
          // }
          continue;
        }

        defvec_local[tag_to_snp()[random_tag_index]] = 1;
        ld_matrix_csr_.extract_tag_row(TagIndex(random_tag_index), &ld_matrix_row);
        auto iter_end = ld_matrix_row.end();
        int num_changes = 0;
        for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
          const int tag_index = snp_to_tag_[iter.index()];
          if (tag_index == -1) continue;
          const float r2_value = iter.r2();  // here we are interested in r2 (hvec is irrelevant)        
          if (r2_value < r2_threshold) continue;
          if (processed_tag_indices[tag_index]) continue;
          processed_tag_indices[tag_index] = 1;         // mark as processed, and
          non_processed_tag_indices.erase(tag_index);   // remove from the set
          num_changes++;
        }
        if (num_changes == 0) {
          LOG << " set_weights_randprune is stuck, num_changes=0. Cancel random pruning iteration " << prune_i;
          break;
        }
      }

      for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
        int snp_index = tag_to_snp()[tag_index];
        if (!defvec_local[snp_index]) continue;  // skip SNPs that did not pass random prunning

        if (use_w_ld) {  // compute 1/w_ld towards other tag SNPs passed random prunning
          float w_ld = 0.0f;
          ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), &ld_matrix_row);
          auto iter_end = ld_matrix_row.end();
          
          for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
            const int snp_index = iter.index();
            if (!defvec_local[snp_index]) continue;
            w_ld += iter.r2();
          }
          weight_random_pruning_local[tag_index] += 1.0f / w_ld;
        } else {
          weight_random_pruning_local[tag_index] += 1.0f;
        }
      }
    }

#pragma omp critical
    weight_random_pruning += weight_random_pruning_local;
  }

  weights_.clear(); weights_.resize(num_tag_, 0.0f);
  for (int i = 0; i < num_tag_; i++)
    weights_[i] = weight_random_pruning[i] / static_cast<float>(n);

  LOG << "<" << ss.str() << "), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::retrieve_zvec(int trait, int length, float* buffer) {
  if (length != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  std::vector<float>& zvec(*get_zvec(trait));
  if (zvec.size() != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec.size() != num_tag_"));
  LOG << " retrieve_zvec()";
  for (int i = 0; i < num_tag_; i++) buffer[i] = zvec[i];
  return 0;
}

int64_t BgmgCalculator::retrieve_nvec(int trait, int length, float* buffer) {
  if (length != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  std::vector<float>& nvec(*get_nvec(trait));
  if (nvec.size() != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("nvec.size() != num_tag_"));
  LOG << " retrieve_nvec()";
  for (int i = 0; i < num_tag_; i++) buffer[i] = nvec[i];
  return 0;
}

int64_t BgmgCalculator::retrieve_causalbetavec(int trait, int length, float* buffer) {
  if (length != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  const std::vector<float>& causalbetavec(*get_causalbetavec(trait));
  if (causalbetavec.size() != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("causalbetavec.size() != num_snp_"));
  LOG << " retrieve_causalbetavec()";
  for (int i = 0; i < num_snp_; i++) buffer[i] = causalbetavec[i];
  return 0;
}

int64_t BgmgCalculator::retrieve_mafvec(int length, float* buffer) {
  if (length != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (mafvec_.size() != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("mafvec_.size() != num_tag_"));
  LOG << " retrieve_mafvec()";
  for (int i = 0; i < num_snp_; i++) buffer[i] = mafvec_[i];
  return 0;
}

int64_t BgmgCalculator::retrieve_weights(int length, float* buffer) {
  if (length != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (weights_.size() != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("weights_.size() != num_tag_"));
  LOG << " retrieve_weights()";
  for (int i = 0; i < num_tag_; i++) buffer[i] = weights_[i];
  return 0;
}

int64_t BgmgCalculator::set_chrnumvec(int num_snp, const int* chrlabel) {
  if (!chrnumvec_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("can not set chrnumvec twice"));
  for (int i = 1; i < num_snp; i++) {
    if (chrlabel[i] < chrlabel[i - 1]) BGMG_THROW_EXCEPTION(::std::runtime_error("chrnumvec must be sorted"));
  }

  LOG << ">set_chrnumvec(" << num_snp << "); ";
  check_num_snp(num_snp);
  chrnumvec_.assign(chrlabel, chrlabel + num_snp);
  LOG << "<set_chrnumvec(" << num_snp << "); ";

  ld_matrix_csr_.init_chunks();
  return 0;
}

int64_t BgmgCalculator::retrieve_chrnumvec(int length, int* buffer) {
  if (length != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("wrong buffer size"));
  if (chrnumvec_.size() != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("chrnumvec_.size() != num_snp_"));
  LOG << " retrieve_chrnumvec()";
  for (int i = 0; i < num_snp_; i++) buffer[i] = chrnumvec_[i];
  return 0;
}

int64_t BgmgCalculator::num_ld_r_tag(int tag_index) {
  CHECK_TAG_INDEX((*this), tag_index);
  return ld_matrix_csr_.num_ld_r_tag(TagIndex(tag_index));
}

int64_t BgmgCalculator::retrieve_ld_r_tag(int tag_index, int length, int* snp_index, float* r) {
  if (length != num_ld_r_tag(tag_index)) BGMG_THROW_EXCEPTION(::std::runtime_error("length does not match num_ld_r_tag"));
  LOG << " retrieve_ld_r_tag(tag_index=" << tag_index << ")";
  
  LdMatrixRow ld_matrix_row;
  ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), &ld_matrix_row);
  auto iter_end = ld_matrix_row.end();
  int r_index = 0;
  for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
    snp_index[r_index] = iter.index();
    r[r_index] = iter.r();
    r_index++;
  }

  return 0;
}

int64_t BgmgCalculator::num_ld_r_chr(int chr_label) {
  int64_t retval = 0;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (chrnumvec_[tag_to_snp_[tag_index]] != chr_label) continue;
    retval += num_ld_r_tag(tag_index);
  }

  return retval;
}

int64_t BgmgCalculator::retrieve_ld_r_chr(int chr_label, int64_t length, int* tag_index_vector, int* snp_index_vector, float* r) {
  if (length != num_ld_r_chr(chr_label)) BGMG_THROW_EXCEPTION(::std::runtime_error("length does not match num_ld_r_chr"));
  LOG << " retrieve_ld_r_chr(chr_label=" << chr_label << ")";
  
  int r_index = 0;
  LdMatrixRow ld_matrix_row;
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (chrnumvec_[tag_to_snp_[tag_index]] != chr_label) continue;
    ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), &ld_matrix_row);
    auto iter_end = ld_matrix_row.end();
    for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
      tag_index_vector[r_index] = tag_index;
      snp_index_vector[r_index] = iter.index();
      r[r_index] = iter.r();
      r_index++;
    }
  }

  return 0;
}

int64_t BgmgCalculator::init(std::string bim_file, std::string frq_file, std::string chr_labels, std::string trait1_file, std::string trait2_file, std::string exclude, std::string extract, std::string exclude_ranges) {
  if (!trait2_file.empty() && trait1_file.empty())
    BGMG_THROW_EXCEPTION(std::runtime_error("trait2_file can be provided only together with trait1_file"));
  LOG << ">init(bim_file=" << bim_file << ", frq_file=" << frq_file << ", chr_labels=" << chr_labels << ", trait1_file=" << trait1_file << ", trait2_file=" << trait2_file << ", exclude=" << exclude << ", extract=" << extract << ", exclude_ranges=" << exclude_ranges << "); ";
  SimpleTimer timer(-1);

  std::vector<std::string> chr_labels_vector, bim_files, frq_files;

  if (chr_labels.empty()) {
    for (int i = 1; i <= 22; i++)
      chr_labels_vector.push_back(boost::lexical_cast<std::string>(i));
  } else {
    const std::string separators = " ,;\t\n\r";
    boost::trim_if(chr_labels, boost::is_any_of(separators));
    boost::split(chr_labels_vector, chr_labels, boost::is_any_of(separators), boost::token_compress_on);
  }

  if (bim_file.find("@") != std::string::npos) {
    for (auto chrlabel : chr_labels_vector) {
      bim_files.push_back(bim_file);
      boost::replace_all(bim_files.back(), "@", chrlabel);
    }
  }
  else {
    bim_files.push_back(bim_file);
  }

  if (frq_file.find("@") != std::string::npos) {
    for (auto chrlabel : chr_labels_vector) {
      frq_files.push_back(frq_file);
      boost::replace_all(frq_files.back(), "@", chrlabel);
    }
  }
  else if (!frq_file.empty()) {
    frq_files.push_back(frq_file);
  }

  for (auto& bim_file : bim_files) {
    if (!boost::filesystem::exists(bim_file)) {
      std::stringstream ss; ss << "ERROR: input file " << bim_file << " does not exist";
      BGMG_THROW_EXCEPTION(std::runtime_error(ss.str()));
    }
  }

  for (auto& frq_file : frq_files) {
    if (!boost::filesystem::exists(frq_file)) {
      std::stringstream ss; ss << "ERROR: input file " << frq_file << " does not exist";
      BGMG_THROW_EXCEPTION(std::runtime_error(ss.str()));
    }
  }

  if (!trait1_file.empty() && !boost::filesystem::exists(trait1_file))
    BGMG_THROW_EXCEPTION(std::runtime_error(trait1_file + " does not exist"));
  if (!trait2_file.empty() && !boost::filesystem::exists(trait2_file))
    BGMG_THROW_EXCEPTION(std::runtime_error(trait2_file + " does not exist"));

  bim_file_.clear();
  bim_file_.read(bim_files);
  bim_file_.find_snp_to_index_map();

  FrqFile frq_file_object;
  if (!frq_file.empty()) {
    frq_file_object.read(bim_file_, frq_files);
    frq_file_object.align_to_reference(bim_file_);  // raise an error if any of reference variants is not present in frq files.
  }

  std::vector<int> defvec(bim_file_.size(), 1);

  SumstatFile trait1_file_object;
  if (!trait1_file.empty()) {
    trait1_file_object.read(bim_file_, trait1_file, allow_ambiguous_snps_);
    for (int i = 0; i < bim_file_.size(); i++) {
      if (!std::isfinite(trait1_file_object.zscore()[i])) defvec[i] = 0;
      if (!std::isfinite(trait1_file_object.sample_size()[i])) defvec[i] = 0;
    }

    LOG << " constrain analysis to " << std::accumulate(defvec.begin(), defvec.end(), 0) << " tag variants (due to trait1_file='" << trait1_file << "')";
  }

  SumstatFile trait2_file_object;
  if (!trait2_file.empty()) {
    trait2_file_object.read(bim_file_, trait2_file, allow_ambiguous_snps_);
    for (int i = 0; i < bim_file_.size(); i++) {
      if (!std::isfinite(trait2_file_object.zscore()[i])) defvec[i] = 0;
      if (!std::isfinite(trait2_file_object.sample_size()[i])) defvec[i] = 0;
    }

    LOG << " constrain analysis to " << std::accumulate(defvec.begin(), defvec.end(), 0) << " tag variants (due to trait2_file='" << trait2_file << "')";
  }

  if (!extract.empty()) apply_extract(extract, bim_file_, &defvec);
  if (!exclude.empty()) apply_exclude(exclude, bim_file_, &defvec);
  if (!exclude_ranges.empty()) apply_exclude_ranges(exclude_ranges, bim_file_, &defvec);

  // Find tag indices
  std::vector<int> tag_indices;
  for (int i = 0; i < bim_file_.size(); i++)
    if (defvec[i] || use_complete_tag_indices_) tag_indices.push_back(i);

  // Initialize bgmg_calculator, e.i.
  // - set_tag_indices
  // - set_chrnumvec
  // - set_mafvec (if frq file is available)
  // - set_zvec, set_nvec (for each trait, if they are available)
  set_tag_indices(defvec.size(), tag_indices.size(), &tag_indices[0]);
  set_chrnumvec(bim_file_.size(), &bim_file_.chr_label()[0]);
  if (!frq_file.empty())
    set_mafvec(bim_file_.size(), &frq_file_object.frq()[0]);

  if (!trait1_file.empty()) {
    std::vector<float> zvec(tag_indices.size(), 0), nvec(tag_indices.size(), 0);
    for (int i = 0; i < tag_indices.size(); i++) {
      const int snp_index = tag_indices[i];
      zvec[i] = defvec[snp_index] ? trait1_file_object.zscore()[snp_index] : NAN;
      nvec[i] = defvec[snp_index] ? trait1_file_object.sample_size()[snp_index] : NAN;
    }
    set_zvec(1, tag_indices.size(), &zvec[0]);
    set_nvec(1, tag_indices.size(), &nvec[0]);
  }

  if (!trait2_file.empty()) {
    std::vector<float> zvec(tag_indices.size(), 0), nvec(tag_indices.size(), 0);
    for (int i = 0; i < tag_indices.size(); i++) {
      const int snp_index = tag_indices[i];
      zvec[i] = defvec[snp_index] ? trait2_file_object.zscore()[snp_index] : NAN;
      nvec[i] = defvec[snp_index] ? trait2_file_object.sample_size()[snp_index] : NAN;
    }
    set_zvec(2, tag_indices.size(), &zvec[0]);
    set_nvec(2, tag_indices.size(), &nvec[0]);
  }

  LOG << "<init(bim_file=" << bim_file << 
    ", frq_file=" << frq_file << 
    ", chr_labels=" << chr_labels << 
    ", trait1_file=" << trait1_file << 
    ", trait2_file=" << trait2_file << 
    ", exclude=" << exclude << 
    ", extract=" << extract <<
    ");  elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::set_weights_hardprune(float hardprune_subset, float hardprune_r2, float hardprune_maf, bool use_w_ld) {
  std::stringstream ss;
  ss << "set_weights_hardprune(subset=" << hardprune_subset << ", r2=" << hardprune_r2 << ", maf=" << hardprune_maf << ", use_w_ld=" << (use_w_ld ? "yes" : "no") << ", seed=" << seed_ << ")";
  LOG << ">" << ss.str();
  SimpleTimer timer(-1);

  std::vector<int> defvec(num_snp_, 0);
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if ((!zvec1_.empty()) && !std::isfinite(zvec1_[tag_index])) continue;
    if ((!nvec1_.empty()) && !std::isfinite(nvec1_[tag_index])) continue;
    if ((!zvec2_.empty()) && !std::isfinite(zvec2_[tag_index])) continue;
    if ((!nvec2_.empty()) && !std::isfinite(nvec2_[tag_index])) continue;
    defvec[tag_to_snp()[tag_index]] = 1;
  }
  LOG << " constrain analysis to " << std::accumulate(defvec.begin(), defvec.end(), 0) << " tag variants (due to undefined zvec or nvec)";

  if (hardprune_subset > 0) {
    if (hardprune_subset < 1) {
      hardprune_subset = int(hardprune_subset * float(num_snp_));
      LOG << " --hardprune-subset less than 1.0 is interpreted as a fraction, thus " << hardprune_subset << " SNPs will be selected";
    }

    // hardprune_subset acts all on SNPs from the reference, regardless of which ones are well defined
    std::vector<int> defvec_subset(num_snp_, 0);
    SubsetSampler subset_sampler((seed_ > 0) ? seed_ : (seed_ - 1), 1, num_snp_);
    uint32_t count = subset_sampler.sample_shuffle_count((uint32_t)hardprune_subset);
    for (int i = (num_snp_-count); i<num_snp_; i++) {
      int snp_index = subset_sampler.data()[i];
      defvec_subset[snp_index] = 1;
    }

    for (int i = 0; i < num_snp_; i++) {
      if (!defvec_subset[i]) {
        defvec[i] = 0;
      }
    }
    LOG << " constrain analysis to " << std::accumulate(defvec.begin(), defvec.end(), 0) << " tag variants (due to --hardprune-subset=" << hardprune_subset << ")";
  }

  if (hardprune_maf > 0) {
    if (mafvec_.size() != num_snp_) BGMG_THROW_EXCEPTION(::std::runtime_error("mafvec_ is incompatible with num_snp_"));
    apply_maf_threshold(hardprune_maf, mafvec_, &defvec);
  }

  if ((hardprune_r2 > 0) &&  (hardprune_r2 < 1)) {
    std::vector<float> buffer(num_tag_, NAN);
    SubsetSampler subset_sampler2((seed_ > 0) ? seed_ : (seed_ - 1), 1, num_tag_);
    subset_sampler2.sample_shuffle_count((uint32_t)num_tag_);
    for (int tag_index = 0; tag_index<num_tag_; tag_index++) {
      if (!defvec[tag_to_snp()[tag_index]]) continue;
      buffer[tag_index] = static_cast<float>(subset_sampler2.data()[tag_index]);
    }
    this->perform_ld_clump(hardprune_r2, num_tag_, &buffer[0]);
    for (int tag_index = 0; tag_index<num_tag_; tag_index++) {
      if (!std::isfinite(buffer[tag_index])) {
        defvec[tag_to_snp()[tag_index]] = 0;
      }
    }
    LOG << " constrain analysis to " << std::accumulate(defvec.begin(), defvec.end(), 0) << " tag variants (due to --hardprune-r2=" << hardprune_r2 << ")";
  }

  weights_.clear(); weights_.resize(num_tag_, 0.0f);

  if (use_w_ld) { // compute LD score towards other tag SNPs
    LdMatrixRow ld_matrix_row;
    for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
      if (!defvec[tag_to_snp()[tag_index]]) continue;
      ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), &ld_matrix_row);
      auto iter_end = ld_matrix_row.end();
      float w_ld = 0;
      for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
        const int snp_index = iter.index();
        if (!defvec[snp_index]) continue;
        w_ld += iter.r2();
      }
      weights_[tag_index] = 1.0f / w_ld;
    }
  } else {
    for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
      weights_[tag_index] = static_cast<float>(defvec[tag_to_snp()[tag_index]]);
    }
  }

  LOG << "<" << ss.str() << ", elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::read_trait_file(int trait_index, std::string trait_file, std::string exclude, std::string extract, std::string exclude_ranges) {
  if ((trait_index != 1) && (trait_index != 2)) BGMG_THROW_EXCEPTION(::std::runtime_error("trait must be 1 or 2"));
  if (trait_file.empty()) BGMG_THROW_EXCEPTION(std::runtime_error("trait_file can not be an empty string"));
  
  LOG << ">read_trait_file(trait_index=" << trait_index << ", trait_file=" << trait_file << ", exclude=" << exclude << ", extract=" << extract << ", exclude_ranges=" << exclude_ranges << "); ";
  SimpleTimer timer(-1);

  SumstatFile trait_file_object;
  trait_file_object.read(bim_file_, trait_file, allow_ambiguous_snps_);
  
  std::vector<int> defvec(bim_file_.size(), 1);
  if (!extract.empty()) apply_extract(extract, bim_file_, &defvec);
  if (!exclude.empty()) apply_exclude(exclude, bim_file_, &defvec);
  if (!exclude_ranges.empty()) apply_exclude_ranges(exclude_ranges, bim_file_, &defvec);

  std::vector<float> zvec(num_tag_, NAN), nvec(num_tag_, NAN);
  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    const int snp_index = tag_to_snp_[tag_index];
    const float zscore = trait_file_object.zscore()[snp_index];
    const float sample_size = trait_file_object.sample_size()[snp_index];
    if (std::isfinite(zscore) && std::isfinite(sample_size) && defvec[snp_index]) {
      zvec[tag_index] = zscore;
      nvec[tag_index] = sample_size;
    }
  }
  set_zvec(trait_index, num_tag_, &zvec[0]);
  set_nvec(trait_index, num_tag_, &nvec[0]);

  LOG << "<read_trait_file(trait_index=" << trait_index << ", trait_file=" << trait_file << ", exclude=" << exclude << ", extract=" << extract << ", exclude_ranges=" << exclude_ranges << ");  elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}

int64_t BgmgCalculator::convert_plink_ld(std::string plink_ld_gz, std::string plink_ld_bin) {
  PlinkLdFile plink_ld_file(bim_file_, plink_ld_gz);
  plink_ld_file.save_as_binary(plink_ld_bin);
  return 0;
}

int64_t BgmgCalculator::num_ld_r_tag_range(int tag_index_from, int tag_index_to) {
  return retrieve_ld_r_tag_range(tag_index_from, tag_index_to, -1, nullptr, nullptr, nullptr);
}

int64_t BgmgCalculator::retrieve_ld_r_tag_range(int tag_index_from, int tag_index_to, int length, int* tag_index_vector, int* snp_index_vector, float* r) {
  // length < 0 indicate that we just calculate num_ld_r_tag_range (i.e. the buffer size needed to copy out the result)
  CHECK_TAG_INDEX((*this), tag_index_from);
  CHECK_TAG_INDEX((*this), tag_index_to);
  int64_t num_r = 0;
  LdMatrixRow ld_matrix_row;
  for (int tag_index = tag_index_from; tag_index < tag_index_to; tag_index++) {
    ld_matrix_csr_.extract_tag_row(TagIndex(tag_index), &ld_matrix_row);
    auto iter_end = ld_matrix_row.end();
    for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
      const int snp_index = iter.index();
      num_r++;
      if (length < 0) continue;
      if (num_r > length) BGMG_THROW_EXCEPTION(::std::runtime_error("insufficient length for retrieve_ld_r_snp_range"));
      tag_index_vector[num_r - 1] = tag_index;
      snp_index_vector[num_r - 1] = snp_index;
      r[num_r - 1] = iter.r(); // here we are interested in r (hvec is irrelevant)
    }
  }
  LOG << ((length < 0) ? " num_ld_r_tag_range(" : " retrieve_ld_r_tag_range(from=") << tag_index_from << ", to=" << tag_index_to << "), return " << num_r;
  return (length < 0) ? num_r : 0;
}

int64_t BgmgCalculator::retrieve_fixed_effect_delta(int trait_index, int length, float* delta) {
  check_num_tag(length);
  std::valarray<float> fixed_effect_delta(0.0, num_tag_);
  calc_fixed_effect_delta_from_causalbetavec(trait_index, &fixed_effect_delta);
  for (int i = 0; i < num_tag_; i++) delta[i] = fixed_effect_delta[i];
}

void BgmgCalculator::calc_fixed_effect_delta_from_causalbetavec(int trait_index, std::valarray<float>* delta) {
  if (delta->size() != num_tag_) BGMG_THROW_EXCEPTION(::std::runtime_error("calc_fixed_effect_delta_from_causalbetavec expect delta to be already initialized"));
  *delta = 0.0f;

  const std::vector<float>& causalbetavec(*get_causalbetavec(trait_index));
  if (causalbetavec.empty()) return;
  check_num_snp(causalbetavec.size());

  LOG << ">calc_fixed_effect_delta_from_causalbetavec(trait_index=" << trait_index << ")";
  SimpleTimer timer(-1);

  std::vector<float> sqrt_hvec;
  find_hvec(*this, &sqrt_hvec);
  for (int i = 0; i < sqrt_hvec.size(); i++) sqrt_hvec[i] = sqrt(sqrt_hvec[i]);

#pragma omp parallel
  {
    LdMatrixRow ld_matrix_row;
    std::valarray<float> delta_local(0.0, num_tag_);

    // many entries in causalbetavec are expected to be zero, therefore static scheduler may give an unbalanced load
    // however it's fairly short operation anyway, so we don't bother too much.
#pragma omp for schedule(static)  // or dynamic?
    for (int snp_index = 0; snp_index < num_snp_; snp_index++) {
      if (causalbetavec[snp_index] == 0.0f) continue;
      ld_matrix_csr_.extract_snp_row(SnpIndex(snp_index), &ld_matrix_row);
      auto iter_end = ld_matrix_row.end();
      for (auto iter = ld_matrix_row.begin(); iter < iter_end; iter++) {
        const int tag_index = iter.index();
        const float r_value = iter.r();
        delta_local[tag_index] += r_value * sqrt_hvec[snp_index] * causalbetavec[snp_index];
      }
    }
#pragma omp critical
    (*delta) += delta_local;
  }

  const std::vector<float>& nvec(*get_nvec(trait_index));
  for (int i = 0; i < nvec.size(); i++) (*delta)[i] *= sqrt(nvec[i]);

  LOG << "<calc_fixed_effect_delta_from_causalbetavec(trait_index=" << trait_index  << "), elapsed time " << timer.elapsed_ms() << "ms";
}

void BgmgCalculator::find_z_minus_fixed_effect_delta(int trait_index, std::vector<float>* z_minus_fixed_effect_delta) {
  std::vector<float>& zvec(*get_zvec(trait_index));
  if (zvec.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("zvec is not set"));

  std::valarray<float> fixed_effect_delta(0.0, num_tag_);
  calc_fixed_effect_delta_from_causalbetavec(trait_index, &fixed_effect_delta);

  z_minus_fixed_effect_delta->resize(num_tag_);
  for (int i = 0; i < num_tag_; i++)
    z_minus_fixed_effect_delta->at(i) = std::isfinite(zvec[i]) ? (zvec[i] - fixed_effect_delta[i]) : zvec[i];
}

int BgmgCalculator::find_deftag_indices(const float* weights, std::vector<int>* deftag_indices) {
  if (weights == nullptr) {
    if (weights_.empty()) BGMG_THROW_EXCEPTION(::std::runtime_error("weights are not set"));
    weights = &weights_[0];
  }

  for (int tag_index = 0; tag_index < num_tag_; tag_index++) {
    if (weights[tag_index] == 0) continue;
    deftag_indices->push_back(tag_index);
  }
  return deftag_indices->size();
}

void BgmgCalculator::dump(BgmgSerializer* serializer) {
  serializer->dump_value(&num_snp_);
  serializer->dump_value(&num_tag_);
  serializer->dump_vector(&tag_to_snp_);
  serializer->dump_vector(&snp_to_tag_);
  serializer->dump_vector(&is_tag_);
  serializer->dump_vector(&chrnumvec_);

  if (serializer->is_load()) ld_matrix_csr_.init_chunks();
  ld_matrix_csr_.dump(serializer);

  serializer->dump_vector(&zvec1_);
  serializer->dump_vector(&nvec1_);
  serializer->dump_vector(&zvec2_);
  serializer->dump_vector(&nvec2_);
  serializer->dump_vector(&weights_);
  serializer->dump_vector(&mafvec_);

  serializer->dump_vector(&causalbetavec1_);
  serializer->dump_vector(&causalbetavec2_);

  serializer->dump_value(&k_max_);
  // serializer->dump_value(&k_max_pdf_);
  serializer->dump_value(&seed_);
  serializer->dump_value(&use_complete_tag_indices_);
  serializer->dump_value(&disable_snp_to_tag_map_);
  // serializer->dump_value(&allow_ambiguous_snps_);
  serializer->dump_value(&r2_min_);
  serializer->dump_value(&z1max_);
  serializer->dump_value(&z2max_);

  serializer->dump_value(&cost_calculator_);
  serializer->dump_value(&cubature_abs_error_);
  serializer->dump_value(&cubature_rel_error_);
  serializer->dump_value(&cubature_max_evals_);
  serializer->dump_value(&aux_option_);
  serializer->dump_value(&ld_format_version_);
  serializer->dump_value(&retrieve_ld_sum_type_);
  bim_file_.dump(serializer);
}

int64_t BgmgCalculator::save(const char* file) {
  LOG << ">save(filename=" << file << ")";  
  SimpleTimer timer(-1);
  BgmgSerializer serializer(file, SerializeDirection_Save);
  dump(&serializer);
  LOG << "<save(filename=" << file << "), elapsed time " << timer.elapsed_ms() << "ms";    
  return 0;
}

int64_t BgmgCalculator::load(const char* file) {
  LOG << ">load(filename=" << file << ")";
  SimpleTimer timer(-1);
  BgmgSerializer serializer(file, SerializeDirection_Load);
  dump(&serializer);
  LOG << "<load(filename=" << file << "), elapsed time " << timer.elapsed_ms() << "ms";
  return 0;
}
