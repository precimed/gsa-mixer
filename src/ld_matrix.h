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

#include <string>
#include <valarray>

#include "ld_matrix_csr.h"

void generate_ld_matrix_from_bed_file(std::string bfile, float r2min, float ldscore_r2min, int ld_window, float ld_window_kb, std::string out_file);

void save_ld_matrix(LdMatrixCsrChunk& chunk,
                    std::vector<float>& freqvec,
                    std::vector<float>& ld_r2_sum,
                    std::vector<float>& ld_r2_sum_adjust_for_hvec,
                    std::string filename);

void load_ld_matrix(std::string filename,
                    LdMatrixCsrChunk* chunk,
                    std::vector<float>* freqvec,
                    std::vector<float>* ld_r2_sum,
                    std::vector<float>* ld_r2_sum_adjust_for_hvec);

void load_ld_matrix_version0(std::string filename,
                             std::vector<int>* snp_index,
                             std::vector<int>* snp_other_index,
                             std::vector<float>* r2);