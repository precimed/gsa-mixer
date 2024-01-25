#include <stdint.h>
#include <stddef.h>

#if defined(_MSC_VER)
#define __restrict__ __restrict
#endif

/***
* These functions compute fast successive differences, and recover the original
* values from the fast successive differences (i.e., they compute a prefix sum)
* using fast SIMD instructions.
*
* Reference :
* Daniel Lemire, Nathan Kurz, Leonid Boytsov, SIMD Compression and the Intersection of Sorted
* Integers, Software: Practice and Experience (to appear)
* http://arxiv.org/abs/1401.6399
*/

#ifdef __cplusplus
extern "C" {
#endif

// write to output the successive differences of input (input[0]-starting_point, input[1]-input[2], ...)
// there are "length" values in input and output
// input and output must be distinct
// it can make sense to set to zero by default
void compute_deltas(const uint32_t * __restrict__ input, size_t length, uint32_t * __restrict__ output, uint32_t starting_point);



// write to buffer the successive differences of buffer (buffer[0]-starting_point, buffer[1]-buffer[2], ...)
// there are "length" values in buffer
// it can make sense to set to zero by default
void compute_deltas_inplace(uint32_t * buffer, size_t length, uint32_t starting_point);


// write to output the successive differences of input (input[0]-starting_point, input[1]-input[2], ...)
// there are "length" values in input and output
// input and output must be distinct
void compute_prefix_sum(const uint32_t * __restrict__ input, size_t length, uint32_t * __restrict__ output, uint32_t starting_point);


// write to buffer the successive differences of buffer (buffer[0]-starting_point, buffer[1]-buffer[2], ...)
// there are "length" values in buffer
void compute_prefix_sum_inplace(uint32_t * buffer, size_t length, uint32_t starting_point);

#ifdef __cplusplus
}
#endif