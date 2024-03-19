#ifndef _BLAS_H_
#define _BLAS_H_


#ifndef _BLAS_UNIT_LEN_
#define _BLAS_UNIT_LEN_ 4
#endif

#ifndef _BLAS_U32_H_
#define _BLAS_U32_H_

#ifndef _GF16_H_
#define _GF16_H_

#include <stdint.h>

// gf4 := gf2[x]/x^2+x+1
static inline uint8_t gf4_mul_2(uint8_t a) {
    uint8_t r = a << 1;
    r ^= (a >> 1) * 7;
    return r;
}

static inline uint8_t gf4_mul_3(uint8_t a) {
    uint8_t msk = (a - 2) >> 1;
    return (msk & (a * 3)) | ((~msk) & (a - 1));
}

static inline uint8_t gf4_mul(uint8_t a, uint8_t b) {
    uint8_t r = a * (b & 1);
    return r ^ (gf4_mul_2(a) * (b >> 1));
}

static inline uint8_t gf4_squ(uint8_t a) {
    return a ^ (a >> 1);
}

static inline uint8_t gf4_inv(uint8_t a) {
    return a ^ (a >> 1);
}

static inline uint32_t gf4v_mul_2_u32(uint32_t a) {
    uint32_t bit0 = a & 0x55555555;
    uint32_t bit1 = a & 0xaaaaaaaa;
    return (bit0 << 1) ^ bit1 ^ (bit1 >> 1);
}

static inline uint32_t gf4v_mul_3_u32(uint32_t a) {
    uint32_t bit0 = a & 0x55555555;
    uint32_t bit1 = a & 0xaaaaaaaa;
    return (bit0 << 1) ^ bit0 ^ (bit1 >> 1);
}

static inline uint32_t gf4v_mul_u32(uint32_t a, uint8_t b) {
    uint32_t bit0_b = ((uint32_t) 0) - ((uint32_t)(b & 1));
    uint32_t bit1_b = ((uint32_t) 0) - ((uint32_t)((b >> 1) & 1));
    return (a & bit0_b) ^ (bit1_b & gf4v_mul_2_u32(a));
}

static inline uint32_t _gf4v_mul_u32_u32(uint32_t a0, uint32_t a1, uint32_t b0, uint32_t b1) {
    uint32_t c0 = a0 & b0;
    uint32_t c2 = a1 & b1;
    uint32_t c1_ = (a0 ^ a1) & (b0 ^ b1);
    return ((c1_ ^ c0) << 1) ^ c0 ^ c2;
}

static inline uint32_t gf4v_mul_u32_u32(uint32_t a, uint32_t b) {
    uint32_t a0 = a & 0x55555555;
    uint32_t a1 = (a >> 1) & 0x55555555;
    uint32_t b0 = b & 0x55555555;
    uint32_t b1 = (b >> 1) & 0x55555555;

    return _gf4v_mul_u32_u32(a0, a1, b0, b1);
}

static inline uint32_t gf4v_squ_u32(uint32_t a) {
    uint32_t bit1 = a & 0xaaaaaaaa;
    return a ^ (bit1 >> 1);
}

static inline uint8_t gf16_is_nonzero(uint8_t a) {
    unsigned a4 = a & 0xf;
    unsigned r = ((unsigned) 0) - a4;
    r >>= 4;
    return r & 1;
}

// gf16 := gf4[y]/y^2+y+x
static inline uint8_t gf16_mul(uint8_t a, uint8_t b) {
    uint8_t a0 = a & 3;
    uint8_t a1 = (a >> 2);
    uint8_t b0 = b & 3;
    uint8_t b1 = (b >> 2);
    uint8_t a0b0 = gf4_mul(a0, b0);
    uint8_t a1b1 = gf4_mul(a1, b1);
    uint8_t a0b1_a1b0 = gf4_mul(a0 ^ a1, b0 ^ b1) ^ a0b0 ^ a1b1;
    uint8_t a1b1_x2 = gf4_mul_2(a1b1);
    return ((a0b1_a1b0 ^ a1b1) << 2) ^ a0b0 ^ a1b1_x2;
}

static inline uint8_t gf16_squ(uint8_t a) {
    uint8_t a0 = a & 3;
    uint8_t a1 = (a >> 2);
    a1 = gf4_squ(a1);
    uint8_t a1squ_x2 = gf4_mul_2(a1);
    return (a1 << 2) ^ a1squ_x2 ^ gf4_squ(a0);
}

static inline uint8_t gf16_inv(uint8_t a) {
    uint8_t a2 = gf16_squ(a);
    uint8_t a4 = gf16_squ(a2);
    uint8_t a8 = gf16_squ(a4);
    uint8_t a6 = gf16_mul(a4, a2);
    return gf16_mul(a8, a6);
}

static inline uint8_t gf16_mul_4(uint8_t a) {
    return (((a << 2) ^ a) & (8 + 4)) ^ gf4_mul_2(a >> 2);
}

static inline uint8_t gf16_mul_8(uint8_t a) {
    uint8_t a0 = a & 3;
    uint8_t a1 = a >> 2;
    return (gf4_mul_2(a0 ^ a1) << 2) | gf4_mul_3(a1);
}

// gf16 := gf4[y]/y^2+y+x
static inline uint32_t gf16v_mul_u32(uint32_t a, uint8_t b) {
    uint32_t axb0 = gf4v_mul_u32(a, b);
    uint32_t axb1 = gf4v_mul_u32(a, b >> 2);
    uint32_t a0b1 = (axb1 << 2) & 0xcccccccc;
    uint32_t a1b1 = axb1 & 0xcccccccc;
    uint32_t a1b1_2 = a1b1 >> 2;

    return axb0 ^ a0b1 ^ a1b1 ^ gf4v_mul_2_u32(a1b1_2);
}

static inline uint32_t _gf16v_mul_u32_u32(uint32_t a0, uint32_t a1, uint32_t a2, uint32_t a3, uint32_t b0, uint32_t b1, uint32_t b2, uint32_t b3) {
    uint32_t c0 = _gf4v_mul_u32_u32(a0, a1, b0, b1);
    uint32_t c1_ = _gf4v_mul_u32_u32(a0 ^ a2, a1 ^ a3, b0 ^ b2, b1 ^ b3);

    uint32_t c2_0 = a2 & b2;
    uint32_t c2_2 = a3 & b3;
    uint32_t c2_1_ = (a2 ^ a3) & (b2 ^ b3);
    uint32_t c2_r0 = c2_0 ^ c2_2;
    uint32_t c2_r1 = c2_0 ^ c2_1_;
    //uint32_t c2 = c2_r0^(c2_r1<<1);
    // GF(4) x2: (bit0<<1)^bit1^(bit1>>1);
    return ((c1_ ^ c0) << 2) ^ c0 ^ (c2_r0 << 1) ^ c2_r1 ^ (c2_r1 << 1);
}

static inline uint32_t gf16v_mul_u32_u32(uint32_t a, uint32_t b) {
    uint32_t a0 = a & 0x11111111;
    uint32_t a1 = (a >> 1) & 0x11111111;
    uint32_t a2 = (a >> 2) & 0x11111111;
    uint32_t a3 = (a >> 3) & 0x11111111;
    uint32_t b0 = b & 0x11111111;
    uint32_t b1 = (b >> 1) & 0x11111111;
    uint32_t b2 = (b >> 2) & 0x11111111;
    uint32_t b3 = (b >> 3) & 0x11111111;

    return _gf16v_mul_u32_u32(a0, a1, a2, a3, b0, b1, b2, b3);
}

static inline uint8_t gf256v_reduce_u32(uint32_t a) {
    uint16_t *aa = (uint16_t *) (&a);
    uint16_t r = aa[0] ^ aa[1];
    uint8_t *rr = (uint8_t *) (&r);
    return rr[0] ^ rr[1];
}

static inline uint8_t gf16v_reduce_u32(uint32_t a) {
    uint8_t r256 = gf256v_reduce_u32(a);
    return (r256 & 0xf) ^ (r256 >> 4);
}

static inline uint32_t gf16v_squ_u32(uint32_t a) {
    uint32_t a2 = gf4v_squ_u32(a);

    return a2 ^ gf4v_mul_2_u32((a2 >> 2) & 0x33333333);
}

static inline uint32_t gf16v_mul_4_u32(uint32_t a) {
    uint32_t a1 = a & 0xcccccccc;
    uint32_t a0 = (a << 2) & 0xcccccccc;
    return a0 ^ a1 ^ gf4v_mul_2_u32(a1 >> 2);
}

static inline uint32_t gf16v_mul_8_u32(uint32_t a) {
    uint32_t a1 = a & 0xcccccccc;
    uint32_t a0 = (a << 2) & 0xcccccccc;
    return gf4v_mul_2_u32(a0 ^ a1) | gf4v_mul_3_u32(a1 >> 2);
}

static inline uint8_t gf256_is_nonzero(uint8_t a) {
    unsigned a8 = a;
    unsigned r = ((unsigned) 0) - a8;
    r >>= 8;
    return r & 1;
}

// gf256 := gf16[X]/X^2+X+xy
static inline uint8_t gf256_mul(uint8_t a, uint8_t b) {
    uint8_t a0 = a & 15;
    uint8_t a1 = (a >> 4);
    uint8_t b0 = b & 15;
    uint8_t b1 = (b >> 4);
    uint8_t a0b0 = gf16_mul(a0, b0);
    uint8_t a1b1 = gf16_mul(a1, b1);
    uint8_t a0b1_a1b0 = gf16_mul(a0 ^ a1, b0 ^ b1) ^ a0b0 ^ a1b1;
    uint8_t a1b1_x8 = gf16_mul_8(a1b1);
    return ((a0b1_a1b0 ^ a1b1) << 4) ^ a0b0 ^ a1b1_x8;
}

static inline uint8_t gf256_mul_gf16(uint8_t a, uint8_t gf16_b) {
    uint8_t a0 = a & 15;
    uint8_t a1 = (a >> 4);
    uint8_t b0 = gf16_b & 15;
    uint8_t a0b0 = gf16_mul(a0, b0);
    uint8_t a1b0 = gf16_mul(a1, b0);
    return a0b0 ^ (a1b0 << 4);
}

static inline uint8_t gf256_squ(uint8_t a) {
    uint8_t a0 = a & 15;
    uint8_t a1 = (a >> 4);
    a1 = gf16_squ(a1);
    uint8_t a1squ_x8 = gf16_mul_8(a1);
    return (a1 << 4) ^ a1squ_x8 ^ gf16_squ(a0);
}

static inline uint8_t gf256_inv(uint8_t a) {
    // 128+64+32+16+8+4+2 = 254
    uint8_t a2 = gf256_squ(a);
    uint8_t a4 = gf256_squ(a2);
    uint8_t a8 = gf256_squ(a4);
    uint8_t a4_2 = gf256_mul(a4, a2);
    uint8_t a8_4_2 = gf256_mul(a4_2, a8);
    uint8_t a64_ = gf256_squ(a8_4_2);
    a64_ = gf256_squ(a64_);
    a64_ = gf256_squ(a64_);
    uint8_t a64_2 = gf256_mul(a64_, a8_4_2);
    uint8_t a128_ = gf256_squ(a64_2);
    return gf256_mul(a2, a128_);
}

static inline uint32_t gf256v_mul_u32(uint32_t a, uint8_t b) {
    uint32_t axb0 = gf16v_mul_u32(a, b);
    uint32_t axb1 = gf16v_mul_u32(a, b >> 4);
    uint32_t a0b1 = (axb1 << 4) & 0xf0f0f0f0;
    uint32_t a1b1 = axb1 & 0xf0f0f0f0;
    uint32_t a1b1_4 = a1b1 >> 4;

    return axb0 ^ a0b1 ^ a1b1 ^ gf16v_mul_8_u32(a1b1_4);
}

static inline uint32_t gf256v_squ_u32(uint32_t a) {
    uint32_t a2 = gf16v_squ_u32(a);
    uint32_t ar = (a2 >> 4) & 0x0f0f0f0f;

    return a2 ^ gf16v_mul_8_u32(ar);
}

static inline uint32_t gf256v_mul_gf16_u32(uint32_t a, uint8_t gf16_b) {
    return gf16v_mul_u32(a, gf16_b);
}


// gf256 := gf16[X]/X^2+X+xy
static inline uint32_t gf256v_mul_0x10_u32(uint32_t a) {
    uint32_t a0 = a&0x0f0f0f0f;
    uint32_t a1 = a&0xf0f0f0f0;
    uint32_t a1x = gf16v_mul_8_u32(a1>>4);
    return (a0<<4)^a1^a1x;
}

static inline uint32_t gf256v_mul_0x20_u32(uint32_t a) {
    uint32_t a0 = gf4v_mul_2_u32(a&0x0f0f0f0f);
    uint32_t a1 = gf4v_mul_2_u32(a&0xf0f0f0f0);
    uint32_t a1x = gf16v_mul_8_u32(a1>>4);
    return (a0<<4)^a1^a1x;
}

static inline uint32_t gf256v_mul_0x40_u32(uint32_t a) {
    uint32_t a0 = gf16v_mul_4_u32(a&0x0f0f0f0f);
    uint32_t a1 = gf16v_mul_4_u32(a&0xf0f0f0f0);
    uint32_t a1x = gf16v_mul_8_u32(a1>>4);
    return (a0<<4)^a1^a1x;
}

static inline uint32_t gf256v_mul_0x80_u32(uint32_t a) {
    uint32_t a0 = gf16v_mul_8_u32(a&0x0f0f0f0f);
    uint32_t a1 = gf16v_mul_8_u32(a&0xf0f0f0f0);
    uint32_t a1x = gf16v_mul_8_u32(a1>>4);
    return (a0<<4)^a1^a1x;
}

#endif // _GF16_H_

#include <string.h>
#include <stdint.h>

static inline void _gf256v_add_u32(uint8_t *accu_b, const uint8_t *a, unsigned _num_byte) {
    unsigned n_u32 = _num_byte >> 2;
    for (unsigned i = 0; i < n_u32; i++) {
      uint32_t bx;
      uint32_t ax;
      memcpy(&bx,accu_b+4*i,4);
      memcpy(&ax,a+4*i,4);
      bx ^= ax;
      memcpy(accu_b+4*i,&bx,4);
    }

    unsigned rem = _num_byte & 3;
    if( !rem ) return;
    a += (n_u32 << 2);
    accu_b += (n_u32 << 2);
    for (unsigned i = 0; i < rem; i++) accu_b[i] ^= a[i];
}

static inline void _gf256v_conditional_add_u32(uint8_t *accu_b, uint8_t condition, const uint8_t *a, unsigned _num_byte) {
    uint32_t pr_u32 = ((uint32_t) 0) - ((uint32_t) condition);
    uint8_t pr_u8 = pr_u32 & 0xff;

    unsigned n_u32 = _num_byte >> 2;
    for (unsigned i = 0; i < n_u32; i++) {
      uint32_t bx;
      uint32_t ax;
      memcpy(&bx,accu_b+4*i,4);
      memcpy(&ax,a+4*i,4);
      bx ^= ax&pr_u32;
      memcpy(accu_b+4*i,&bx,4);
    }

    unsigned rem = _num_byte & 3;
    if( !rem ) return;
    a += (n_u32 << 2);
    accu_b += (n_u32 << 2);
    for (unsigned i = 0; i < rem; i++) accu_b[i] ^= (a[i] & pr_u8);
}

static inline void _gf16v_mul_scalar_u32(uint8_t *a, uint8_t gf16_b, unsigned _num_byte) {
    unsigned n_u32 = _num_byte >> 2;
    for (unsigned i = 0; i < n_u32; i++) {
      uint32_t ax;
      memcpy(&ax,a+4*i,4);
      ax = gf16v_mul_u32(ax, gf16_b);
      memcpy(a+4*i,&ax,4);
    }

    unsigned rem = _num_byte & 3;
    if( !rem ) return;
    union tmp_32 {
        uint8_t u8[4];
        uint32_t u32;
    } t;
    a += (n_u32 << 2);
    for (unsigned i = 0; i < rem; i++) t.u8[i] = a[i];
    t.u32 = gf16v_mul_u32(t.u32, gf16_b);
    for (unsigned i = 0; i < rem; i++) a[i] = t.u8[i];
}

static inline void _gf256v_mul_scalar_u32(uint8_t *a, uint8_t b, unsigned _num_byte) {
    unsigned n_u32 = _num_byte >> 2;
    for (unsigned i = 0; i < n_u32; i++) {
      uint32_t ax;
      memcpy(&ax,a+4*i,4);
      ax = gf256v_mul_u32(ax, b);
      memcpy(a+4*i,&ax,4);
    }

    unsigned rem = _num_byte & 3;
    if( !rem ) return;
    union tmp_32 {
        uint8_t u8[4];
        uint32_t u32;
    } t;
    a += (n_u32 << 2);
    for (unsigned i = 0; i < rem; i++) t.u8[i] = a[i];
    t.u32 = gf256v_mul_u32(t.u32, b);
    for (unsigned i = 0; i < rem; i++) a[i] = t.u8[i];
}

static inline void _gf16v_madd_u32(uint8_t *accu_c, const uint8_t *a, uint8_t gf16_b, unsigned _num_byte) {
    unsigned n_u32 = _num_byte >> 2;
    for (unsigned i = 0; i < n_u32; i++) {
      uint32_t ax;
      uint32_t cx;
      memcpy(&ax,a+4*i,4);
      memcpy(&cx,accu_c+4*i,4);
      cx ^= gf16v_mul_u32(ax, gf16_b);
      memcpy(accu_c+4*i,&cx,4);
    }

    unsigned rem = _num_byte & 3;
    if( !rem ) return;
    union tmp_32 {
        uint8_t u8[4];
        uint32_t u32;
    } t;
    accu_c += (n_u32 << 2);
    a += (n_u32 << 2);
    for (unsigned i = 0; i < rem; i++) t.u8[i] = a[i];
    t.u32 = gf16v_mul_u32(t.u32, gf16_b);
    for (unsigned i = 0; i < rem; i++) accu_c[i] ^= t.u8[i];
}

static inline void _gf256v_madd_u32(uint8_t *accu_c, const uint8_t *a, uint8_t gf256_b, unsigned _num_byte) {
    unsigned n_u32 = _num_byte >> 2;
    for (unsigned i = 0; i < n_u32; i++) {
      uint32_t ax;
      uint32_t cx;
      memcpy(&ax,a+4*i,4);
      memcpy(&cx,accu_c+4*i,4);
      cx ^= gf256v_mul_u32(ax, gf256_b);
      memcpy(accu_c+4*i,&cx,4);
    }

    unsigned rem = _num_byte & 3;
    if( !rem ) return;
    union tmp_32 {
        uint8_t u8[4];
        uint32_t u32;
    } t;
    accu_c += (n_u32 << 2);
    a += (n_u32 << 2);
    for (unsigned i = 0; i < rem; i++) t.u8[i] = a[i];
    t.u32 = gf256v_mul_u32(t.u32, gf256_b);
    for (unsigned i = 0; i < rem; i++) accu_c[i] ^= t.u8[i];
}


#endif // _BLAS_U32_H_

#define gf16v_mul_scalar _gf16v_mul_scalar_u32
#define gf16v_madd _gf16v_madd_u32

#define gf256v_add _gf256v_add_u32
#define gf256v_mul_scalar _gf256v_mul_scalar_u32
#define gf256v_madd _gf256v_madd_u32

#define gf256v_conditional_add _gf256v_conditional_add_u32


#ifndef _BLAS_MATRIX_H_
#define _BLAS_MATRIX_H_

#include <stdint.h>

#ifdef  __cplusplus
extern  "C" {
#endif

void gf16mat_prod(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);

void gf256mat_prod(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);

unsigned gf16mat_solve_linear_eq_32x32(uint8_t *sol, const uint8_t *inp_mat, const uint8_t *c_terms );

unsigned gf256mat_solve_linear_eq_48x48(uint8_t *sol, const uint8_t *inp_mat, const uint8_t *c_terms );

unsigned gf256mat_solve_linear_eq_64x64(uint8_t *sol, const uint8_t *inp_mat, const uint8_t *c_terms );

unsigned gf16mat_inv_32x32(uint8_t *inv_a, const uint8_t *a );

unsigned gf256mat_inv_32x32(uint8_t *inv_a, const uint8_t *a );

unsigned gf256mat_inv_36x36(uint8_t *inv_a, const uint8_t *a );

static inline uint8_t gf256v_get_ele(const uint8_t *a, unsigned i) { return a[i]; }

#ifdef  __cplusplus
}
#endif
#endif  // _BLAS_MATRIX_H
#endif // _BLAS_H_
#include <stdint.h>


static inline
unsigned idx_of_trimat( unsigned i_row , unsigned j_col , unsigned dim )
{
    return (dim + dim - i_row + 1 )*i_row/2 + j_col - i_row;
}

static inline uint8_t gf16v_get_ele(const uint8_t *a, unsigned i) {
    uint8_t r = a[i >> 1];
    uint8_t r0 = r&0xf;
    uint8_t r1 = r>>4;
    uint8_t m = (uint8_t)(-(i&1));
    return (r1&m)|((~m)&r0);
}

static inline
unsigned idx_of_2trimat( unsigned i_row , unsigned j_col , unsigned n_var )
{
   if( i_row > j_col ) return idx_of_trimat(j_col,i_row,n_var);
   else return idx_of_trimat(i_row,j_col,n_var);
}

void gf256v_set_zero(uint8_t *b, unsigned _num_byte) {
    gf256v_add(b, b, _num_byte);
}

unsigned gf256v_is_zero(const uint8_t *a, unsigned _num_byte) {
    uint8_t r = 0;
    while( _num_byte-- ) { r |= a[0]; a++; }
    return (0 == r);
}

void UpperTrianglize( unsigned char * btriC , const unsigned char * bA , unsigned Awidth, unsigned size_batch )
{
    unsigned char * runningC = btriC;
    unsigned Aheight = Awidth;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<i;j++) {
            unsigned idx = idx_of_trimat(j,i,Aheight);
            gf256v_add( btriC + idx*size_batch , bA + size_batch*(i*Awidth+j) , size_batch );
        }
        gf256v_add( runningC , bA + size_batch*(i*Awidth+i) , size_batch*(Aheight-i) );
        runningC += size_batch*(Aheight-i);
    }
}

void batch_trimat_madd_gf16( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Awidth = Bheight;
    unsigned Aheight = Awidth;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                if(k<i) continue;
                gf16v_madd( bC , & btriA[ (k-i)*size_batch ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
        btriA += (Aheight-i)*size_batch;
    }
}

void batch_trimat_madd_gf256( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Awidth = Bheight;
    unsigned Aheight = Awidth;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                if(k<i) continue;
                gf256v_madd( bC , & btriA[ (k-i)*size_batch ] , gf256v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
        btriA += (Aheight-i)*size_batch;
    }
}

void batch_trimatTr_madd_gf16( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Aheight = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                if(i<k) continue;
                gf16v_madd( bC , & btriA[ size_batch*(idx_of_trimat(k,i,Aheight)) ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}

void batch_trimatTr_madd_gf256( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Aheight = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                if(i<k) continue;
                gf256v_madd( bC , & btriA[ size_batch*(idx_of_trimat(k,i,Aheight)) ] , gf256v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}

void batch_2trimat_madd_gf16( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Aheight = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                if(i==k) continue;
                gf16v_madd( bC , & btriA[ size_batch*(idx_of_2trimat(i,k,Aheight)) ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}

void batch_2trimat_madd_gf256( unsigned char * bC , const unsigned char* btriA ,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Aheight = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                if(i==k) continue;
                gf256v_madd( bC , & btriA[ size_batch*(idx_of_2trimat(i,k,Aheight)) ] , gf256v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}

void batch_matTr_madd_gf16( unsigned char * bC , const unsigned char* A_to_tr , unsigned Aheight, unsigned size_Acolvec, unsigned Awidth,
        const unsigned char* bB, unsigned Bwidth, unsigned size_batch )
{
    unsigned Atr_height = Awidth;
    unsigned Atr_width  = Aheight;
    for(unsigned i=0;i<Atr_height;i++) {
        for(unsigned j=0;j<Atr_width;j++) {
            gf16v_madd( bC , & bB[ j*Bwidth*size_batch ] , gf16v_get_ele( &A_to_tr[size_Acolvec*i] , j ) , size_batch*Bwidth );
        }
        bC += size_batch*Bwidth;
    }
}

void batch_matTr_madd_gf256( unsigned char * bC , const unsigned char* A_to_tr , unsigned Aheight, unsigned size_Acolvec, unsigned Awidth,
        const unsigned char* bB, unsigned Bwidth, unsigned size_batch )
{
    unsigned Atr_height = Awidth;
    unsigned Atr_width  = Aheight;
    for(unsigned i=0;i<Atr_height;i++) {
        for(unsigned j=0;j<Atr_width;j++) {
            gf256v_madd( bC , & bB[ j*Bwidth*size_batch ] , gf256v_get_ele( &A_to_tr[size_Acolvec*i] , j ) , size_batch*Bwidth );
        }
        bC += size_batch*Bwidth;
    }
}

void batch_bmatTr_madd_gf16( unsigned char *bC , const unsigned char *bA_to_tr, unsigned Awidth_before_tr,
        const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch )
{
    const unsigned char *bA = bA_to_tr;
    unsigned Aheight = Awidth_before_tr;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                gf16v_madd( bC , & bA[ size_batch*(i+k*Aheight) ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}

void batch_bmatTr_madd_gf256( unsigned char *bC , const unsigned char *bA_to_tr, unsigned Awidth_before_tr,
        const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch )
{
    const unsigned char *bA = bA_to_tr;
    unsigned Aheight = Awidth_before_tr;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                gf256v_madd( bC , & bA[ size_batch*(i+k*Aheight) ] , gf256v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
    }
}

void batch_mat_madd_gf16( unsigned char * bC , const unsigned char* bA , unsigned Aheight,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Awidth = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                gf16v_madd( bC , & bA[ k*size_batch ] , gf16v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
        bA += (Awidth)*size_batch;
    }
}

void batch_mat_madd_gf256( unsigned char * bC , const unsigned char* bA , unsigned Aheight,
        const unsigned char* B , unsigned Bheight, unsigned size_Bcolvec , unsigned Bwidth, unsigned size_batch )
{
    unsigned Awidth = Bheight;
    for(unsigned i=0;i<Aheight;i++) {
        for(unsigned j=0;j<Bwidth;j++) {
            for(unsigned k=0;k<Bheight;k++) {
                gf256v_madd( bC , & bA[ k*size_batch ] , gf256v_get_ele( &B[j*size_Bcolvec] , k ) , size_batch );
            }
            bC += size_batch;
        }
        bA += (Awidth)*size_batch;
    }
}