static void cswap4x(gfe *x, gfe *y, int b)
{
  crypto_uint32 db = -b;
  crypto_int32 t;
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<5;j++) {
      t = x[i].v[j] ^ y[i].v[j];
      t &= db;
      x[i].v[j] ^= t;
      y[i].v[j] ^= t;
    }
}