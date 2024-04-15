static void cmov(unsigned char *r, const unsigned char *x, size_t len, unsigned char b)
{
  size_t i;
  b = -b;
  for(i=0;i<len;i++)
    r[i] ^= b & (x[i] ^ r[i]);
}