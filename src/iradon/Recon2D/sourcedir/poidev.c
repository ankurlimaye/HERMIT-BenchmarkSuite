#include <math.h>

extern float ran1(int *);
extern float gammln(float);

float poidev(float xm, int *idum) {
  static float sq, alxm, g, oldm = (-1.0);
  float em, t, y;

  if (xm < 12.0) {
    if (xm != oldm) {
      oldm = xm;
      g = exp(-xm);
    }
    em = -1;
    t = 1.0;
    do {
      em += 1.0;
      t *= ran1(idum);
    } while (t > g);
  } else {
    if (xm != oldm) {
      oldm = xm;
      sq = sqrt(2.0 * xm);
      alxm = log(xm);
      g = xm * alxm - gammln(xm + 1.0);
    }
    do {
      do {
        y = tan(PI * ran1(idum));
        em = sq * y + xm;
      } while (em < 0.0);
      em = floor(em);
      t = 0.9 * (1.0 + y * y) * exp(em * alxm - gammln(em + 1.0) - g);
    } while (ran1(idum) > t);
  }
  return em;
}

